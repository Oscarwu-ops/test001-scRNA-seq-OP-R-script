###############################################################################
#  scRNA_OCP_full_pipeline.R  --  CI-friendly, auto-extract, OCP/OC DEG ready
#
#  用法：
#    Rscript scRNA_OCP_full_pipeline.R --data_dir . --out_dir docs
#
#  環境變數（也可用 CLI 取代）：
#    DATA_DIR    = "."      # 原始檔案根目錄；會自動遞迴解壓並尋找 *_matrix*.mtx.gz
#    OUTPUT_DIR  = "docs"   # 所有輸出放這裡（圖、CSV、RDS、log、summary）
#    QUICK       = "1"      # "1"=抽樣≤5000細胞做快速驗證；"0"=全量
#    STRICT_FAIL = "0"      # "1"=遇到關鍵錯誤直接 exit 1；"0"=寫 _error.txt 後正常結束
#    DO_ENRICH   = "0"      # "1"=做 GO 富集（需 Bioconductor；建議正式產物時才開）
###############################################################################

## ========== 0. 前置：安裝必要 CRAN、讀參數、工具 ==========
quiet_install <- function(pkgs) {
  need <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
}
quiet_install(c("optparse","jsonlite"))

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

opt <- OptionParser()
opt <- add_option(opt, c("-d","--data_dir"), type="character",
                  default = Sys.getenv("DATA_DIR","."), help="root with 10x files")
opt <- add_option(opt, c("-o","--out_dir"),  type="character",
                  default = Sys.getenv("OUTPUT_DIR","docs"), help="output directory")
args <- parse_args(opt)

root_dir <- normalizePath(args$data_dir, mustWork = TRUE)
out_dir  <- normalizePath(args$out_dir,  mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, FALSE, TRUE)

STRICT_FAIL <- identical(Sys.getenv("STRICT_FAIL","0"), "1")
QUICK       <- identical(Sys.getenv("QUICK","0"), "1")
DO_ENRICH   <- identical(Sys.getenv("DO_ENRICH","0"), "1") && !QUICK  # 只在非快速模式時做富集

logi <- function(...) message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...))
append_file <- function(path, ...) cat(paste0(..., collapse=""), file = path, sep = "", append = TRUE)
fail_now <- function(msg) {
  append_file(file.path(out_dir, "_error.txt"), msg, "\n")
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  if (STRICT_FAIL) quit(status = 1, save = "no") else return(invisible(NULL))
}

set.seed(42)

## ========== 1. 主流程套件 ==========
# 只裝 CRAN；Bioc（clusterProfiler/org.Hs.eg.db）僅在 DO_ENRICH=1 時安裝
quiet_install(c("Seurat","SeuratObject","matrixStats","clue","patchwork",
                "readr","dplyr","ggplot2","future","ggrepel"))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(matrixStats); library(clue)
  library(patchwork); library(readr); library(dplyr); library(ggplot2)
  library(future); library(ggrepel)
})

# 多核（GitHub Actions 多為 2 核）
plan(multisession, workers = max(1, parallel::detectCores(logical = TRUE) - 0))

## ========== 2. 找 10x 檔案；必要時自動解壓 ==========
find_mtx <- function(root) {
  list.files(root, pattern = "_matrix.*\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
}

mat_files <- find_mtx(root_dir)
logi("Initial scan: found ", length(mat_files), " *_matrix*.mtx.gz")

auto_extract_archives <- function(path_root, max_rounds = 5) {
  rounds <- 0
  repeat {
    rounds <- rounds + 1
    z <- list.files(path_root, recursive = TRUE, full.names = TRUE,
                    pattern = "(\\.tar$)|(\\.tar\\.gz$)|(\\.tgz$)|(\\.zip$)")
    if (!length(z) || rounds > max_rounds) break
    logi("Auto-extract round ", rounds, " — archives: ", length(z))
    for (a in z) {
      ext <- tools::file_ext(a); dst <- dirname(a)
      logi("  -> ", a)
      if (ext %in% c("tar","gz","tgz")) utils::untar(a, exdir = dst, tar = "internal")
      else if (ext == "zip") utils::unzip(a, exdir = dst)
      # 可選：unlink(a)
    }
  }
}

if (length(mat_files) == 0) {
  logi("No matrix files yet; try auto-extracting any archives under: ", root_dir)
  auto_extract_archives(root_dir, max_rounds = 5)
  mat_files <- find_mtx(root_dir)
  logi("After auto-extract: found ", length(mat_files), " *_matrix*.mtx.gz")
}

if (length(mat_files) == 0) {
  append_file(file.path(out_dir, "_error.txt"),
              "No _matrix.mtx.gz under ", root_dir, "\nTop-level files:\n",
              paste0(list.files(root_dir, full.names = TRUE), collapse = "\n"), "\n")
  write_json(list(n_matrices = 0), file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  if (STRICT_FAIL) quit(status = 1, save = "no") else quit(status = 0, save = "no")
}

# 記錄前幾個找到的矩陣路徑
cat("First few matrices:\n", paste0(utils::head(mat_files, 10), collapse = "\n"), "\n",
    file = file.path(out_dir, "_matrices_found.txt"))

## ========== 3. 讀入與合併 ==========
pick1 <- function(dirp, pat) {
  fs <- list.files(dirp, pat, full.names = TRUE)
  if (length(fs)) fs[order(!grepl("\\.gz$", fs))][1] else NA_character_
}

read_one <- function(mtx_fp) {
  dirp   <- dirname(mtx_fp)
  prefix <- sub("_matrix.*\\.mtx\\.gz$", "", basename(mtx_fp))
  bc   <- pick1(dirp, paste0("^", prefix, ".*barcode[s]*.*tsv(\\.gz)?$"))
  feat <- pick1(dirp, paste0("^", prefix, ".*(features|genes).*tsv(\\.gz)?$"))
  if (is.na(bc) || is.na(feat)) {
    append_file(file.path(out_dir, "_error.txt"),
                "Missing barcodes/features for prefix: ", prefix, " in ", dirp, "\n",
                "dir listing:\n", paste0(list.files(dirp, full.names = TRUE), collapse = "\n"), "\n")
    return(NULL)
  }
  CreateSeuratObject(ReadMtx(mtx = mtx_fp, cells = bc, features = feat),
                     project = basename(dirname(mtx_fp)))
}

logi("Loading matrices into Seurat objects…")
objs <- lapply(mat_files, read_one)
objs <- Filter(Negate(is.null), objs)
if (!length(objs)) fail_now("All matrices failed pairing with barcodes/features; see _error.txt")

obj  <- if (length(objs) == 1) objs[[1]] else Reduce(function(x,y) merge(x, y), objs)

## ========== 4. QC ==========
obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
logi("QC filter: cells kept ", ncol(obj), "/", n0)
if (ncol(obj) < 200 && STRICT_FAIL) fail_now(paste0("Too few cells after QC: ", ncol(obj)))

## ========== 5. 快速模式抽樣（可選） ==========
if (QUICK) {
  keep <- min(5000, ncol(obj))
  if (ncol(obj) > keep) {
    set.seed(42)
    obj <- subset(obj, cells = sample(colnames(obj), keep))
    logi("QUICK mode: downsampled to ", keep, " cells")
  } else {
    logi("QUICK mode: total cells ≤ ", keep, ", skip downsampling")
  }
}

## ========== 6. SCTransform + PCA/UMAP/聚類 ==========
# 若已安裝 Bioc 的 glmGamPoi，則使用（更快、省記憶體）；否則使用預設
method_sel <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else NULL
if (is.null(method_sel)) {
  logi("SCTransform without glmGamPoi (install Bioc glmGamPoi to speed up).")
  obj <- SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE, verbose = FALSE)
} else {
  logi("SCTransform with glmGamPoi.")
  obj <- SCTransform(obj, vst.flavor = "v2", method = "glmGamPoi",
                     conserve.memory = TRUE, verbose = FALSE)
}

logi("PCA / UMAP / neighbors / clustering…")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE) |>
       RunUMAP(dims = 1:20, verbose = FALSE) |>
       FindNeighbors(dims = 1:20, verbose = FALSE) |>
       FindClusters(resolution = 0.6, verbose = FALSE)

## ========== 7. UMAP（OCP/OC 標註 + 圖上說明）==========
OCP_pos <- c("CSF1R","IRF8")
OC_pos  <- c("ACP5","CTSK","ATP6V0D2")

DefaultAssay(obj) <- "SCT"
legend_title <- "標準化表現量 (SCT；低→高)"

fmt_plots <- function(plist, legend_title) {
  lapply(plist, function(p) {
    p +
      labs(color = legend_title) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 9)
      )
  })
}

p_ocp <- FeaturePlot(
  obj, OCP_pos, combine = FALSE,
  cols = c("grey85","firebrick"),
  min.cutoff = "q05", max.cutoff = "q95"
)
p_oc  <- FeaturePlot(
  obj, OC_pos,  combine = FALSE,
  cols = c("grey85","steelblue"),
  min.cutoff = "q05", max.cutoff = "q95"
)

p_ocp <- fmt_plots(p_ocp, legend_title)
p_oc  <- fmt_plots(p_oc,  legend_title)

p_umap <- patchwork::wrap_plots(p_ocp, nrow = 1) /
          patchwork::wrap_plots(p_oc,  nrow = 1)

p_umap <- p_umap + patchwork::plot_annotation(
  title   = "OCP / OC 標記基因 UMAP",
  caption = paste(
    "色軸：該基因在每個細胞的標準化表現量（SCT），由低(灰)到高(色)。",
    "灰色點：該基因在該細胞幾乎無表現或低於最小截斷（min.cutoff）。"
  )
)

ggsave(file.path(plot_dir, "UMAP_markers.png"), p_umap, width = 10, height = 6.5, dpi = 300)

## ========== 8. 自動 OCP/OC 標註 + DEG + 火山圖 ==========
# 模組分數 + 分位數門檻：上 65% 判 OCP，下 35% 判 OC，其他為 Other
obj <- AddModuleScore(
  object   = obj,
  features = list(OCP_pos, OC_pos),
  name     = c("OCPscore","OCscore")
)
ocp_col <- grep("^OCPscore", colnames(obj@meta.data), value = TRUE)[1]
oc_col  <- grep("^OCscore",  colnames(obj@meta.data), value = TRUE)[1]

delta <- obj@meta.data[[ocp_col]] - obj@meta.data[[oc_col]]
qhi   <- as.numeric(quantile(delta, 0.65, na.rm = TRUE))
qlo   <- as.numeric(quantile(delta, 0.35, na.rm = TRUE))
label <- ifelse(delta >= qhi, "OCP", ifelse(delta <= qlo, "OC", "Other"))
obj$short_cluster <- factor(label, levels = c("OCP","OC","Other"))
Idents(obj) <- "short_cluster"

n_OCP <- sum(Idents(obj) == "OCP")
n_OC  <- sum(Idents(obj) == "OC")
logi(sprintf("Auto labels → OCP=%d, OC=%d, Other=%d", n_OCP, n_OC, sum(Idents(obj)=="Other")))
write.csv(table(obj$short_cluster), file.path(out_dir, "auto_labels_counts.csv"))

# DEG：OCP vs OC（Wilcoxon）
deg_path <- file.path(out_dir, "DE_OCP_vs_OC.csv")
if (n_OCP >= 100 && n_OC >= 100) {
  logi("FindMarkers OCP vs OC (wilcox)…")
  de <- FindMarkers(
    object = obj,
    ident.1 = "OCP",
    ident.2 = "OC",
    test.use = "wilcox",
    logfc.threshold = 0,      # 先不過濾，交給火山圖門檻
    min.pct = 0.1
  )
  de <- de[order(de$p_val_adj, -abs(de$avg_log2FC)), , drop = FALSE]
  write.csv(de, deg_path)
  logi("✓ DEG written: ", deg_path)

  # 火山圖
  vol <- transform(
    de,
    gene = rownames(de),
    neglog10padj = -log10(pmax(p_val_adj, .Machine$double.eps)),
    sig = (p_val_adj < 0.05) & (abs(avg_log2FC) > 0.25)
  )
  up_lab <- vol[vol$avg_log2FC > 0 & vol$sig, ]
  dn_lab <- vol[vol$avg_log2FC < 0 & vol$sig, ]
  up_lab <- head(up_lab[order(-up_lab$avg_log2FC, up_lab$p_val_adj), ], 10)
  dn_lab <- head(dn_lab[order(dn_lab$avg_log2FC,  dn_lab$p_val_adj), ], 10)
  lab_genes <- unique(c(up_lab$gene, dn_lab$gene))

  p_vol <- ggplot(vol, aes(x = avg_log2FC, y = neglog10padj)) +
    geom_point(aes(alpha = sig), size = 1.2) +
    scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4)) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05),    linetype = 2, linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = subset(vol, gene %in% lab_genes),
      aes(label = gene),
      size = 3, max.overlaps = 30, box.padding = 0.4
    ) +
    labs(x = "avg_log2FC (OCP vs OC)", y = "-log10(adjusted p-value)",
         title = "DEG Volcano: OCP vs OC") +
    theme_minimal(base_size = 12)

  vol_path <- file.path(plot_dir, "DE_volcano_OCP_vs_OC.png")
  ggsave(vol_path, p_vol, width = 8.5, height = 6.0, dpi = 300)
  logi("✓ Volcano plot: ", vol_path)

  # （可選）GO 富集：只在 DO_ENRICH=1 且有顯著上調基因時執行
  if (DO_ENRICH) {
    up_genes <- vol$gene[vol$sig & vol$avg_log2FC > 0]
    if (length(up_genes)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(c("clusterProfiler","org.Hs.eg.db","AnnotationDbi"),
                           update = FALSE, ask = FALSE)
      suppressPackageStartupMessages({
        library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
      })
      entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, up_genes, "ENTREZID", "SYMBOL", multiVals = "first")
      entrez <- na.omit(unname(entrez))
      if (length(entrez)) {
        ego <- clusterProfiler::enrichGO(entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
        write.csv(as.data.frame(ego), file.path(out_dir, "GO_up_OCP.csv"), row.names = FALSE)
        logi("✓ GO enrichment saved: ", file.path(out_dir, "GO_up_OCP.csv"))
      } else {
        logi("No ENTREZ IDs mapped; skip GO.")
      }
    } else {
      logi("No significant up genes for enrichment; skip GO.")
    }
  }
} else {
  logi("⚠ OCP/OC too few for DEG (need ≥100 each by default). Skipped DEG/Volcano.")
}

## ========== 9. 技術下采樣 ARI ==========
logi("Technical half-subsampling (5x) for ARI…")
ari_vec <- replicate(5, {
  cells_half <- sample(colnames(obj), floor(0.5 * ncol(obj)))
  obj_half <- subset(obj, cells = cells_half) |>
              FindNeighbors(dims = 1:20, verbose = FALSE) |>
              FindClusters(resolution = 0.6, verbose = FALSE)
  adjustedRandIndex(Idents(obj)[cells_half], Idents(obj_half))
})
write.csv(data.frame(ARI = ari_vec), file.path(out_dir, "tech_subsample_ARI.csv"), row.names = FALSE)

## ========== 10.（可選）Leave-One-Donor-Out ==========
if ("donor" %in% colnames(obj@meta.data)) {
  logi("LODO across donors…")
  donors <- unique(obj$donor)
  lodo <- sapply(donors, function(d) {
    train <- subset(obj, subset = donor != d)
    test  <- subset(obj, subset = donor == d)
    train <- RunPCA(train, npcs = 30, verbose = FALSE) |>
             FindNeighbors(dims = 1:20, verbose = FALSE) |>
             FindClusters(resolution = 0.6, verbose = FALSE)
    anchors <- FindTransferAnchors(reference = train, query = test, dims = 1:20)
    test    <- MapQuery(anchorset = anchors, query = test,
                        reference = train, refdata = list(cluster = "seurat_clusters"))
    mean(test$predicted.cluster == test$seurat_clusters)
  })
  write.csv(data.frame(donor = names(lodo), acc = as.numeric(lodo)),
            file.path(out_dir, "LODO_accuracy.csv"), row.names = FALSE)
}

## ========== 11. 摘要、保存 ==========
summary <- list(
  n_matrices     = length(mat_files),
  cells_before   = n0,
  cells_after    = ncol(obj),
  n_genes        = nrow(obj),
  n_clusters     = length(levels(Idents(obj))),
  ari_mean       = if (exists("ari_vec")) mean(ari_vec) else NA_real_,
  quick_mode     = QUICK,
  used_glmGamPoi = !is.null(method_sel),
  n_OCP          = if (exists("n_OCP")) n_OCP else NA_integer_,
  n_OC           = if (exists("n_OC"))  n_OC  else NA_integer_
)
write_json(summary, file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)

saveRDS(obj, file.path(out_dir, "combined_seurat.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

logi("Done. Outputs at: ", out_dir)
