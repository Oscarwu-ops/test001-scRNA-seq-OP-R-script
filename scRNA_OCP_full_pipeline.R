###############################################################################
#  scRNA_OCP_full_pipeline.R  --  CI-friendly, auto-extract, fast defaults
#
#  用法 (GitHub Actions / 本機皆可)：
#    Rscript scRNA_OCP_full_pipeline.R --data_dir . --out_dir docs
#
#  重要環境變數 (也可改用 CLI 參數)：
#    DATA_DIR   = "."     # 原始檔案所在根目錄；會遞迴解壓並尋找 *_matrix*.mtx.gz
#    OUTPUT_DIR = "docs"  # 所有輸出放這裡
#    QUICK      = "1"     # "1"=快速模式(抽樣≤5000細胞)；"0"=全量
#    STRICT_FAIL= "0"     # "1"=遇到不合理數據(如細胞太少)直接 exit 1；"0"=只寫 _error.txt
#    DO_ENRICH  = "0"     # "1"=做 GO 富集 (會安裝 Bioc；建議只在正式產物時開啟)
#
#  主要流程：
#   1) 自動(必要時)把 GSE*_RAW.tar / *.tar.gz / *.zip 解到 data_dir 之下
#   2) 收集 10x 三件組；讀入並合併
#   3) 基本 QC（nFeature_RNA、percent.mt）
#   4) SCTransform（若有 glmGamPoi 則用以加速，否則採預設）
#   5) PCA → UMAP → 鄰居圖 → 聚類
#   6) 繪製 OCP/OC markers UMAP、技術下采樣 ARI、(可選) LODO、(可選) DE/GO
#   7) 寫出 _summary.json、run artifacts（圖、csv、rds、sessionInfo）
###############################################################################

## ========== 0. 前置：參數與小工具 ==========
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
  # 另外輸出 sessionInfo 便於診斷
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  if (STRICT_FAIL) quit(status = 1, save = "no") else return(invisible(NULL))
}

set.seed(42)

## ========== 1. 套件 ==========
# 主流程的 CRAN 套件
quiet_install(c("Seurat","SeuratObject","matrixStats","clue","patchwork","readr","dplyr","ggplot2","future"))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(matrixStats); library(clue)
  library(patchwork); library(readr); library(dplyr); library(ggplot2); library(future)
})

# 多核（GitHub Actions 多半 2 核）
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
    # 找到所有巢狀壓縮
    z <- list.files(path_root, recursive = TRUE, full.names = TRUE,
                    pattern = "(\\.tar$)|(\\.tar\\.gz$)|(\\.tgz$)|(\\.zip$)")
    if (!length(z) || rounds > max_rounds) break
    logi("Auto-extract round ", rounds, " — archives: ", length(z))
    for (a in z) {
      ext <- tools::file_ext(a); dst <- dirname(a)
      logi("  -> ", a)
      if (ext %in% c("tar","gz","tgz")) utils::untar(a, exdir = dst, tar = "internal")
      else if (ext == "zip") utils::unzip(a, exdir = dst)
      # 如需節省空間可刪除原壓縮：try(unlink(a), silent = TRUE)
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
  # 輸出目錄快照，幫助你在 Actions 上判斷
  append_file(file.path(out_dir, "_error.txt"),
              "No _matrix.mtx.gz under ", root_dir, "\nTop-level files:\n",
              paste0(list.files(root_dir, full.names = TRUE), collapse = "\n"), "\n")
  if (STRICT_FAIL) quit(status = 1, save = "no") else {
    # 仍寫出空的 summary 與 sessionInfo 方便你看
    write_json(list(n_matrices = 0), file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)
    sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
    quit(status = 0, save = "no")
  }
}

# 記錄前幾個找到的矩陣
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
if (!length(objs)) fail_now("All matrices failed to pair with barcodes/features; see _error.txt")

obj  <- if (length(objs) == 1) objs[[1]] else Reduce(function(x,y) merge(x, y), objs)

## ========== 4. QC ==========
obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
logi("QC filter: cells kept ", ncol(obj), "/", n0)
if (ncol(obj) < 200 && STRICT_FAIL) fail_now(paste0("Too few cells after QC: ", ncol(obj)))

## ========== 5. 快速模式抽樣 (可選) ==========
if (QUICK) {
  keep <- min(5000, ncol(obj))  # 視資料量調整上限
  if (ncol(obj) > keep) {
    set.seed(42)
    obj <- subset(obj, cells = sample(colnames(obj), keep))
    logi("QUICK mode: downsampled to ", keep, " cells")
  } else {
    logi("QUICK mode: total cells ≤ ", keep, ", skip downsampling")
  }
}

## ========== 6. SCTransform + PCA/UMAP/聚類 ==========
# 若有 glmGamPoi 則使用（更快、省記憶體）；沒有就用預設
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

## ========== 7. 繪圖：OCP / OC markers ==========
OCP_pos <- c("CSF1R","IRF8")
OC_pos  <- c("ACP5","CTSK","ATP6V0D2")
p <- (FeaturePlot(obj, OCP_pos, combine = FALSE) |> patchwork::wrap_plots()) /
     (FeaturePlot(obj, OC_pos,  combine = FALSE) |> patchwork::wrap_plots())
ggsave(file.path(plot_dir, "UMAP_markers.png"), p, width = 9, height = 6, dpi = 300)

## ========== 8. 技術下采樣 ARI (五次) ==========
logi("Technical half-subsampling (5x) for ARI…")
quiet_install("clue"); suppressPackageStartupMessages(library(clue))
ari_vec <- replicate(5, {
  cells_half <- sample(colnames(obj), floor(0.5 * ncol(obj)))
  obj_half <- subset(obj, cells = cells_half) |>
              FindNeighbors(dims = 1:20, verbose = FALSE) |>
              FindClusters(resolution = 0.6, verbose = FALSE)
  adjustedRandIndex(Idents(obj)[cells_half], Idents(obj_half))
})
write.csv(data.frame(ARI = ari_vec), file.path(out_dir, "tech_subsample_ARI.csv"), row.names = FALSE)

## ========== 9. (可選) LODO ==========
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

## ========== 10. (可選) 差異表達與 GO 富集 ==========
# 假設你會在外部（或上一步）把 OCP/OC label 寫到 meta 欄位 short_cluster
if ("short_cluster" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "short_cluster"
  if (all(c("OCP","OC") %in% levels(Idents(obj)))) {
    logi("FindMarkers OCP vs OC (wilcox)…")
    de <- FindMarkers(obj, ident.1 = "OCP", ident.2 = "OC", test.use = "wilcox")
    write.csv(de, file.path(out_dir, "DE_OCP_vs_OC.csv"))
    # 可選 GO 富集：只在 DO_ENRICH=1 且有顯著基因時執行
    if (DO_ENRICH && nrow(de) > 0) {
      up <- rownames(de)[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05]
      if (length(up)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(c("clusterProfiler","org.Hs.eg.db","AnnotationDbi"),
                             update = FALSE, ask = FALSE)
        suppressPackageStartupMessages({
          library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
        })
        entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, up, "ENTREZID", "SYMBOL", multiVals = "first")
        entrez <- na.omit(unname(entrez))
        if (length(entrez)) {
          ego <- clusterProfiler::enrichGO(entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
          write.csv(as.data.frame(ego), file.path(out_dir, "GO_up_OCP.csv"), row.names = FALSE)
        }
      }
    }
  } else {
    logi("short_cluster exists but lacks both OCP & OC levels; skip DE/GO.")
  }
} else {
  logi("No 'short_cluster' in metadata; skip DE/GO.")
}

## ========== 11. 摘要、保存與 session ==========
summary <- list(
  n_matrices   = length(mat_files),
  cells_before = n0,
  cells_after  = ncol(obj),
  n_genes      = nrow(obj),
  n_clusters   = length(levels(Idents(obj))),
  ari_mean     = if (exists("ari_vec")) mean(ari_vec) else NA_real_,
  quick_mode   = QUICK,
  used_glmGamPoi = !is.null(method_sel)
)
write_json(summary, file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)

saveRDS(obj, file.path(out_dir, "combined_seurat.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

logi("Done. Outputs at: ", out_dir)
