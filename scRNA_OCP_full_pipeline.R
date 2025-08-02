###############################################################################
#  scRNA_OCP_full_pipeline.R
#  - QC → (可選) Doublet → (可選) 批次整合 → SCTransform → PCA/UMAP/Cluster
#  - OCP/OC 模組分數自動標註 → 單細胞 DE + 火山圖
#  - (可選) Pseudobulk DE（以捐贈者為樣本）→ 火山圖
#  - (可選) Slingshot 軌跡
#  - QUICK 確認模式：先快跑（≤5k）達標即全量
#
#  用法：
#    Rscript scRNA_OCP_full_pipeline.R --data_dir . --out_dir docs
#
#  可用環境變數（如未指定則使用括號內預設）：
#    DATA_DIR="."                 原始檔根目錄（含 *_RAW.tar 或 10x 三件組）
#    OUTPUT_DIR="docs"            輸出根目錄
#    QUICK="1"                    1=先快跑確認（達標後跑全量）；0=直接全量
#    MIN_OCP="100"                QUICK 通過門檻：OCP 細胞數
#    MIN_OC="100"                 QUICK 通過門檻：OC 細胞數
#    STRICT_FAIL="0"              1=遇到錯誤即 exit 1；0=寫 _error.txt 後繼續
#
#    DO_DOUBLET="0"               1=啟用 doublet 檢測（優先用 scDblFinder，否則 MAD 問題值）
#    DO_HARMONY="0"               1=用 Harmony 整合（需指定 BATCH_KEY），否則嘗試 Seurat anchors
#    DO_ANCHOR="0"                1=用 Seurat anchors 整合（需指定 BATCH_KEY）
#    BATCH_KEY=""                 批次欄位名稱（建議 'donor' 或 'sample' 或 'batch'）；留空則自動探測
#
#    DO_PSEUDOBULK="1"            1=若有 donor 欄位就做 pseudobulk DE（以捐贈者為單位）
#    PSEUDO_MIN_CELLS="50"        每個 donor × 群組的最小細胞數（不足則略過該 donor）
#
#    DO_ENRICH="0"                1=對 OCP 上調基因做 GO 富集（只在「全量」時執行；需 Bioc）
#    DO_TRAJECTORY="0"            1=Slingshot 軌跡（需 Bioc）
###############################################################################

## ===== 0. 安裝與載入基礎套件、參數 =====
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

STRICT_FAIL    <- identical(Sys.getenv("STRICT_FAIL","0"), "1")
QUICK_CONFIRM  <- identical(Sys.getenv("QUICK","1"), "1")
MIN_OCP        <- as.integer(Sys.getenv("MIN_OCP","100"))
MIN_OC         <- as.integer(Sys.getenv("MIN_OC","100"))
DO_DOUBLET     <- identical(Sys.getenv("DO_DOUBLET","0"), "1")
DO_HARMONY     <- identical(Sys.getenv("DO_HARMONY","0"), "1")
DO_ANCHOR      <- identical(Sys.getenv("DO_ANCHOR","0"), "1")
BATCH_KEY      <- Sys.getenv("BATCH_KEY","")
DO_PSEUDOBULK  <- identical(Sys.getenv("DO_PSEUDOBULK","1"), "1")
PSEUDO_MIN     <- as.integer(Sys.getenv("PSEUDO_MIN_CELLS","50"))
DO_ENRICH      <- identical(Sys.getenv("DO_ENRICH","0"), "1")
DO_TRAJECTORY  <- identical(Sys.getenv("DO_TRAJECTORY","0"), "1")

logi <- function(...) message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...))
append_file <- function(path, ...) cat(paste0(..., collapse=""), file = path, sep = "", append = TRUE)
fail_now <- function(msg) {
  append_file(file.path(out_dir, "_error.txt"), msg, "\n")
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  if (STRICT_FAIL) quit(status = 1, save = "no") else return(invisible(NULL))
}

set.seed(42)

## ===== 1. 主流程套件 =====
quiet_install(c("Seurat","SeuratObject","matrixStats","clue","patchwork",
                "readr","dplyr","ggplot2","ggrepel","future"))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(matrixStats); library(clue)
  library(patchwork); library(readr); library(dplyr); library(ggplot2); library(ggrepel)
  library(future)
})
plan(multisession, workers = max(1, parallel::detectCores(logical = TRUE) - 0))

## ===== 2. 找 10x 檔案；必要時自動解壓 =====
find_mtx <- function(root) {
  list.files(root, pattern = "_matrix.*\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
}
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
    }
  }
}

mat_files <- find_mtx(root_dir)
if (length(mat_files) == 0) {
  logi("No matrix files yet; try auto-extract…")
  auto_extract_archives(root_dir, max_rounds = 5)
  mat_files <- find_mtx(root_dir)
}
if (length(mat_files) == 0) {
  append_file(file.path(out_dir, "_error.txt"),
              "No _matrix.mtx.gz under ", root_dir, "\nTop-level files:\n",
              paste0(list.files(root_dir, full.names = TRUE), collapse = "\n"), "\n")
  write_json(list(n_matrices = 0), file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  if (STRICT_FAIL) quit(status = 1, save = "no") else quit(status = 0, save = "no")
}
cat("First matrices:\n", paste0(utils::head(mat_files, 10), collapse = "\n"), "\n",
    file = file.path(out_dir, "_matrices_found.txt"))

## ===== 3. 讀入與合併 =====
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
logi("Loading matrices…")
objs <- Filter(Negate(is.null), lapply(mat_files, read_one))
if (!length(objs)) fail_now("All matrices failed pairing barcodes/features; see _error.txt")
obj_full <- if (length(objs) == 1) objs[[1]] else Reduce(function(x,y) merge(x, y), objs)

## ===== 4. QC（+ 可選：Doublet）=====
obj_full[["percent.mt"]] <- PercentageFeatureSet(obj_full, "^MT-")
cells_before <- ncol(obj_full)
obj_full <- subset(obj_full, subset = nFeature_RNA > 500 & percent.mt < 20)
logi("QC filter: cells kept ", ncol(obj_full), "/", cells_before)
if (ncol(obj_full) < 200 && STRICT_FAIL) fail_now("Too few cells after QC")

# (可選) Doublet 檢測：優先 scDblFinder，否則 MAD 門檻
if (DO_DOUBLET) {
  logi("Doublet detection requested…")
  dbl_removed <- 0
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  try({
    BiocManager::install(c("SingleCellExperiment","scDblFinder"), update = FALSE, ask = FALSE)
    suppressPackageStartupMessages({
      library(SingleCellExperiment); library(scDblFinder)
    })
    sce <- as.SingleCellExperiment(obj_full)
    set.seed(42)
    sce <- scDblFinder(sce)
    obj_full$doublet <- as.logical(colData(sce)$scDblFinder.class == "doublet")
    dbl_removed <- sum(obj_full$doublet, na.rm = TRUE)
    obj_full <- subset(obj_full, cells = colnames(obj_full)[!obj_full$doublet])
    logi("scDblFinder removed doublets: ", dbl_removed)
  }, silent = TRUE)
  if (!"doublet" %in% colnames(obj_full@meta.data)) {
    # 後備：MAD-based 高度可疑 doublet（高 nCount & 高 nFeature）
    logi("scDblFinder unavailable; using MAD-based heuristic.")
    nUMI <- obj_full$nCount_RNA; nGene <- obj_full$nFeature_RNA
    thr_umi  <- median(nUMI)  + 3*mad(nUMI);  thr_gene <- median(nGene) + 3*mad(nGene)
    obj_full$doublet_mad <- (nUMI > thr_umi) & (nGene > thr_gene)
    dbl_removed <- sum(obj_full$doublet_mad, na.rm = TRUE)
    obj_full <- subset(obj_full, cells = colnames(obj_full)[!obj_full$doublet_mad])
    logi("MAD-heuristic removed cells: ", dbl_removed)
  }
  write.csv(table(dbl_removed = dbl_removed), file.path(out_dir, "doublet_removed_summary.csv"))
}

## ===== 5. 批次鍵推斷（若未提供）=====
detect_batch_key <- function(obj) {
  cand <- c("donor","sample","batch","orig.ident","patient","subject")
  hit <- cand[cand %in% colnames(obj@meta.data)][1]
  if (length(hit)) hit else ""
}
if (BATCH_KEY == "") BATCH_KEY <- detect_batch_key(obj_full)

## ===== 6. 核心分析函式（支援：Harmony / Anchors / 無整合）=====
run_core <- function(obj, out_dir_here, label_suffix = "", do_enrich = FALSE) {
  dir.create(out_dir_here, recursive = TRUE, showWarnings = FALSE)
  plot_here <- file.path(out_dir_here, "plots"); dir.create(plot_here, FALSE, TRUE)

  OCP_pos <- c("CSF1R","IRF8")
  OC_pos  <- c("ACP5","CTSK","ATP6V0D2")

  # 6.1 批次整合（可選）
  if (DO_HARMONY && BATCH_KEY %in% colnames(obj@meta.data)) {
    logi("[", label_suffix, "] Harmony integration by ", BATCH_KEY)
    quiet_install("harmony")
    suppressPackageStartupMessages(library(harmony))
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- RunHarmony(obj, group.by.vars = BATCH_KEY, dims.use = 1:30, verbose = FALSE)
    red  <- "harmony"
    dr_dims <- 1:20
  } else if (DO_ANCHOR && BATCH_KEY %in% colnames(obj@meta.data)) {
    logi("[", label_suffix, "] Seurat anchors integration by ", BATCH_KEY)
    objs <- SplitObject(obj, split.by = BATCH_KEY)
    objs <- lapply(objs, \(x) SCTransform(x, vst.flavor="v2", conserve.memory=TRUE, verbose=FALSE))
    feats <- SelectIntegrationFeatures(objs, nfeatures = 3000)
    objs  <- PrepSCTIntegration(objs, anchor.features = feats, verbose = FALSE)
    anchors <- FindIntegrationAnchors(objs, normalization.method = "SCT",
                                      anchor.features = feats, verbose = FALSE)
    obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
    DefaultAssay(obj) <- "integrated"
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    red  <- "pca"; dr_dims <- 1:20
  } else {
    # 無整合：直接 SCTransform
    method_sel <- if (requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else NULL
    if (is.null(method_sel)) {
      logi("[", label_suffix, "] SCTransform (no glmGamPoi)")
      obj <- SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE, verbose = FALSE)
    } else {
      logi("[", label_suffix, "] SCTransform (glmGamPoi)")
      obj <- SCTransform(obj, vst.flavor = "v2", method = "glmGamPoi",
                         conserve.memory = TRUE, verbose = FALSE)
    }
    DefaultAssay(obj) <- "SCT"
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    red  <- "pca"; dr_dims <- 1:20
  }

  # 6.2 UMAP/Cluster
  obj <- RunUMAP(obj, reduction = red, dims = dr_dims, verbose = FALSE) |>
         FindNeighbors(reduction = red, dims = dr_dims, verbose = FALSE) |>
         FindClusters(resolution = 0.6, verbose = FALSE)

  # 6.3 UMAP（帶色軸解釋）
  DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
  legend_title <- "標準化表現量 (SCT；低→高)"
  fmt_plots <- function(plist, legend_title) {
    lapply(plist, function(p) {
      p + labs(color = legend_title) +
        theme(legend.position = "right",
              legend.title = element_text(size = 11),
              legend.text  = element_text(size = 9))
    })
  }
  p_ocp <- FeaturePlot(obj, OCP_pos, combine = FALSE,
                       cols = c("grey85","firebrick"), min.cutoff = "q05", max.cutoff = "q95")
  p_oc  <- FeaturePlot(obj, OC_pos,  combine = FALSE,
                       cols = c("grey85","steelblue"), min.cutoff = "q05", max.cutoff = "q95")
  p_ocp <- fmt_plots(p_ocp, legend_title); p_oc <- fmt_plots(p_oc, legend_title)
  p_umap <- patchwork::wrap_plots(p_ocp, nrow = 1) / patchwork::wrap_plots(p_oc, nrow = 1)
  p_umap <- p_umap + patchwork::plot_annotation(
    title   = paste0("OCP / OC 標記基因 UMAP", label_suffix),
    caption = "色軸：標準化表現量（SCT），由低(灰)到高(色)。灰點＝低/無表現（或低於 min.cutoff）。"
  )
  ggsave(file.path(plot_here, paste0("UMAP_markers", label_suffix, ".png")),
         p_umap, width = 10, height = 6.5, dpi = 300)

  # 6.4 OCP/OC 自動標註（模組分數 + 分位數）
  obj <- AddModuleScore(obj, features = list(OCP_pos, OC_pos), name = c("OCPscore","OCscore"))
  ocp_col <- grep("^OCPscore", colnames(obj@meta.data), value = TRUE)[1]
  oc_col  <- grep("^OCscore",  colnames(obj@meta.data), value = TRUE)[1]
  delta <- obj@meta.data[[ocp_col]] - obj@meta.data[[oc_col]]
  qhi   <- as.numeric(quantile(delta, 0.65, na.rm = TRUE))
  qlo   <- as.numeric(quantile(delta, 0.35, na.rm = TRUE))
  obj$short_cluster <- factor(ifelse(delta >= qhi, "OCP",
                              ifelse(delta <= qlo, "OC", "Other")),
                              levels = c("OCP","OC","Other"))
  Idents(obj) <- "short_cluster"
  tbl <- table(obj$short_cluster)
  write.csv(tbl, file.path(out_dir_here, paste0("auto_labels_counts", label_suffix, ".csv")))
  n_OCP <- unname(tbl["OCP"]); n_OC <- unname(tbl["OC"]); n_OCP[is.na(n_OCP)] <- 0; n_OC[is.na(n_OC)] <- 0
  logi("[", label_suffix, "] auto labels: OCP=", n_OCP, "  OC=", n_OC)

  # 6.5 單細胞 DE（Wilcoxon）+ 火山圖
  deg_ok <- FALSE
  if (n_OCP >= MIN_OCP && n_OC >= MIN_OC) {
    de <- FindMarkers(obj, ident.1 = "OCP", ident.2 = "OC",
                      test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)
    de <- de[order(de$p_val_adj, -abs(de$avg_log2FC)), , drop = FALSE]
    write.csv(de, file.path(out_dir_here, paste0("DE_OCP_vs_OC_singlecell", label_suffix, ".csv")))
    deg_ok <- TRUE

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
      ggrepel::geom_text_repel(data = subset(vol, gene %in% lab_genes),
                               aes(label = gene), size = 3, max.overlaps = 30, box.padding = 0.4) +
      labs(x = "avg_log2FC (OCP vs OC)", y = "-log10(adj. p)", title = paste0("Single-cell DEG", label_suffix)) +
      theme_minimal(base_size = 12)
    ggsave(file.path(plot_here, paste0("DE_volcano_singlecell", label_suffix, ".png")),
           p_vol, width = 8.5, height = 6.0, dpi = 300)
  } else {
    logi("[", label_suffix, "] DE skipped: not enough OCP/OC (need ≥", MIN_OCP, " each).")
  }

  # 6.6 （可選）Pseudobulk DE（以 donor 為樣本；無 Bioc 依賴）
  pseudo_ok <- FALSE
  if (DO_PSEUDOBULK && "donor" %in% colnames(obj@meta.data)) {
    logi("[", label_suffix, "] Pseudobulk DE by donor…")
    # 取原始 counts 以 donor×label 彙總成「樣本」
    counts <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
    meta   <- obj@meta.data %>% dplyr::mutate(label = as.character(short_cluster))
    # 只留 OCP/OC
    keep_cells <- rownames(meta)[meta$label %in% c("OCP","OC")]
    counts <- counts[, keep_cells, drop = FALSE]
    meta   <- meta[keep_cells, , drop = FALSE]
    # 篩掉 donor×label 細胞數不足者
    meta$dl <- paste(meta$donor, meta$label, sep = "|")
    tab <- table(meta$dl)
    ok_dl <- names(tab)[tab >= PSEUDO_MIN]
    meta <- meta[meta$dl %in% ok_dl, , drop = FALSE]
    counts <- counts[, rownames(meta), drop = FALSE]
    if (ncol(counts) > 0) {
      # 對每個 donor×label 做 colSums 當作「樣本」的 pseudobulk counts
      dl_levels <- unique(meta$dl)
      pseudo_mat <- sapply(dl_levels, function(k) {
        cols <- rownames(meta)[meta$dl == k]
        Matrix::rowSums(counts[, cols, drop = FALSE])
      })
      colnames(pseudo_mat) <- dl_levels
      # 轉回 data.frame：為每個 donor 只取 OCP vs OC 一對
      info <- do.call(rbind, strsplit(colnames(pseudo_mat), "\\|"))
      colnames(info) <- c("donor","label")
      info <- as.data.frame(info)
      # 對每個基因：以 donor 為配對單位，做配對 Wilcoxon（若同一 donor 同時有 OCP 與 OC）
      donors <- unique(info$donor)
      # 先把 counts CPM 後 log1p
      lib <- Matrix::colSums(pseudo_mat)
      cpm <- t(t(pseudo_mat) / pmax(lib, 1)) * 1e6
      logcpm <- log1p(cpm)
      # 找同時有兩組的 donor
      good_donors <- donors[sapply(donors, function(d) {
        sum(info$donor == d & info$label == "OCP") == 1 &&
        sum(info$donor == d & info$label == "OC")  == 1
      })]
      if (length(good_donors) >= 3) {  # 至少 3 個 donor
        ocp_cols <- paste0(good_donors, "|OCP")
        oc_cols  <- paste0(good_donors, "|OC")
        ocp_mat <- logcpm[, ocp_cols, drop = FALSE]
        oc_mat  <- logcpm[, oc_cols,  drop = FALSE]
        # 配對 Wilcoxon：每基因在 good_donors 間比較 OCP vs OC
        pvals <- apply(ocp_mat - oc_mat, 1, function(v) {
          tryCatch(wilcox.test(v, mu = 0, paired = TRUE)$p.value, error = function(e) 1)
        })
        fc <- rowMeans(ocp_mat) - rowMeans(oc_mat)  # log 差近似 logFC
        res <- data.frame(avg_logFC_pb = fc, p_val = pvals,
                          p_val_adj = p.adjust(pvals, method = "BH"))
        res <- res[order(res$p_val_adj, -abs(res$avg_logFC_pb)), , drop = FALSE]
        write.csv(res, file.path(out_dir_here, paste0("DE_OCP_vs_OC_pseudobulk", label_suffix, ".csv")))
        pseudo_ok <- TRUE

        vol2 <- transform(
          res,
          gene = rownames(res),
          neglog10padj = -log10(pmax(p_val_adj, .Machine$double.eps)),
          sig = (p_val_adj < 0.05) & (abs(avg_logFC_pb) > 0.25)
        )
        up2 <- vol2[vol2$avg_logFC_pb > 0 & vol2$sig, ]
        dn2 <- vol2[vol2$avg_logFC_pb < 0 & vol2$sig, ]
        up2 <- head(up2[order(-up2$avg_logFC_pb, up2$p_val_adj), ], 10)
        dn2 <- head(dn2[order(dn2$avg_logFC_pb,  dn2$p_val_adj), ], 10)
        lab2 <- unique(c(rownames(up2), rownames(dn2)))
        p_vol2 <- ggplot(vol2, aes(x = avg_logFC_pb, y = neglog10padj)) +
          geom_point(aes(alpha = sig), size = 1.2) +
          scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4)) +
          geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, linewidth = 0.3) +
          geom_hline(yintercept = -log10(0.05),    linetype = 2, linewidth = 0.3) +
          ggrepel::geom_text_repel(data = subset(vol2, gene %in% lab2),
                                   aes(label = gene), size = 3, max.overlaps = 30, box.padding = 0.4) +
          labs(x = "avg_logFC (pseudo; OCP-OC)", y = "-log10(adj. p)",
               title = paste0("Pseudobulk DEG", label_suffix)) +
          theme_minimal(base_size = 12)
        ggsave(file.path(plot_here, paste0("DE_volcano_pseudobulk", label_suffix, ".png")),
               p_vol2, width = 8.5, height = 6.0, dpi = 300)
      } else {
        logi("[", label_suffix, "] Pseudobulk skipped: donors with both OCP & OC < 3.")
      }
    } else {
      logi("[", label_suffix, "] Pseudobulk skipped: no cells after filters.")
    }
  }

  # 6.7 （可選）Slingshot 軌跡
  if (DO_TRAJECTORY) {
    logi("[", label_suffix, "] Slingshot trajectory…")
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    ok_traj <- try({
      BiocManager::install(c("SingleCellExperiment","slingshot"), update = FALSE, ask = FALSE)
      suppressPackageStartupMessages({ library(SingleCellExperiment); library(slingshot) })
      sce <- as.SingleCellExperiment(obj)
      # 以 UMAP 為空間，並以 OCP → OC 當作大致方向（起點群給 OCP）
      clus <- as.factor(obj$short_cluster)
      reducedDim(sce, "UMAP") <- Embeddings(obj, "umap")
      sce <- slingshot(sce, clusterLabels = clus, reducedDim = "UMAP",
                       start.clus = "OCP", end.clus = "OC")
      # 繪圖
      um <- as.data.frame(Embeddings(obj, "umap"))
      colnames(um) <- c("UMAP1","UMAP2")
      um$label <- obj$short_cluster
      p_traj <- ggplot(um, aes(UMAP1, UMAP2, color = label)) +
        geom_point(size = 0.4, alpha = 0.6) +
        theme_minimal() + labs(title = paste0("Slingshot trajectory", label_suffix))
      # 將曲線疊上
      crv <- slingCurves(sce)
      for (i in seq_along(crv)) {
        dfc <- as.data.frame(crv[[i]]$s[curv[[i]]$ord, ]) # 有些版本是 crv[[i]]$s
      }
      # 為避免版本差異，改用 lines() 輸出 PNG
      png(file.path(plot_here, paste0("trajectory_slingshot", label_suffix, ".png")),
          width = 900, height = 650, res = 120)
      plot(um$UMAP1, um$UMAP2, col = scales::hue_pal()(3)[as.numeric(um$label)],
           pch = 16, cex = 0.4, xlab = "UMAP1", ylab = "UMAP2",
           main = paste0("Slingshot trajectory", label_suffix))
      lines(SlingshotDataSet(sce), lwd = 2, col = "black")
      dev.off()
      TRUE
    }, silent = TRUE)
    if (inherits(ok_traj, "try-error") || isFALSE(ok_traj)) {
      logi("[", label_suffix, "] Slingshot failed or unavailable; skipped.")
    }
  }

  # 6.8（可選）GO 富集：僅在全量且你開 DO_ENRICH 時做（用單細胞 DE 上調）
  if (do_enrich && deg_ok) {
    up <- rownames(de)[de$p_val_adj < 0.05 & de$avg_log2FC > 0.25]
    if (length(up)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(c("clusterProfiler","org.Hs.eg.db","AnnotationDbi"),
                           update = FALSE, ask = FALSE)
      suppressPackageStartupMessages({
        library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
      })
      eg <- AnnotationDbi::mapIds(org.Hs.eg.db, up, "ENTREZID", "SYMBOL", multiVals = "first")
      eg <- na.omit(unname(eg))
      if (length(eg)) {
        ego <- enrichGO(eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
        write.csv(as.data.frame(ego), file.path(out_dir_here, paste0("GO_up_OCP", label_suffix, ".csv")),
                  row.names = FALSE)
      }
    }
  }

  # 摘要與保存
  summary <- list(
    cells_after    = ncol(obj),
    n_genes        = nrow(obj),
    n_clusters     = length(levels(Idents(obj))),
    n_OCP          = as.integer(n_OCP),
    n_OC           = as.integer(n_OC),
    deg_singlecell = deg_ok,
    pseudobulk     = pseudo_ok,
    batch_method   = if (DO_HARMONY && BATCH_KEY %in% colnames(obj@meta.data)) "Harmony"
                     else if (DO_ANCHOR && BATCH_KEY %in% colnames(obj@meta.data)) "Anchors"
                     else "None",
    batch_key      = ifelse(BATCH_KEY == "", NA, BATCH_KEY)
  )
  write_json(summary, file.path(out_dir_here, paste0("_summary", label_suffix, ".json")),
             auto_unbox = TRUE, pretty = TRUE)
  saveRDS(obj, file.path(out_dir_here, paste0("combined_seurat", label_suffix, ".rds")))
  sink(file.path(out_dir_here, paste0("sessionInfo", label_suffix, ".txt"))); print(sessionInfo()); sink()

  list(obj = obj, summary = summary, n_OCP = n_OCP, n_OC = n_OC,
       out_dir = out_dir_here, plot_dir = plot_here)
}

## ===== 7. QUICK 確認 → 通過自動跑全量 =====
if (QUICK_CONFIRM) {
  logi("QUICK confirm: run quick subset first…")
  quick_out <- file.path(out_dir, "quick_check")
  keep <- min(5000, ncol(obj_full))
  obj_quick <- if (ncol(obj_full) > keep) subset(obj_full, cells = sample(colnames(obj_full), keep)) else obj_full
  res_quick <- run_core(obj_quick, out_dir_here = quick_out, label_suffix = "_quick", do_enrich = FALSE)
  if (res_quick$n_OCP >= MIN_OCP && res_quick$n_OC >= MIN_OC) {
    logi("Quick PASSED (OCP=", res_quick$n_OCP, ", OC=", res_quick$n_OC, "). Run FULL…")
    res_full <- run_core(obj_full, out_dir_here = out_dir, label_suffix = "", do_enrich = DO_ENRICH)
    logi("Done. FULL outputs at: ", res_full$out_dir, " (quick preview at: ", res_quick$out_dir, ")")
  } else {
    msg <- paste0("Quick NOT passed: OCP=", res_quick$n_OCP, ", OC=", res_quick$n_OC,
                  " (need ≥", MIN_OCP, " each). See quick_check outputs.")
    append_file(file.path(out_dir, "_error.txt"), msg, "\n")
    logi(msg)
    if (STRICT_FAIL) quit(status = 1, save = "no")
  }
} else {
  logi("QUICK=0 → run FULL directly…")
  res_full <- run_core(obj_full, out_dir_here = out_dir, label_suffix = "", do_enrich = DO_ENRICH)
  logi("Done. FULL outputs at: ", res_full$out_dir)
}
