###############################################################################
#  scRNA_OCP_full_pipeline.R  --  GitHub Actions friendly version
#  用法：
#    本地端: Rscript scRNA_OCP_full_pipeline.R --data_dir data --out_dir docs
#    Actions: 由 workflow 設定 DATA_DIR / OUTPUT_DIR 或直接用預設
###############################################################################

## ========== 0. 參數解析 ==========
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE))
    install.packages("optparse", repos = "https://cloud.r-project.org")
  library(optparse)
})

opt <- OptionParser()
opt <- add_option(opt, c("-d", "--data_dir"), type = "character",
                  default = Sys.getenv("DATA_DIR", "."),
                  help = "root directory that contains GSE folders [default: %default]")
opt <- add_option(opt, c("-o", "--out_dir"),  type = "character",
                  default = Sys.getenv("OUTPUT_DIR", "docs"),
                  help = "output directory [default: %default]")
args <- parse_args(opt)

root_dir <- normalizePath(args$data_dir, mustWork = TRUE)
out_dir  <- normalizePath(args$out_dir,  mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, FALSE, TRUE)

## ========== 1. 套件安裝 & 載入 ==========
install_cran <- function(pkgs) {
  need <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
}
install_bioc <- function(pkgs) {
  if (!length(pkgs)) return()
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(pkgs, update = FALSE, ask = FALSE)
}

install_cran(c("Seurat", "SeuratObject", "matrixStats", "clue",
               "patchwork", "readr", "dplyr", "ggplot2"))
install_bioc(c("clusterProfiler", "org.Hs.eg.db"))

suppressPackageStartupMessages({
  library(Seurat)
  library(matrixStats)
  library(clue)
  library(patchwork)
  library(readr)
  library(dplyr)     # 含 %>% 管線
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(42)

## ========== 2. 讀取 10x 矩陣 ==========
# 會在 root_dir 之下遞迴尋找 *_matrix*.mtx.gz，並依同資料夾中的 barcodes/features 配對
mat_files <- list.files(root_dir, pattern = "_matrix.*\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
if (length(mat_files) == 0) stop("No _matrix.mtx.gz files found under ", root_dir)

read_one <- function(mtx_fp) {
  dirp   <- dirname(mtx_fp)
  prefix <- sub("_matrix.*\\.mtx\\.gz$", "", basename(mtx_fp))
  pick1  <- function(pat) {
    fs <- list.files(dirp, pat, full.names = TRUE)
    if (length(fs)) fs[order(!grepl("\\.gz$", fs))][1] else NA_character_
  }
  bc   <- pick1(paste0("^", prefix, ".*barcode[s]*.*tsv(\\.gz)?$"))
  feat <- pick1(paste0("^", prefix, ".*(features|genes).*tsv(\\.gz)?$"))
  if (is.na(bc) || is.na(feat)) stop("Missing barcodes/features for ", prefix)
  CreateSeuratObject(ReadMtx(mtx = mtx_fp, cells = bc, features = feat),
                     project = basename(dirname(mtx_fp)))
}

message(">> Loading matrices …")
objs <- lapply(mat_files, read_one)
obj  <- if (length(objs) == 1) objs[[1]] else Reduce(function(x, y) merge(x, y), objs)

## ========== 3. 基礎 QC ==========
obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)

## ========== 4. SCTransform ==========
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)

## ========== 5. PCA / UMAP / Clustering ==========
obj <- RunPCA(obj, npcs = 30, verbose = FALSE) |>
       RunUMAP(dims = 1:20) |>
       FindNeighbors(dims = 1:20) |>
       FindClusters(resolution = 0.6)

## 5.1 UMAP + 標記圖
OCP_pos <- c("CSF1R", "IRF8"); OC_pos <- c("ACP5", "CTSK", "ATP6V0D2")
p <- (FeaturePlot(obj, OCP_pos, combine = FALSE) %>% wrap_plots()) /
     (FeaturePlot(obj, OC_pos,  combine = FALSE) %>% wrap_plots())
ggsave(file.path(plot_dir, "UMAP_markers.png"), p, width = 9, height = 6, dpi = 300)

## ========== 6. 技術下采樣 ARI ==========
ari_vec <- replicate(5, {
  cells_half <- sample(colnames(obj), floor(0.5 * ncol(obj)))
  obj_half <- subset(obj, cells = cells_half) |>
              FindNeighbors(dims = 1:20) |>
              FindClusters(resolution = 0.6)
  adjustedRandIndex(Idents(obj)[cells_half], Idents(obj_half))
})
write.csv(data.frame(ARI = ari_vec),
          file.path(out_dir, "tech_subsample_ARI.csv"), row.names = FALSE)

## ========== 7. (可選) Leave-One-Donor-Out ==========
if ("donor" %in% colnames(obj@meta.data)) {
  donors <- unique(obj$donor)
  lodo <- sapply(donors, function(d) {
    train <- subset(obj, subset = donor != d)
    test  <- subset(obj, subset = donor == d)
    train <- RunPCA(train, npcs = 30, verbose = FALSE) |>
             FindNeighbors(dims = 1:20) |>
             FindClusters(resolution = 0.6)
    anchors <- FindTransferAnchors(reference = train, query = test, dims = 1:20)
    test    <- MapQuery(anchorset = anchors, query = test,
                        reference = train, refdata = list(cluster = "seurat_clusters"))
    mean(test$predicted.cluster == test$seurat_clusters)
  })
  write.csv(data.frame(donor = names(lodo), acc = as.numeric(lodo)),
            file.path(out_dir, "LODO_accuracy.csv"), row.names = FALSE)
}

## ========== 8. 差異表達：OCP vs OC ==========
# 提醒：請確認你已用 marker/比例或手動方式把群集命名成 OCP / OC
# 這裡示範：若剛好 cluster 0/1 對應 OCP/OC，則進行 DE；否則跳過（避免誤用）
if (!"short_cluster" %in% colnames(obj@meta.data)) {
  # 你可以在這裡依需要修改對應
  # obj$short_cluster <- recode(Idents(obj), `0` = "OCP", `1` = "OC")
  message("No 'short_cluster' labels found; skip DE. Please set Idents to 'OCP'/'OC' before rerun.")
} else {
  Idents(obj) <- "short_cluster"
}

if (all(c("OCP", "OC") %in% levels(Idents(obj)))) {
  de <- FindMarkers(obj, ident.1 = "OCP", ident.2 = "OC", test.use = "wilcox")
  write.csv(de, file.path(out_dir, "DE_OCP_vs_OC.csv"))
}

## ========== 9. GO 富集（可選） ==========
if (exists("de") && is.data.frame(de) && nrow(de) > 0) {
  up <- rownames(de)[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05]
  if (length(up)) {
    entrez <- mapIds(org.Hs.eg.db, up, "ENTREZID", "SYMBOL", multiVals = "first")
    entrez <- na.omit(unname(entrez))
    if (length(entrez)) {
      ego <- enrichGO(entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP")
      write.csv(as.data.frame(ego), file.path(out_dir, "GO_up_OCP.csv"), row.names = FALSE)
    }
  }
}

## ========== 10. 存檔 & 紀錄 ==========
saveRDS(obj, file.path(out_dir, "combined_seurat.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

message("Pipeline finished. Outputs saved to: ", out_dir)
