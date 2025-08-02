###############################################################################
#  scRNA_OCP_full_pipeline.R  --  CI-friendly with diagnostics (non-strict)
#  Run:
#    Rscript scRNA_OCP_full_pipeline.R --data_dir . --out_dir docs
###############################################################################

## ========== 0. CLI options ==========
quiet_install <- function(pkgs) {
  need <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
}
quiet_install(c("optparse","jsonlite"))
suppressPackageStartupMessages({
  library(optparse); library(jsonlite)
})

opt <- OptionParser()
opt <- add_option(opt, c("-d","--data_dir"), type="character",
                  default=Sys.getenv("DATA_DIR","."), help="root with 10x files")
opt <- add_option(opt, c("-o","--out_dir"),  type="character",
                  default=Sys.getenv("OUTPUT_DIR","docs"), help="output dir")
args <- parse_args(opt)

root_dir <- normalizePath(args$data_dir, mustWork = TRUE)
out_dir  <- normalizePath(args$out_dir,  mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, FALSE, TRUE)

STRICT_FAIL <- identical(Sys.getenv("STRICT_FAIL","0"), "1")

logi <- function(...) message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...))
append_file <- function(path, ...) cat(paste0(..., collapse = ""), file = path, sep = "", append = TRUE)

## ========== 1. packages ==========
quiet_install(c("Seurat","SeuratObject","matrixStats","clue","patchwork","readr","dplyr","ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler","org.Hs.eg.db"), update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(matrixStats); library(clue); library(patchwork)
  library(readr);  library(dplyr);       library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(42)

## ========== 2. locate 10x files ==========
logi("Scanning for *_matrix*.mtx.gz under: ", root_dir)
mat_files <- list.files(root_dir, pattern = "_matrix.*\\.mtx\\.gz$", recursive = TRUE, full.names = TRUE)
logi("Found ", length(mat_files), " matrix file(s).")
if (length(mat_files) > 0) {
  cat("First few matrices:\n", paste0(head(mat_files, 10), collapse = "\n"), "\n",
      file = file.path(out_dir, "_matrices_found.txt"))
} else {
  append_file(file.path(out_dir, "_error.txt"), "No _matrix.mtx.gz under ", root_dir, "\n")
  append_file(file.path(out_dir, "_error.txt"), "Top-level files:\n",
              paste0(list.files(root_dir, full.names = TRUE), collapse = "\n"), "\n")
  quit(status = if (STRICT_FAIL) 1 else 0, save = "no")
}

pick1 <- function(dirp, pat) {
  fs <- list.files(dirp, pat, full.names = TRUE)
  if (length(fs)) fs[order(!grepl("\\.gz$", fs))][1] else NA_character_
}

read_one <- function(mtx_fp) {
  dirp <- dirname(mtx_fp)
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
if (!length(objs)) {
  append_file(file.path(out_dir, "_error.txt"), "All matrices failed to pair with barcodes/features.\n")
  quit(status = if (STRICT_FAIL) 1 else 0, save = "no")
}
obj  <- if (length(objs) == 1) objs[[1]] else Reduce(function(x,y) merge(x, y), objs)

## ========== 3. QC ==========
obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
logi("QC filter: cells kept ", ncol(obj), "/", n0)

## ========== 4. normalize + reduce ==========
logi("SCTransform…")
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)

logi("PCA/UMAP/Clustering…")
obj <- RunPCA(obj, npcs = 30, verbose = FALSE) |>
       RunUMAP(dims = 1:20) |>
       FindNeighbors(dims = 1:20) |>
       FindClusters(resolution = 0.6)

## ========== 5. plots ==========
OCP_pos <- c("CSF1R","IRF8"); OC_pos <- c("ACP5","CTSK","ATP6V0D2")
p <- (FeaturePlot(obj, OCP_pos, combine = FALSE) %>% wrap_plots()) /
     (FeaturePlot(obj, OC_pos,  combine = FALSE) %>% wrap_plots())
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plot_dir, "UMAP_markers.png"), p, width = 9, height = 6, dpi = 300)

## ========== 6. subsampling ARI ==========
logi("Technical half-subsampling (5x) for ARI…")
ari_vec <- replicate(5, {
  cells_half <- sample(colnames(obj), floor(0.5 * ncol(obj)))
  obj_half <- subset(obj, cells = cells_half) |>
              FindNeighbors(dims = 1:20) |>
              FindClusters(resolution = 0.6)
  adjustedRandIndex(Idents(obj)[cells_half], Idents(obj_half))
})
write.csv(data.frame(ARI = ari_vec), file.path(out_dir, "tech_subsample_ARI.csv"), row.names = FALSE)

## ========== 7. optional LODO ==========
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

## ========== 8. summary & session ==========
summary <- list(
  n_matrices   = length(mat_files),
  cells_before = n0,
  cells_after  = ncol(obj),
  n_genes      = nrow(obj),
  n_clusters   = length(levels(Idents(obj))),
  ari_mean     = if (exists("ari_vec")) mean(ari_vec) else NA_real_
)
write_json(summary, file.path(out_dir, "_summary.json"), auto_unbox = TRUE, pretty = TRUE)

saveRDS(obj, file.path(out_dir, "combined_seurat.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

logi("Done. Outputs at: ", out_dir)


