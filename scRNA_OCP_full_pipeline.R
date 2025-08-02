###############################################################################
#  scRNA_OCP_full_pipeline.R  --  CI-friendly with diagnostics
#  Run:
#    Rscript scRNA_OCP_full_pipeline.R --data_dir data --out_dir docs
###############################################################################

## ========== 0. CLI options ==========
quiet_install <- function(pkgs) {
  need <- pkgs[!sapply(pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")
}
quiet_install(c("optparse"))
suppressPackageStartupMessages(library(optparse))

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

## ========== helpers ==========
die <- function(...) {
  msg <- paste0(format(Sys.time(), "%F %T"), " ERROR: ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, file = file.path(out_dir, "_error.txt"), sep = "\n")
  # also dump sessionInfo
  sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
  quit(status = 1, save = "no")
}

logi <- function(...) message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ...))

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
if (length(mat_files) == 0) {
  # print a directory snapshot to help debugging
  cat("No matrix files found. dir snapshot:\n", file = file.path(out_dir, "_error.txt"))
  sink(file.path(out_dir, "_error.txt"), append = TRUE)
  print(list.dirs(root_dir, full.names = TRUE, recursive = FALSE))
  print(head(list.files(root_dir, recursive = TRUE, full.names = TRUE), 200))
  sink()
  die("No _matrix.mtx.gz files found under ", root_dir,
      " — check LFS checkout or DATA_DIR path.")
}
logi("Found ", length(mat_files), " matrices.")

read_one <- function(mtx_fp) {
  dirp <- dirname(mtx_fp)
  prefix <- sub("_matrix.*\\.mtx\\.gz$", "", basename(mtx_fp))
  pick1 <- function(pat) {
    fs <- list.files(dirp, pat, full.names = TRUE)
    if (length(fs)) fs[order(!grepl("\\.gz$", fs))][1] else NA_character_
  }
  bc   <- pick1(paste0("^", prefix, ".*barcode[s]*.*tsv(\\.gz)?$"))
  feat <- pick1(paste0("^", prefix, ".*(features|genes).*tsv(\\.gz)?$"))
  if (is.na(bc) || is.na(feat))
    die("Missing barcodes/features for prefix: ", prefix, " in ", dirp)
  CreateSeuratObject(ReadMtx(mtx = mtx_fp, cells = bc, features = feat),
                     project = basename(dirname(mtx_fp)))
}

## ========== 3. load + QC ==========
logi("Loading matrices into Seurat objects…")
objs <- lapply(mat_files, read_one)
obj  <- if (length(objs) == 1) objs[[1]] else Reduce(function(x,y) merge(x,y), objs)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^MT-")
n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA > 500 & percent.mt < 20)
logi("QC filter: cells kept ", ncol(obj), "/", n0)
if (ncol(obj) < 500) die("Too few cells after QC (", ncol(obj), ").")

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
ggsave(file.path(plot_dir, "UMAP_markers.png"), p, width = 9, height = 6, dpi = 300)

## ========== 6. technical subsampling ==========
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
  logi("LODO…")
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

## ========== 8. DE (if labels provided) ==========
if ("short_cluster" %in% colnames(obj@meta.data)) {
  Idents(obj) <- "short_cluster"
  if (all(c("OCP","OC") %in% levels(Idents(obj)))) {
    logi("FindMarkers OCP vs OC (wilcox)…")
    de <- FindMarkers(obj, ident.1 = "OCP", ident.2 = "OC", test.use = "wilcox")
    write.csv(de, file.path(out_dir, "DE_OCP_vs_OC.csv"))
    if (nrow(de) > 0) {
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
  } else {
    logi("short_cluster exists but lacks OCP/OC levels; skip DE.")
  }
} else {
  logi("No 'short_cluster' in meta; skip DE.")
}

## ========== 9. save and session info ==========
saveRDS(obj, file.path(out_dir, "combined_seurat.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
logi("Done. Outputs at: ", out_dir)

