###############################################################################
# OCP/OC Quick‑Check WITHOUT Seurat (Matrix + uwot + irlba + ggplot2)
# 2025-07-30
###############################################################################

## 0. 载入 & 安装依赖包 ##
pkgs <- c("Matrix","readr","dplyr","ggplot2","uwot","irlba","matrixStats")
for(p in pkgs){
  if(!requireNamespace(p,quietly=TRUE)){
    install.packages(p, repos="https://cloud.r-project.org")
  }
  library(p, character.only=TRUE)
}

## 1. 路径设置 ##
root_dir <- "D:/OP GEO collection"   # ← 修改为你的路径
plot_dir <- file.path(root_dir, "QC_plots")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

## 2. 读 10X 样本辅助函数 ##
read_10x_sample <- function(sample_dir){
  mtx  <- list.files(sample_dir, pattern="matrix.*\\.mtx(\\.gz)?$", full.names=TRUE)
  bc   <- list.files(sample_dir, pattern="barcodes.*\\.tsv(\\.gz)?$", full.names=TRUE)
  feat <- list.files(sample_dir, pattern="(features|genes).*\\.tsv(\\.gz)?$", full.names=TRUE)
  if(length(mtx)!=1||length(bc)!=1||length(feat)!=1) stop("缺少 mtx/barcodes/features: ", sample_dir)
  mat   <- Matrix::readMM(mtx)
  genes <- readr::read_tsv(feat,   col_names=FALSE)$X2
  cells <- readr::read_tsv(bc,     col_names=FALSE)$X1
  rownames(mat) <- genes; colnames(mat) <- cells
  mat
}

## 3. 样本处理 + UMAP ##
process_sample <- function(mtx_file){
  sample_dir <- dirname(mtx_file)
  sample_id  <- basename(sample_dir)
  message("Processing ", sample_id)
  mat <- read_10x_sample(sample_dir)
  # QC
  nFeature  <- Matrix::colSums(mat>0)
  pct_mt    <- Matrix::colSums(mat[grep("^MT-",rownames(mat)),])/Matrix::colSums(mat)*100
  keep      <- which(nFeature>500 & pct_mt<20)
  mat       <- mat[, keep, drop=FALSE]
  if(ncol(mat)<500){ message(" ⚠ 少于500细胞，跳过"); return() }
  # Normalize->log1p(CPM)
  libsize <- Matrix::colSums(mat)
  mat_cpm  <- t(t(mat)/libsize*1e4)
  logmat   <- log1p(mat_cpm)
  # HVG
  hvg      <- names(sort(matrixStats::rowVars(as.matrix(logmat)), dec=TRUE))[1:2000]
  # PCA + UMAP
  pca <- irlba::prcomp_irlba(t(logmat[hvg,]), n=20)
  umap <- uwot::umap(pca$x, n_neighbors=30, min_dist=0.3)
  df <- data.frame(UMAP1=umap[,1], UMAP2=umap[,2], cell=rownames(pca$x))
  # Plot OCP / OC
  for(set in list(OCP=c("CSF1R","IRF8"), OC=c("ACP5","CTSK","ATP6V0D2"))){
    expr <- rowMeans(as.matrix(logmat[intersect(set,rownames(logmat)), , drop=FALSE]))
    df$expr <- expr[df$cell]
    p <- ggplot(df, aes(UMAP1,UMAP2,color=expr)) +
         geom_point(size=0.6) + theme_minimal() +
         scale_color_viridis_c() +
         ggtitle(paste(sample_id, names(set)))
    out <- file.path(plot_dir, paste0(sample_id,"_",names(set),"_UMAP.png"))
    ggsave(out, p, width=6, height=5, dpi=300)
    message(" ✓ saved ", basename(out))
  }
}

# 4. 执行
datasets <- c("GSE169396","GSE120221","GSE185381")
all_mtx  <- unlist(lapply(datasets, function(ds)
  list.files(file.path(root_dir, ds), pattern="matrix.*mtx", full.names=TRUE)))
lapply(all_mtx, process_sample)
message("🎉 完成！PNG in ", plot_dir)
