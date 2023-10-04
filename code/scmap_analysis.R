library(dplyr);library(tibble);library(tidyr)
library(SingleCellExperiment);library(scmap) #(http://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html) Here we only test mapCluster 
library(Seurat) 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

ref <- as.SingleCellExperiment(pbmcsca)
rowData(ref)$feature_symbol <- rownames(ref)
ref <- selectFeatures(ref, suppress_plot = T)
ref <- indexCluster(ref,cluster_col = "CellType")


files <- list.files(pattern = 'Rdata')
for(file in files){
  cat('[ ',file,' ]\n')
  seurat_obj <- NULL
  load(file)
  cat('-- scmap --\n')
  cat(nrow(seurat_obj@assays$RNA@counts),'\n')
  t1 <- Sys.time()
  sce <- NormalizeData(object = seurat_obj)
  sce <- as.SingleCellExperiment(sce)
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  sce <- sce[,colData(sce)$celltype != "Unassigned"]
  sce <- sce[,colData(sce)$celltype != "unassigned"]
  sce <- selectFeatures(sce, suppress_plot = F)
  scmapCluster_results <- scmapCluster(#https://github.com/hemberg-lab/scmap/issues/10
    projection = sce, 
    index_list = list(ref = metadata(ref)$scmap_cluster_index))
  t2 <- Sys.time()
  scmap.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  scmap_res <- cbind.data.frame(scmap = scmapCluster_results$combined_labs, seurat_obj@meta.data[,c('celltype', 'cluster')] %>%
                                  dplyr::rename(benchmark = 'celltype'))
  save(scmap_res, scmapCluster_results, scmap.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/SingleR/scmap/',file, '_scmap.Rdata'))
}
