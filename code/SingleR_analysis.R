library(SingleR);library(dplyr);library(tidyr)
hpca.se <- HumanPrimaryCellAtlasData()
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


setwd("/research/bsi/projects/staff_analysis/m216453/scmapping/data/TabulaMuris/")
files <- list.files(pattern = 'RData$')
files1 <- files
name <- tolower(((gsub('_droplet.RData|_smartseq2.RData','',files))))
name[grep('fat',name)] = "adipose tissue"
name[grep("limb_muscle",name)] = "muscle"
name[grep("brain_myeloid",name)] = "brain"
name[grep("brain_non-myeloid",name)] = "brain"
name[grep("heart_and_aorta",name)] = "heart"
name[grep("large_intestine",name)] = "Intestine"
name[grep("marrow",name)] = "bone marrow"
name[grep("mammary_gland",name)] = "mammary gland"
name[grep("trachea",name)] = "lung"
name[grep("tongue",name)] = 'Tonsil'
names(files1) <- firstup(name)
tissues <- unique(cellmatch$tissue)
ff <- c()
for(file in files){
  cat(file,'\n')
  load(file)
  meta.dat <- seurat_obj@meta.data %>% dplyr::select(cluster.cell, cell_ontology_class) %>% dplyr::rename(celltype = cell_ontology_class)
  rownames(meta.dat) <- NULL
  gc()
  cat('-- SingleR --\n')
  t1 <- Sys.time()
  pred.hesc <- SingleR(test =dt, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  t2 <- Sys.time()
  SingleR.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  SingleR_res <- cbind.data.frame(SingleR = pred.hesc$labels, meta.dat %>% dplyr::rename(benchmark = 'celltype'))
  save(SingleR_res, pred.hesc, SingleR.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/SingleR/',file, '_SingleR.Rdata'))
}




setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/data/GenomeBiology2019_Tamin/cleaned_Lu/')
files <- c("Baron Human.RData","Baron Mouse.RData","Muraro.RData","Segerstolpe.RData", "Xin.RData")
for(i in 1:5){
  file <- files[i]
  cat(file,'\n')
  load(file)
  meta.dat <- seurat_obj@meta.data %>% 
    dplyr::select(cluster, celltype) 
  rownames(meta.dat) <- NULL
  gc()
  cat('-- SingleR --\n')
  t1 <- Sys.time()
  pred.hesc <- SingleR(test =dt, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  t2 <- Sys.time()
  SingleR.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  SingleR_res <- cbind.data.frame(SingleR = pred.hesc$labels, meta.dat %>% dplyr::rename(benchmark = 'celltype'))
  save(SingleR_res, pred.hesc, SingleR.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/SingleR/',file, '_SingleR.Rdata'))
  }






setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/')
files <- list.files(pattern = 'Rdata')
for(file in files){
  cat('[ ',file,' ]\n')
  seurat_obj <- NULL
  load(file)
  cat('[ ',file,' ]\n')
  seurat_obj <- meta.dat <- NULL
  load(file)
  meta.dat <- seurat_obj@meta.data %>% dplyr::select(celltype, cluster)
  gc()
  cat('-- SingleR --\n')
  t1 <- Sys.time()
  pred.hesc <- SingleR(test =dt, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
  t2 <- Sys.time()
  SingleR.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  SingleR_res <- cbind.data.frame(SingleR = pred.hesc$labels, meta.dat %>% dplyr::rename(benchmark = 'celltype'))
  save(SingleR_res, pred.hesc, SingleR.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/SingleR/',file, '_SingleR.Rdata'))
  }

