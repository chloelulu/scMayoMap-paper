setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/')
files <- list.files(pattern = 'Rdata')
library(SCINA)
library(scMayoMap)
library(tibble)
library(dplyr)
signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="blood",c(2,grep('^blood',colnames(scMayoMapDatabase)))]
rownames(signatures) <- NULL
signatures <- signatures %>% column_to_rownames('gene')
colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
signatures[1:4,1:4]
signatures <- apply(signatures, 2, function(x) names(which(x!=0)))
SCINA.time <- c()
names(SCINA.time) <- files
res.list <- list()
for(file in files){
  cat('[ ',file,' ]\n')
  seurat.markers <- seurat_obj <- NULL
  load(file)
  cat('-- SCINA --\n')
  exp <- seurat_obj@assays$RNA@counts
  t1 <- Sys.time()
  results <- SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                   convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                   rm_overlap=F, allow_unknown=TRUE, log_file='/research/bsi/projects/staff_analysis/m216453/scmapping/SCINA.log')
  t2 <- Sys.time()
  SCINA.time[file] <- difftime(t2,t1, units = 'secs')
  SCINA.res <- cbind.data.frame(seurat_obj@meta.data %>% dplyr::select(celltype, cluster),SCINA =results$cell_labels)
  res.list[[file]] <- SCINA.res
}

SCINA.res <- do.call(rbind, res.list)
write.csv(SCINA.res, file = '/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/SCINA_res.csv')







setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/')
files <- c("Xin.RData","Baron Human.RData","Baron Mouse.RData","Segerstolpe.RData","Muraro.RData")
signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="pancreas",c(2,grep('^pancreas',colnames(scMayoMapDatabase)))]
rownames(signatures) <- NULL
signatures <- signatures %>% column_to_rownames('gene')
colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
signatures[1:4,1:4]
signatures <- apply(signatures, 2, function(x) names(which(x!=0)))

SCINA.time <- c()
res.list <- list()
for(file in files[3:5]){
  cat('[ ',file,' ]\n')
  seurat.markers <- seurat_obj <- NULL
  load(file)
  cat('-- SCINA --\n')
  
  exp <- seurat_obj@assays$RNA@counts
  results <- NULL
  t1 <- Sys.time()
  results <- SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                   convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                   rm_overlap=F, allow_unknown=TRUE, log_file='/research/bsi/projects/staff_analysis/m216453/scmapping/SCINA.log')
  t2 <- Sys.time()
  SCINA.time[file] <- difftime(t2,t1, units = 'secs')
  SCINA.res <- cbind.data.frame(seurat_obj@meta.data %>% dplyr::select(celltype, cluster),SCINA =results$cell_labels)
  res.list[[file]] <- SCINA.res
}

SCINA.res <- do.call(rbind, res.list)
write.csv(SCINA.res, file = '/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/SCINA_res.csv')








setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/data/TabulaMuris/')
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
name[grep("Intestine",name)] = "gastrointestinal tract"
name[grep("tongue",name)] = "other"
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
names(files1) <- firstup(name)
scType.DB <- c('Immune system','Pancreas','Liver','Eye','Kidney','Brain','Lung','Adrenal','Heart','Intestine','Muscle','Placenta','Spleen','Stomach','Thymus')
files2 <- names(files1)[!(names(files1) %in% scType.DB)]
files3 <- files[which(!(names(files1) %in% scType.DB))]

files <- files[which(names(files1) %in% scType.DB)]
files1 <- files1[which(names(files1) %in% scType.DB)]


SCINA.time <- c()
res.list <- list()
for(file in files){
  cat('[ ',file,' ]\n')
  seurat.markers <- seurat_obj <- NULL
  load(file)
  cat('-- SCINA --\n')
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  if(file %in% files3){
    tissue = scType.DB
  }else{
    tissue = names(files1[files1 %in% file])
  }
  # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
  gs_list = gene_sets_prepare(db_, tissue)
  
  tissue.sub <- names(files1[which(file %in% files1)])
  signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue==tissue.sub,c(2,grep(tissue.sub,colnames(scMayoMapDatabase)))]
  rownames(signatures) <- NULL
  signatures <- signatures %>% column_to_rownames('gene')
  colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
  signatures[1:4,1:4]
  signatures <- apply(signatures, 2, function(x) names(which(x!=0)))
  
  exp <- seurat_obj@assays$originalexp@data[1:4,1:4]
  results <- SCINA.res<- NULL
  
  try({
    t1 <- Sys.time()
    results <- SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                   convergence_rate = 0.999, sensitivity_cutoff = 0.9, 
                   rm_overlap=F, allow_unknown=TRUE, log_file='/research/bsi/projects/staff_analysis/m216453/scmapping/SCINA.log')
    t2 <- Sys.time()
    SCINA.time[file] <- difftime(t2,t1, units = 'secs')
    SCINA.res <- cbind.data.frame(seurat_obj@meta.data %>% dplyr::select(celltype, cluster),SCINA =results$cell_labels)
    })
  
  res.list[[file]] <- SCINA.res
}

SCINA.res <- do.call(rbind, res.list)
SCINA.res$dataset <- gsub('.RData.*','',rownames(SCINA.res))
table(SCINA.res$dataset)
write.csv(SCINA.res, file = '/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/SCINA_res.csv')


