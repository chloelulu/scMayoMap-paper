library(scCATCH);library(dplyr);library(tidyr); library(scMayoMap)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
## === scCATCH with internal DB ====
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
  cat('----- ',file,' -----\n')
  tissue <- names(files1[files %in% file])
  t1 <- Sys.time()
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  rownames(dt) <- toupper(rownames(dt))
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = "Human", geneinfo = geneinfo) # We are not able to detect markers by Mouse database, thus we try human
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cluster.cell))
  scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = "Human", marker = cellmatch, tissue = "Blood")
  scCATCH_obj <- findcelltype(object = scCATCH_obj)
  t2 <- Sys.time()
  scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
  save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH/',file, '_scCATCH.Rdata'))
}


files <- ff
ff2 <- c()
for(file in ff){
  cat('----- ',file,' -----\n')
  load(file)
  meta.dat <- seurat_obj@meta.data %>%
    dplyr::select(cluster.cell, cell_ontology_class) %>%
    dplyr::rename(celltype = cell_ontology_class)
  rownames(meta.dat) <- NULL
  gc()
  tissue <- names(files1[files1 %in% file])
  cat('-- scCATCH --\n')
  t1 <- Sys.time()
  dt <- scCATCH_obj <- NULL
  dt <- as.matrix(seurat_obj@assays$originalexp@counts)
  rownames(dt) <- toupper(rownames(dt))
  species <- "Human"
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = species, geneinfo = geneinfo) # We are not able to detect markers by Mouse database, thus we try human
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cell_ontology_class))
  try(scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = species, marker = cellmatch, tissue = tissue))
  if(nrow(scCATCH_obj@marker)>0){
    scCATCH_obj <- findcelltype(object = scCATCH_obj)
    t2 <- Sys.time()
    scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
    scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
    save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH/',file, '_',species,'_scCATCH.Rdata'))
  }else{
    ff2 <- c(ff2, file)
  }
}




setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/data/GenomeBiology2019_Tamin/cleaned_Lu/')
files <- c("Baron Human.RData","Baron Mouse.RData","Muraro.RData","Segerstolpe.RData", "Xin.RData")
for(i in 1:5){
  file <- files[i]
  cat(file,'\n')
  load(file)
  gc()
  tissue <- "Pancreas"
  # cat('-- scCATCH --\n')
  t1 <- Sys.time()
  dt <- scCATCH_obj <- NULL
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  if(file =="Muraro.RData"){
    rownames(dt) <- toupper(gsub('--.*','',rownames(dt)))
  }
  if(length(grep('Mouse',file))==1){
    species <- "Mouse"
  }else{
    species <- "Human"
  }
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = species, geneinfo = geneinfo) # We are not able to detect markers by Mouse database, thus we try human
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$celltype))
  try(scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = species, marker = cellmatch, tissue = tissue))
  if(nrow(scCATCH_obj@marker)>0){
    scCATCH_obj <- findcelltype(object = scCATCH_obj)
    t2 <- Sys.time()
    scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
    scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
    save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH/',file, '_',species,'_scCATCH.Rdata'))
  }
}





setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/')
files <- list.files(pattern = 'Rdata')
for(file in files){
  cat('[ ',file,' ]\n')
  seurat_obj <- NULL
  load(file)
  gc()
  cat('-- scCATCH --\n')
  t1 <- Sys.time()
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  rownames(dt) <- toupper(rownames(dt))
  species = "Human"
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = species, geneinfo = geneinfo) # We are not able to detect markers by Mouse database, thus we try human
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cluster))
  scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = species, marker = cellmatch, tissue = "Blood")
  scCATCH_obj <- findcelltype(object = scCATCH_obj)
  t2 <- Sys.time()
  scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
  save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH/',file, '_',species,'_scCATCH.Rdata'))
}





## ==== scCATCH with scMayoDB ====
## create scMayoDB to scCATCH database format [ref:https://github.com/ZJUFanLab/scCATCH/wiki/Add-more-marker-genes-to-cellmatch-for-annotation]
custom_marker <- NULL
for(i in unique(scMayoMapDatabase$tissue)){
  tmp <- scMayoMapDatabase[scMayoMapDatabase$tissue ==i,-1] %>% `rownames<-`( NULL )
  colnames(tmp) <- gsub('.*:','',colnames(tmp))
  tmp0 <- tmp %>% pivot_longer(cols = -gene, names_to = "celltype", values_to = "Presence") %>%  
    dplyr::filter(Presence == 1) %>% dplyr::select(-Presence) %>% 
    mutate(species = 'Human',tissue = i,cancer = 'Normal',condition = "Normal cell",subtype1 = NA,subtype2 = NA,subtype3 = NA,resource = NA,pmid = NA) %>% 
    dplyr::select(c('species','tissue','cancer','condition','subtype1','subtype2','subtype3','celltype','gene','resource','pmid'))
  custom_marker <- rbind(custom_marker, tmp0)
}
custom_marker_all <- custom_marker
custom_marker_all$tissue <- 'unknown'

setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/')
files <- list.files(pattern = 'Rdata')
for(file in files){
  cat('[ ',file,' ]\n')
  seurat_obj <- NULL
  load(file)
  gc()
  cat('-- scCATCH --\n')
  t1 <- Sys.time()
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  rownames(dt) <- toupper(rownames(dt))
  species = "Human"
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = species, geneinfo = geneinfo) 
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cluster))
  scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = species, marker = custom_marker, tissue = "blood", if_use_custom_marker=T)
  scCATCH_obj <- findcelltype(object = scCATCH_obj)
  t2 <- Sys.time()
  scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type') %>% mutate(dataset = gsub('.Rdata|.RData','',file))
  save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH_scMayoDB/',file, '_',species,'_scCATCH.Rdata'))
}






setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/data/GenomeBiology2019_Tamin/cleaned_Lu/')
files <- c("Baron Human.RData","Baron Mouse.RData","Muraro.RData","Segerstolpe.RData", "Xin.RData")
tissue <- "pancreas"
for(file in files){
  cat(file,'\n')
  load(file)
  gc()
  t1 <- Sys.time()
  dt <- scCATCH_obj <- NULL
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  dt <- dt[rowSums(dt)>0,]
  if(file =="Muraro.RData"){rownames(dt) <- toupper(gsub('--.*','',rownames(dt)))}
  expr.mtx <- rev_gene(data = dt, data_type = "data", species = 'Human', geneinfo = geneinfo)
  scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$celltype))
  try(scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = species, marker = custom_marker, tissue = tissue, if_use_custom_marker=T))
  if(nrow(scCATCH_obj@marker)>0){
    scCATCH_obj <- findcelltype(object = scCATCH_obj)
    t2 <- Sys.time()
    scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
    scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
    save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH_scMayoDB/',file, '_',species,'_scCATCH.Rdata'))
  }
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
name[grep("large_intestine",name)] = "gastrointestinal tract"
name[grep("marrow",name)] = "bone marrow"
name[grep("mammary_gland",name)] = "mammary gland"
name[grep("trachea",name)] = "lung"
names(files1) <- (name)
ff <- c()
for(file in files){
  cat(file,'\n')
  load(file)
  gc()
  if(file %in% c("Tongue_smartseq2.RData","Tongue_droplet.RData")){
    tissue <- 'unknown'
  }else{
    tissue <- names(files1[files %in% file])
  }
  
  t1 <- Sys.time()
  dt <- as.matrix(seurat_obj@assays$originalexp@counts)
  rownames(dt) <- toupper(rownames(dt))
  dt <- dt[rowSums(dt)>0,]
  # if(file =="Kidney_droplet.RData"){
  #   dt <- dt[,which(seurat_obj$cluster.cell !='1')] # only contains 1 cell in this cluster
  #   expr.mtx <- rev_gene(data = dt, data_type = "data", species = "Human", geneinfo = geneinfo)
  #   scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cluster.cell[which(seurat_obj$cluster.cell !='1')]))
  # }else{
    expr.mtx <- rev_gene(data = dt, data_type = "data", species = "Human", geneinfo = geneinfo)
    scCATCH_obj <- createscCATCH(data = expr.mtx, cluster=paste0('Cluster ',seurat_obj$cluster.cell))
  # }
  
  if(file %in% c("Tongue_smartseq2.RData","Tongue_droplet.RData")){
    scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = "Human", marker = custom_marker_all, tissue = tissue, if_use_custom_marker=T)
  }else{
    scCATCH_obj <- findmarkergene(object = scCATCH_obj, species = "Human", marker = custom_marker, tissue = tissue, if_use_custom_marker=T)
  }
  
  scCATCH_obj <- findcelltype(object = scCATCH_obj)
  t2 <- Sys.time()
  scCATCH.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  scCATCH_res <- scCATCH_obj@celltype[,c('cluster','cell_type')] %>% dplyr::rename(scCATCH = 'cell_type')
  save(scCATCH_obj, scCATCH_res, scCATCH.time, file =paste0('/research/bsi/projects/staff_analysis/m216453/scmapping/tmp/scCATCH_scMayoDB/',file, '_scCATCH.Rdata'))
}














