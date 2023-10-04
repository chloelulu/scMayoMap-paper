lapply(c("dplyr","Seurat","HGNChelper",'tibble','scMayoMap'), library, character.only = T)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
## --------use scType's database --------
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
gs_list = gene_sets_prepare(db_, tissue)
# signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="blood",c(2,grep('^blood',colnames(scMayoMapDatabase)))]
# rownames(signatures) <- NULL
# signatures <- signatures %>% column_to_rownames('gene')
# colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
# gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))
setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/PBMCbench/')
files <- list.files()
scType.res <- list()
scType.time <- c()
for(file in files){
  seurat.markers <- seurat_obj <- NULL
  load(file)
  rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$RNA@counts)
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  scRNAseqData <- seurat_obj[["RNA"]]@scale.data
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  
  if(sum(sctype_scores$ncells) ==0){
    cat(file,'\n')
  }else{
    scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
      full_join(seurat_obj@meta.data %>% select(celltype, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
      dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = celltype, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
  }
}

scType.res1 <- do.call(rbind,scType.res)
save(scType.res1, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType.res1_scTypeDB.RData')




setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/GB')
tissue = "Pancreas" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
gs_list = gene_sets_prepare(db_, tissue)
# signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="pancreas",c(2,grep('^pancreas',colnames(scMayoMapDatabase)))]
# rownames(signatures) <- NULL
# signatures <- signatures %>% column_to_rownames('gene')
# colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
# gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))

files <- list.files()
scType.time <- c()
scType.res <- list()
for(file in files){
  seurat.markers <- seurat_obj <- NULL
  load(file)
  # rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$RNA@counts)
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  scRNAseqData <- seurat_obj[["RNA"]]@scale.data
  rownames(scRNAseqData) <- gsub('\\-\\-.*','',rownames(scRNAseqData))
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  if(sum(sctype_scores$ncells) ==0){
    cat(file,'\n')
  }else{
    scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
      full_join(seurat_obj@meta.data %>% select(celltype, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
      dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = celltype, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
  }
}

scType.res2 <- do.call(rbind,scType.res)
save(scType.res2, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType.res2_scTypeDB.RData')













setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/TabulaMuris/')
# signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="pancreas",c(2,grep('^pancreas',colnames(scMayoMapDatabase)))]
# rownames(signatures) <- NULL
# signatures <- signatures %>% column_to_rownames('gene')
# colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
# gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))
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
# files1 <- files1[which(names(files1) %in% scType.DB)]
scType.res <- list()
scType.time <- c()
for(file in unname(files1)){
  cat(file,'\n')
  if(file %in% files3){
    tissue = scType.DB
    }else{
      tissue = names(files1[files1 %in% file])
  }
   # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
  gs_list = gene_sets_prepare(db_, tissue)
  
  
  seurat.markers <- seurat_obj <- NULL
  load(file)
  # rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$originalexp@counts)
  if(length(grep('^cluster$',colnames(seurat_obj@meta.data)))==0){
    seurat_obj@meta.data$cluster <- seurat_obj@meta.data$cluster.cell 
  }
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj@assays$originalexp@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  gc()
  scRNAseqData <- seurat_obj@assays$originalexp@scale.data
  rownames(scRNAseqData) <- gsub('\\-\\-.*','',rownames(scRNAseqData))
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ]), drop =F]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  
  scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
    full_join(seurat_obj@meta.data %>% select(cell_ontology_class, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
    dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = cell_ontology_class, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
}

length(scType.res)
scType.res3 <- do.call(rbind,scType.res)
save(scType.res3, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType.res3_scTypeDB.RData')




scType_res <- rbind(as.data.frame(scType.res1), as.data.frame(scType.res2),as.data.frame(scType.res3))
save(scType_res, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType_res.RData')
write.csv(scType_res, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType_res.csv', row.names = F)



## --------use scMayoMap's database --------
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
gs_list = gene_sets_prepare(db_, tissue)
signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="blood",c(2,grep('^blood',colnames(scMayoMapDatabase)))]
rownames(signatures) <- NULL
signatures <- signatures %>% column_to_rownames('gene')
colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))
setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/PBMCbench/')
files <- list.files()
scType.res <- list()
scType.time <- c()
for(file in files){

  seurat.markers <- seurat_obj <- NULL
  load(file)
  rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$RNA@counts)
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  scRNAseqData <- seurat_obj[["RNA"]]@scale.data
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  if(sum(sctype_scores$ncells) ==0){
    cat(file,'\n')
  }else{
    scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
      full_join(seurat_obj@meta.data %>% select(celltype, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
      dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = celltype, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
  }
}

scType.res1 <- do.call(rbind,scType.res)
save(scType.res1, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType.res1_scMayoMapDB.RData')




setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/GB')
tissue = "Pancreas" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
gs_list = gene_sets_prepare(db_, tissue)
signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue=="pancreas",c(2,grep('^pancreas',colnames(scMayoMapDatabase)))]
rownames(signatures) <- NULL
signatures <- signatures %>% column_to_rownames('gene')
colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))

files <- list.files()
scType.time <- c()
scType.res <- list()
for(file in files){
  seurat.markers <- seurat_obj <- NULL
  load(file)
  # rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$RNA@counts)
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  scRNAseqData <- seurat_obj[["RNA"]]@scale.data
  rownames(scRNAseqData) <- gsub('\\-\\-.*','',rownames(scRNAseqData))
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  
  if(sum(sctype_scores$ncells) ==0){
    cat(file,'\n')
  }else{
    scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
      full_join(seurat_obj@meta.data %>% select(celltype, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
      dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = celltype, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
  }
}

scType.res2 <- do.call(rbind,scType.res)
save(scType.res2, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType.res2_scMayoMapDB.RData')













setwd('/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/data/TabulaMuris/')
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
names(files1) <- firstup(name)
scType.DB <- c('Immune system','Pancreas','Liver','Eye','Kidney','Brain','Lung','Adrenal','Heart','Intestine','Muscle','Placenta','Spleen','Stomach','Thymus')
files2 <- names(files1)[!(names(files1) %in% scType.DB)]
files3 <- files[which(!(names(files1) %in% scType.DB))]
files <- files[which(names(files1) %in% scType.DB)]
# files1 <- files1[which(names(files1) %in% scType.DB)]
scType.res <- list()
scType.time <- c()
for(file in unname(files1)[29:30]){
  cat(file,'\n')
  if(file %in% files3){
    tissue = scType.DB
  }else{
    tissue = names(files1[files1 %in% file])
  }
  # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
  gs_list = gene_sets_prepare(db_, tissue)
  sub.tissue <- tolower(names(files1[files1 %in% file]))
  if(length(sub.tissue) ==1){
    signatures <- scMayoMapDatabase[scMayoMapDatabase$tissue==sub.tissue,c(2,grep(paste0('^', sub.tissue), colnames(scMayoMapDatabase)))]
    rownames(signatures) <- NULL
    signatures <- signatures %>% column_to_rownames('gene')
    colnames(signatures) <- gsub('.*\\:','',colnames(signatures))
    gs_list$gs_positive <- apply(signatures, 2, function(x) names(which(x!=0)))
  }else{
    signatures <- scMayoMapDatabase[,-1]
    rownames(signatures) <- NULL
    gs_list$gs_positive <- apply(signatures[,-1], 2, function(x) {
      signatures$gene[which(x!=0)]
    })
  }
  seurat.markers <- seurat_obj <- NULL
  load(file)
  # rownames(seurat_obj@meta.data) <- colnames(seurat_obj@assays$originalexp@counts)
  if(length(grep('^cluster$',colnames(seurat_obj@meta.data)))==0){
    seurat_obj@meta.data$cluster <- seurat_obj@meta.data$cluster.cell 
  }
  seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj@assays$originalexp@counts, na.rm = TRUE)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  gc()
  scRNAseqData <- seurat_obj@assays$originalexp@scale.data
  rownames(scRNAseqData) <- gsub('\\-\\-.*','',rownames(scRNAseqData))
  t1 <- Sys.time()
  es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$cluster), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$cluster==cl, ]), drop =F]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$cluster==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  t2 <- Sys.time()
  scType.time[file] <- difftime(t2, t1, units = 'secs')
  
  scType.res[[file]] <- sctype_scores %>% mutate(dataset = file) %>% 
    full_join(seurat_obj@meta.data %>% select(cell_ontology_class, cluster) %>% distinct() %>% tibble::as_tibble()) %>% 
    dplyr::select(-c('scores','ncells')) %>% dplyr::rename(benchmark = cell_ontology_class, scType = type) %>% dplyr::select(dataset, cluster, benchmark, scType)
}

length(scType.res)
scType.res3 <- do.call(rbind,scType.res)
save(scType.res3, scType.time, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType//scType.res3_scMayoMapDB.RData')



scType_res <- rbind(as.data.frame(scType.res1), as.data.frame(scType.res2),as.data.frame(scType.res3))
save(scType_res, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType_res_scMayoMapDB.RData')
write.csv(scType_res, file = '/Users/M216453/Documents/Mayo_project/2022_01_12_scMapping/result/scType/scType_res_scMayoMapDB.csv', row.names = F)






