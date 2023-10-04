suppressPackageStartupMessages({
  library(ExperimentHub)
  library(SingleCellExperiment)
  library(TabulaMurisData)
})
droplet <- TabulaMurisDroplet()
smartseq2 <- TabulaMurisSmartSeq2()

unique(droplet@colData$tissue)

## Split tissue object
for(tissue in unique(droplet@colData$tissue)){
  cat(tissue,'\n')
  sub.droplet <- subset(droplet, , tissue==tissue)
  sub.droplet.seurat <- as.Seurat(sub.droplet, counts = "counts", data = "counts")
  Idents(sub.droplet.seurat) <- "cell_ontology_class"
  save(sub.droplet.seurat, file = paste0('data/',tissue,'_droplet.RData'))
}


for(i in unique(smartseq2@colData$tissue)){
  cat(tissue,'\n')
  sub.smartseq2 <- subset(smartseq2, , tissue==i)
  sub.smartseq2.seurat <- as.Seurat(sub.smartseq2, counts = "counts", data = "counts")
  Idents(sub.smartseq2.seurat) <- "cell_ontology_class"
  save(sub.smartseq2.seurat, file = paste0('data/',tissue,'_smartseq2.RData'))
}



smartseq2 <- TabulaMurisSmartSeq2()
for(i in unique(smartseq2@colData$tissue)){
  sub.smartseq2 <- subset(smartseq2, , tissue== i)
  sub.smartseq2.seurat <- as.Seurat(sub.smartseq2, counts = "counts", data = "counts")
  Idents(sub.smartseq2.seurat) <- seurat_obj@meta.data$cell_ontology_class
  seurat_obj <- sub.smartseq2.seurat
  cat(length(unique(seurat_obj@meta.data$cell_ontology_class)),' celltypes \n')
  seurat_obj <- NormalizeData(object = seurat_obj)
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj)
  seurat_obj <- FindNeighbors(object = seurat_obj)
  seurat_obj@meta.data$cell_ontology_class[is.na(seurat_obj@meta.data$cell_ontology_class)] <- 'Unclassified'
  seurat_obj@meta.data$cluster.cell <- as.character(as.numeric(as.factor(seurat_obj@meta.data$cell_ontology_class)))
  cat(length(unique(seurat_obj@meta.data$cluster.cell)),' clusters \n')
  seurat_obj <- RunTSNE(object = seurat_obj)
  Idents(seurat_obj) <- seurat_obj@meta.data$cluster.cell
  seurat.markers <- FindAllMarkers(seurat_obj, method = 'MAST')
  save(seurat.markers, seurat_obj, file = paste0('data/TabulaMuris/',i,'_smartseq2.RData'))
}



droplet <- TabulaMurisDroplet()
for(i in unique(droplet@colData$tissue)){
  sub.droplet <- subset(droplet, , tissue== i)
  sub.droplet.seurat <- as.Seurat(sub.droplet, counts = "counts", data = "counts")
  Idents(sub.droplet.seurat) <- seurat_obj@meta.data$cell_ontology_class
  seurat_obj <- sub.droplet.seurat
  cat(length(unique(seurat_obj@meta.data$cell_ontology_class)),' celltypes \n')
  seurat_obj <- NormalizeData(object = seurat_obj)
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj)
  seurat_obj <- FindNeighbors(object = seurat_obj)
  seurat_obj@meta.data$cell_ontology_class[is.na(seurat_obj@meta.data$cell_ontology_class)] <- 'Unclassified'
  seurat_obj@meta.data$cluster.cell <- as.character(as.numeric(as.factor(seurat_obj@meta.data$cell_ontology_class)))
  cat(length(unique(seurat_obj@meta.data$cluster.cell)),' clusters \n')
  seurat_obj <- RunTSNE(object = seurat_obj)
  Idents(seurat_obj) <- seurat_obj@meta.data$cluster.cell
  seurat.markers <- FindAllMarkers(seurat_obj, method = 'MAST')
  save(seurat.markers, seurat_obj, file = paste0('data/TabulaMuris/',i,'_droplet.RData'))
}



#### -------------- Arrange data from GenomeBiology2019_Tamin into Seurat format ------ ####
folder <- list.files("data/GenomeBiology2019_Tamin/Intra-dataset/Pancreatic_data/")
for(i in folder){
  cat('------',i,'----- \n')
  fs <- list.files(paste0("data/GenomeBiology2019_Tamin/Intra-dataset/Pancreatic_data/",i))
  expr <- t(read.csv(paste0("data/GenomeBiology2019_Tamin/Intra-dataset/Pancreatic_data/",i,'/',fs[1]), row.names = 1))
  cell.meta <- read.csv(paste0("data/GenomeBiology2019_Tamin/Intra-dataset/Pancreatic_data/",i,'/',fs[2]))
  seurat_obj <- CreateSeuratObject(counts = expr)
  seurat_obj@meta.data$celltype <- cell.meta$x
  save(seurat_obj, file = paste0("data/GenomeBiology2019_Tamin/cleaned_Lu/",i,'.RData'))
}


data <- read.csv('10Xv2/10Xv2_pbmc1.csv') %>% column_to_rownames('X')
labels <- read.csv("10Xv2/10Xv2_pbmc1Labels.csv")
a1 <- (grep('^pbmc1_10X_v2',rownames(data)))
a2 <- (grep('^pbmc1_10X_v3',rownames(data)))
a3 <- (grep('^pbmc1_Celseq2',rownames(data)))
a4 <- (grep('^pbmc1_Drop',rownames(data)))
a5 <- (grep('^pbmc1_inDrops',rownames(data)))
a6 <- (grep('^pbmc1_SM2',rownames(data)))
a7 <- (grep('^pbmc1_Seqwell',rownames(data)))
a11 <- (grep('^pbmc2_10X_V2',rownames(data)))

pbmc1_10X_v2 <- data[a1,]
meta.dat <- labels[a1,,drop =F]
save(pbmc1_10X_v2, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_10X_v2.Rdata')

pbmc1_10X_v3 <- data[a2,]
meta.dat <- labels[a2,,drop =F]
save(pbmc1_10X_v3, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_10X_v3.Rdata')

pbmc1_Celseq2 <- data[a3,]
meta.dat <- labels[a3,,drop =F]
save(pbmc1_Celseq2, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_Celseq2.Rdata')

pbmc1_Drop <- data[a4,]
meta.dat <- labels[a4,,drop =F]
save(pbmc1_Drop, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_Drop.Rdata')

pbmc1_inDrops <- data[a5,]
meta.dat <- labels[a5,,drop =F]
save(pbmc1_inDrops, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_inDrops.Rdata')

pbmc1_SM2 <- data[a6,]
meta.dat <- labels[a6,,drop =F]
save(pbmc1_SM2, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_SM2.Rdata')

pbmc1_Seqwell <- data[a7,]
meta.dat <- labels[a7,,drop =F]
save(pbmc1_Seqwell, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc1_Seqwell.Rdata')

pbmc2_10X_v2 <- data[a11,]
meta.dat <- labels[a11,,drop =F]
save(pbmc2_10X_v2, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc2_10X_v2.Rdata')


data <- read.csv('CEL-Seq/CL_pbmc1.csv') %>% column_to_rownames('X')
a33 <- (grep('^pbmc2_Celseq2',rownames(data)))
pbmc2_Celseq2 <- data[a33,]
labels <- read.csv("CEL-Seq/CL_pbmc1Labels.csv")
meta.dat <- labels[a33,, drop =F]
save(pbmc2_Celseq2, meta.dat,file = '../../cleaned_Lu/pbmcbench/pbmc2_Celseq2.Rdata')

data <- read.csv('Drop-Seq/DR_pbmc1.csv') %>% column_to_rownames('X')
a44 <- (grep('^pbmc2_Drop',rownames(data)))
pbmc2_Drop <- data[a44,]
labels <- read.csv("Drop-Seq/DR_pbmc1Labels.csv")
meta.dat <- labels[a44,, drop =F]
save(pbmc2_Drop, meta.dat,file = '../../cleaned_Lu/pbmcbench/pbmc2_Drop.Rdata')

data <- read.csv('inDrop/iD_pbmc1.csv') %>% column_to_rownames('X')
a55 <- (grep('^pbmc2_inDrops',rownames(data)))
pbmc2_inDrops <- data[a55,]
labels <- read.csv("inDrop/iD_pbmc1Labels.csv")
meta.dat <- labels[a55,, drop =F]
save(pbmc2_inDrops, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc2_inDrops.Rdata')

data <- read.csv('Seq-Well/SW_pbmc1.csv')%>% column_to_rownames('X')
a66 <- (grep('^pbmc2_Seqwell',rownames(data)))
pbmc2_Seqwell <- data[a66,]
labels <- read.csv("Seq-Well/SW_pbmc1Labels.csv")
meta.dat <- labels[a66,, drop =F]
save(pbmc2_Seqwell, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc2_Seqwell.Rdata')

data <- read.csv('Smart-Seq2/SM2_pbmc1.csv')%>% column_to_rownames('X')
a77 <- (grep('^pbmc2_SM2',rownames(data)))
pbmc2_SM2 <- data[a77,]
labels <- read.csv("Smart-Seq2/SM2_pbmc1Labels.csv")
meta.dat <- labels[a77,, drop =F]
save(pbmc2_SM2, meta.dat, file = '../../cleaned_Lu/pbmcbench/pbmc2_SM2.Rdata')

#### --------------Prepare DEs for scMapping; Seurat MAST DE analysis ------ ####
## Cast to cluster run. refer to code: run.R and run.sh
setwd("/research/bsi/projects/staff_analysis/m216453/scmapping")
files <- list.files('data/GenomeBiology2019_Tamin/cleaned_Lu/', pattern = 'RData')
for(file in files){
  cat('[',file,']\n')
  try({
    seurat_obj <- seurat.markers <- NULL;
    tmp <- load(paste0('data/GenomeBiology2019_Tamin/cleaned_Lu/',file))
    cat(length(tmp),'\n')
    seurat_obj <- NormalizeData(object = seurat_obj)
    seurat_obj <- FindVariableFeatures(object = seurat_obj)
    seurat_obj <- ScaleData(object = seurat_obj)
    seurat_obj <- RunPCA(object = seurat_obj)
    seurat_obj <- FindNeighbors(object = seurat_obj)
    seurat_obj <- FindClusters(object = seurat_obj)
    seurat_obj <- RunTSNE(object = seurat_obj)
    seurat_obj@meta.data$celltype[is.na(seurat_obj@meta.data$celltype)] <- 'Unclassified'
    seurat_obj@meta.data$cluster <- paste0('Cluster ',as.numeric(as.factor(seurat_obj@meta.data$celltype)))
    Idents(object = seurat_obj) <- seurat_obj@meta.data$cluster
    seurat.markers <- FindAllMarkers(seurat_obj, method ='MAST')
    save(seurat_obj, seurat.markers, file = paste0('result/GenomeBiology2019_Tamin/',file))
    })
}

## PBMC part
setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/data/GenomeBiology2019_Tamin/cleaned_Lu/pbmcbench/')
library(dplyr);library(tibble);library(Seurat)
files <- list.files()
for(file in files){
  tm <- load(file)
  cat(file,': ')
  if(file =="pbmc1_10X_v2.Rdata" ){counts <- t(pbmc1_10X_v2)}
  if(file =="pbmc1_10X_v3.Rdata" ){counts <- t(pbmc1_10X_v3)}
  if(file =="pbmc1_Celseq2.Rdata" ){counts <- t(pbmc1_Celseq2)}
  if(file =="pbmc1_Drop.Rdata" ){counts <- t(pbmc1_Drop)}
  if(file =="pbmc1_SM2.Rdata" ){counts <- t(pbmc1_SM2)}
  if(file =="pbmc1_Seqwell.Rdata" ){counts <- t(pbmc1_Seqwell)}
  if(file =="pbmc1_inDrops.Rdata" ){counts <- t(pbmc1_inDrops)}
  if(file =="pbmc2_10X_v2.Rdata" ){counts <- t(pbmc2_10X_v2)}
  if(file =="pbmc2_Celseq2.Rdata" ){counts <- t(pbmc2_Celseq2)}
  if(file =="pbmc2_Drop.Rdata" ){counts <- t(pbmc2_Drop)}
  if(file =="pbmc2_SM2.Rdata" ){counts <- t(pbmc2_SM2)}
  if(file =="pbmc2_Seqwell.Rdata" ){counts <- t(pbmc2_Seqwell)}
  if(file =="pbmc2_inDrops.Rdata" ){counts <- t(pbmc2_inDrops)}
  
  meta.dat$cluster <- as.numeric(as.factor(meta.dat$x))
  meta.dat <- meta.dat %>% dplyr::rename(celltype =x)
  seurat_obj <- CreateSeuratObject(counts = counts)
  seurat_obj@meta.data <- meta.dat
  seurat_obj <- NormalizeData(object = seurat_obj)
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj)
  seurat_obj <- FindNeighbors(object = seurat_obj)
  seurat_obj <- FindClusters(object = seurat_obj)
  Idents(object = seurat_obj) <- seurat_obj@meta.data$cluster
  seurat.markers <- FindAllMarkers(seurat_obj, method = 'MAST')
  save(seurat.markers, seurat_obj, file = paste0('../../../../result/GenomeBiology2019_Tamin/PBMCbench/',file))
}


