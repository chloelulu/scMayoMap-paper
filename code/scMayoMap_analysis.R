rm(list = ls())
library(dplyr);library(tibble);library(Seurat) 
library(scMapping)
setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/')

files <- list.files(pattern = 'Rdata')
file = files[12]
scMapping_res_sum <- plot.list <- list()
for(file in files){
  cat('[ ',file,' ]\n')
  seurat_obj <- NULL
  load(file)
  cat('-- scMapping --\n')
  t1 <- Sys.time()
  scmapping.res <- scmapping(data = seurat.markers, database = scMappingDatabase, tissue = 'blood')
  plt <- scmapping.plot(scmapping.object = scmapping.res)
  plot.list[[file]] <- plt
  t2 <- Sys.time()
  scmapping.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  annot <- as.data.frame(cbind(apply(scmapping.res$annotation.norm, 1, function(x) names(which.max(x))))) %>% rownames_to_column('cluster')
  colnames(annot)[2] <- 'scmapping'
  meta.dat <- seurat_obj@meta.data %>% dplyr::select(celltype, cluster)
  orig <- (unique(meta.dat) %>% mutate(cluster = paste0('Cluster ',cluster)))
  cell.num <- cbind.data.frame(table(seurat.markers$cluster))
  colnames(cell.num) <- c('cluster','cell.num')
  cell.num$cluster <- paste0('Cluster ',cell.num$cluster)
  scMapping_res <- full_join(orig, scmapping.res$res %>% mutate(cluster = paste0('Cluster ',cluster))%>% dplyr::rename(scMapping='celltype'), by = 'cluster') %>% full_join(cell.num)
  scMapping_res_sum[[file]] <- scMapping_res
  save(scMapping_res, scmapping.res, scmapping.time, file = paste0('annotation_res/',file, '_scmapping_new.Rdata'))

}





setwd('/research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/')
files <- list.files(pattern = 'RData$')
for(file in files){
  cat('[ ',file,' ]\n')
  load(file)
  seurat.markers <- seurat.markers %>% dplyr::rename(`avg_logFC`="avg_log2FC")
  if(file == "Muraro.RData"){seurat.markers$gene <- gsub('--.*','',seurat.markers$gene)}
  write.csv(seurat.markers, file = paste0(gsub('.Rdata','',file),'_seurat.markers.csv'))
  dt <- as.matrix(seurat_obj@assays$RNA@counts)
  if(file=='Muraro.RData'){rownames(dt) <- gsub('\\-\\-.*','',rownames(dt))}
  if(file=='CellBench_CelSeq2.RData'){
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                      values = rownames(dt), mart= mart)
    gene_IDs<- gene_IDs[which(gene_IDs$hgnc_symbol != ""),]
    dt <- dt[rownames(dt)%in% (gene_IDs$ensembl_gene_id),]
    gene_IDs <- gene_IDs[match(rownames(dt), gene_IDs$ensembl_gene_id),]
    rownames(dt) <- gene_IDs$hgnc_symbol
  }
  rownames(dt) <- toupper(rownames(dt))

  if(file %in% c("Xin.RData","Baron Human.RData","Baron Mouse.RData","Segerstolpe.RData","Muraro.RData")){tissue="pancreas"}
  if(file %in% c("CellBench_10x.RData","CellBench_CelSeq2.RData")){tissue="lung"}
  if(file %in% c("AMB.RData","MouseV1.RData","MouseALM.RData","HumanMTG.RData")){tissue="brain"}
  if(file %in% c("Zheng68K.RData","Zheng_sorted.RData")){tissue="brain"}
  if(file %in% c("TM.RData")){tissue = "muscle"}

  cat('-- scMapping --\n')
  try({
  if(file=='Muraro.RData'){seurat.markers$gene <- gsub('\\-\\-.*','',seurat.markers$gene)}
  if(file=='CellBench_CelSeq2.RData'){
    seurat.markers <- seurat.markers[seurat.markers$gene %in% gene_IDs$ensembl_gene_id,]
    gene_IDs1 <- gene_IDs %>% dplyr::rename(gene ='ensembl_gene_id')
    seurat.markers <- gene_IDs1 %>% right_join(seurat.markers) %>% dplyr::select(-gene) %>% dplyr::rename(gene='hgnc_symbol')
  }
  if(file %in% c("Zheng_sorted.RData","Zheng68K.RData")){
    db = database2
  }else{
    db = database
  }
  seurat.markers1 <- seurat.markers
  seurat.markers1$cluster <- as.character(seurat.markers1$cluster)
  seurat.markers1$cluster <- gsub('Cluster ','',seurat.markers1$cluster)
  t1 <- Sys.time()
  scmapping.res <- scMayoMap(data = seurat.markers1, database=db, tissue = tissue,
                              directory = '/research/bsi/projects/staff_analysis/m216453/scmapping/result/',plot =F)
  t2 <- Sys.time()
  scmapping.time <- as.numeric(difftime(t1,t2, units = 'secs'))
  annot <- scmapping.res$top.n[,c(1,5)]
  colnames(annot)[1] <- 'scmapping'
  meta.dat <- seurat_obj@meta.data %>% dplyr::select(celltype, cluster)
  orig <- unique(meta.dat)
  cell.num <- cbind.data.frame(table(seurat.markers1$cluster))
  colnames(cell.num) <- c('cluster','cell.num')
  cell.num$cluster <- paste0('Cluster ',cell.num$cluster)
  scMapping_res <- full_join(orig, annot, by = 'cluster') %>% full_join(cell.num)
  save(scMapping_res, scmapping.time, file = paste0('annotation_res/',file, '_scmapping.Rdata'))
  })
}

