library(scMapping)
scMappingDatabase[1:4,1:4]
dim(scMappingDatabase)
unique(scMayoMapDatabase$tissue)
for(tissue in unique(scMayoMapDatabase$tissue)){
  tmp <- scMayoMapDatabase[scMayoMapDatabase$tissue ==tissue,-1] 
  rownames(tmp) <- NULL
  colnames(tmp) <- gsub('.*:','',colnames(tmp))
  tmp0 <- tmp %>% pivot_longer(cols = -gene, names_to = "CellType", values_to = "Presence") %>%  
    filter(Presence == 1) %>% select(-Presence)
  tmp0 <- tmp0[,c(2,1)]
  write.table(tmp0, paste0('/research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/',tissue,'.tsv'), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

xx <- scMayoMapDatabase[,-1] %>% pivot_longer(cols = -gene, names_to = "CellType", values_to = "Presence") %>%  
  filter(Presence == 1) %>% select(-Presence)
xx <- xx[,c(2,1)]
write.table(xx, paste0('/research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/scMayoMapDatabase.tsv'), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
