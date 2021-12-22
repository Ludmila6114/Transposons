
setwd('~/2021_year/GENES_piRNA_automatic/trimmed/Unmapped_fasta/final_sam/')
#TPM normalization
fc = read.table('all_counts.out', header=TRUE, sep="\t", row.names=1 )
fc$Length_norm <- fc$Length/1000
fc$SRR060673.sam <- fc$SRR060673.sam/sum(fc$SRR060673.sam)
fc$SRR060673.sam <- fc$SRR060673.sam * 1000000 /fc$Length_norm
fc$SRR060687.sam <- fc$SRR060687.sam/sum(fc$SRR060687.sam)
fc$SRR060687.sam <- fc$SRR060687.sam * 1000000/fc$Length_norm
fc$SRR2096026.sam <- fc$SRR2096026.sam/sum(fc$SRR2096026.sam)
fc$SRR2096026.sam <- fc$SRR2096026.sam * 1000000/fc$Length_norm
fc$SRR2096038.sam <- fc$SRR2096038.sam /sum(fc$SRR2096038.sam )
fc$SRR2096038.sam  <- fc$SRR2096038.sam *1000000/fc$Length_norm
fc$SRR060688.sam <- fc$SRR060688.sam/sum(fc$SRR060688.sam)
fc$SRR060688.sam <- fc$SRR060688.sam * 1000000/fc$Length_norm
fc$SRR060663.sam <- fc$SRR060663.sam / sum(fc$SRR060663.sam)
fc$SRR060663.sam <- fc$SRR060663.sam * 1000000/fc$Length_norm
fc$SRR2096050.sam <- fc$SRR2096050.sam/sum(fc$SRR2096050.sam)
fc$SRR2096050.sam <- fc$SRR2096050.sam * 1000000/fc$Length_norm
fc$SRR2096061.sam <- fc$SRR2096061.sam/sum(fc$SRR2096061.sam)
fc$SRR2096061.sam <- fc$SRR2096061.sam * 1000000/fc$Length_norm
fc$SRR060661.sam <- fc$SRR060661.sam/sum(fc$SRR060661.sam)
fc$SRR060661.sam <- fc$SRR060661.sam * 1000000/fc$Length_norm
fc$SRR060654.sam <- fc$SRR060654.sam/sum(fc$SRR060654.sam)
fc$SRR060654.sam <- fc$SRR060654.sam * 1000000/fc$Length_norm
fc$SRR060662.sam <- fc$SRR060662.sam/sum(fc$SRR060662.sam)
fc$SRR060662.sam <- fc$SRR060662.sam*1000000/fc$Length_norm
fc$SRR060656.sam <- fc$SRR060656.sam/sum(fc$SRR060656.sam)
fc$SRR060656.sam <- fc$SRR060656.sam * 1000000/fc$Length_norm
fc$SRR2096072.sam <- fc$SRR2096072.sam/sum(fc$SRR2096072.sam)
fc$SRR2096072.sam <- fc$SRR2096072.sam * 1000000/fc$Length_norm
fc$SRR2096079.sam <- fc$SRR2096079.sam/sum(fc$SRR2096079.sam)
fc$SRR2096079.sam <- fc$SRR2096079.sam*1000000/fc$Length_norm
fc$SRR2096080.sam <- fc$SRR2096080.sam/sum(fc$SRR2096080.sam)
fc$SRR2096080.sam <- fc$SRR2096080.sam*1000000/fc$Length_norm
fc$SRR2096081.sam <- fc$SRR2096081.sam/sum(fc$SRR2096081.sam)
fc$SRR2096081.sam <- fc$SRR2096081.sam * 1000000/fc$Length_norm

data = fc

data$SRR060673.sam <- NULL
data$SRR060687.sam <- NULL
data$SRR2096026.sam <- NULL
data$SRR2096038.sam <- NULL
data$SRR060688.sam <- NULL
data$SRR060663.sam <- NULL
data$SRR2096050.sam <- NULL
data$SRR2096061.sam <- NULL
data$SRR060661.sam <- NULL
data$SRR060654.sam <- NULL
data$SRR060662.sam <- NULL
data$SRR060656.sam <- NULL
data$SRR2096072.sam <- NULL
data$SRR2096079.sam <- NULL
data$SRR2096080.sam <- NULL
data$SRR2096081.sam <- NULL

data$F0_NR_9_Ov <- fc$SRR060673.sam
data$F0_NR_9_Emb <- fc$SRR060687.sam
data$F0_JB_9_Emb <- (fc$SRR2096026.sam + fc$SRR2096038.sam) /2
data$F0_NR_160_Ov <- fc$SRR060688.sam
data$F0_NR_160_Emb <- fc$SRR060663.sam
data$F0_JB_160_Emb <- (fc$SRR2096050.sam + fc$SRR2096061.sam)/2
data$F1_NR_1609_Emb <- fc$SRR060661.sam
data$F1_NR_1609_Ov <- fc$SRR060654.sam
data$F1_NR_9160_Emb <- fc$SRR060662.sam
data$F1_NR_9160_Ov <- fc$SRR060656.sam
data$F2_JB_Dys_F2 <- (fc$SRR2096072.sam + fc$SRR2096079.sam)/2
data$F2_JB_Rec_F2 <- (fc$SRR2096080.sam + fc$SRR2096081.sam)/2

write.csv(data, 'UPD_TPM_counts_ALL_DATA_SMALLRNA.csv', quote = F, row.names = F)


data <- data[data$F0_NR_9_Ov >=100 & data$F0_NR_9_Emb >=100  & data$F0_JB_9_Emb >=100 & data$F0_NR_160_Ov >=100 &
               data$F0_NR_160_Emb >=100 & data$F0_JB_160_Emb >=100 &
               data$F1_NR_1609_Emb >=100 & data$F1_NR_1609_Ov >=100 &
               data$F1_NR_9160_Emb >=100 & data$F1_NR_9160_Ov >=100 &
               data$F2_JB_Dys_F2 >=100  & data$F2_JB_Rec_F2 >=100 ,]

data <- data[,c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)]


ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
data$dmel_ort <- ort[match(rownames(data), ort$V6),]$V1
data$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(data$dmel_ort), 'SYMBOL','FLYBASE'))



GO <- data
library(clusterProfiler)
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

dotplot(ego2)
data$gene_id <- rownames(data)

data$Gene <- NA
data$Gene <- ifelse(data$dmel.sym != 'NULL', data$dmel.sym, as.character(data$gene_id))
data$Gene <- ifelse(is.na(data$Gene), data$gene_id, data$Gene)
data$gene_id <- NULL
data$dmel_ort <- NULL
data$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(data)
library('gplots')
ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 4, space = "Lab",
                       name="Log10(TPM)") + ggtitle('381 genes, log10(TPM), piRNA reads') +
  theme_bw() + xlab('Drosophila virilis strain') + theme(axis.text.x = element_text(angle = 30 , hjust = 1, size = 10))


cdi <- data[data$F0_JB_160_Emb > data$F0_JB_9_Emb,]

df_melted <- melt(cdi)
library('gplots')
ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 2, space = "Lab",
                       name="Log10(TPM)") + ggtitle('50 genes, (160F0 embryo > 9F0 Embryo)\nlog10(TPM), piRNA reads') +
  theme_bw() + xlab('Drosophila virilis strain') + theme(axis.text.x = element_text(angle = 30 , hjust = 1, size = 10))












ggplot(data = data, aes(x = log10(F0_JB_9_Emb), y = log10(F0_JB_160_Emb))) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', col = 'red') +
  xlab('Log10(RPM, 9 strain, depth ~ 60k)') +
  ylab('LOg10(RPM, 160 strain, depth ~ 43k)') +
  ggtitle('piRNA to genes, RPM')
#gene filtration
genes_160 <- featurecounts[featurecounts$P_160 > 200,]
genes_9 <- featurecounts[featurecounts$Nine > 200,]
RPM = intersect(rownames(genes_160), rownames(genes_9))

#56 common genes
df = data.frame(gene_id = rownames(genes_160[rownames(genes_160) %in% rownames(genes_9),]), Strain9 = NA, Strain160 = NA)
for (i in 1:nrow(df)){
  df$Strain9[i] <- ifelse(df$gene_id[i] %in% rownames(genes_9), 
                          genes_9[which(rownames(genes_9) %in% df$gene_id[i]),]$Nine, 0)
  df$Strain160[i] <- ifelse(df$gene_id[i] %in% rownames(genes_160), 
                            genes_160[which(rownames(genes_160) %in% df$gene_id[i]),]$P_160, 0)
  print(i)
}
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
df$dmel_ort <- ort[match(df$gene_id, ort$V6),]$V1
df$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(df$dmel_ort), 'SYMBOL','FLYBASE'))

GO <- df
library(clusterProfiler)
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

dotplot(ego2)

#heatmap
df$Gene <- NA
df$Gene <- ifelse(df$dmel.sym != 'NULL', df$dmel.sym, as.character(df$gene_id))
df$gene_id <- NULL
df$dmel_ort <- NULL
df$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(df)
library('gplots')
ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 4, space = "Lab", 
                       name="Log10(RPM reads)") + ggtitle('56 genes, log10(RPM), piRNA reads') +
  theme_bw() + xlab('Drosophila virilis strain')


#another genes, 9, but not in intersection
df = data.frame(gene_id = rownames(genes_160[!rownames(genes_160) %in% rownames(genes_9),]), Strain9 = NA, Strain160 = NA)
for (i in 1:nrow(df)){
  df$Strain9[i] <- ifelse(df$gene_id[i] %in% rownames(genes_9), 
                          genes_9[which(rownames(genes_9) %in% df$gene_id[i]),]$Nine, 0)
  df$Strain160[i] <- ifelse(df$gene_id[i] %in% rownames(genes_160), 
                            genes_160[which(rownames(genes_160) %in% df$gene_id[i]),]$P_160, 0)
  print(i)
}
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
df$dmel_ort <- ort[match(df$gene_id, ort$V6),]$V1
df$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(df$dmel_ort), 'SYMBOL','FLYBASE'))
#df[226, 'dmel.sym'] = 'Arc1_again'
GO <- df
library(clusterProfiler)
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2)
#heatmap
df$Gene <- NA
df$Gene <- ifelse(df$dmel.sym != 'NULL', df$dmel.sym, as.character(df$gene_id))
df$gene_id <- NULL
df$dmel_ort <- NULL
df$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(df)
library('gplots')

ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "orange", 
                       midpoint = 2, space = "Lab",
                       name="log10(RPM reads)") + ggtitle('log10(RPM)\npiRNA only in 160 strain (>200RPM)') +
  theme_bw() + xlab('Drosophila virilis strain')












#TPM normalization
setwd('R_filtered_piRNA/')
featurecounts = read.table('piRNA_counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(5, 6, 7, 8, 9)]
featurecounts$Length <- featurecounts$Length / 1000
featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam/sum(featurecounts$pirna26.fasta.sam)
featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam/sum(featurecounts$pirna38.fasta.sam)

featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam / sum(featurecounts$pirna50.fasta.sam)
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam / sum(featurecounts$pirna61.fasta.sam)
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam * 1000000 / featurecounts$Length
featurecounts$Nine <- (featurecounts$pirna26.fasta.sam + featurecounts$pirna38.fasta.sam)/2
featurecounts$P_160 <- (featurecounts$pirna50.fasta.sam + featurecounts$pirna61.fasta.sam)/2

ggplot(data = featurecounts, aes(x = log10(Nine), y = log10(P_160))) +
  geom_point() +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', col = 'red') +
  xlab('Log10(TPM, 9 strain, depth ~ 60k)') +
  ylab('LOg10(TPM, 160 strain, depth ~ 43k)') +
  ggtitle('piRNA to genes, TPM')


#gene filtration
genes_160 <- featurecounts[featurecounts$P_160 > 200,]
genes_9 <- featurecounts[featurecounts$Nine > 200,]
TPM = intersect(rownames(genes_160), rownames(genes_9))

#70 common genes
df = data.frame(gene_id = rownames(genes_160[rownames(genes_160) %in% rownames(genes_9),]), Strain9 = NA, Strain160 = NA)
for (i in 1:nrow(df)){
  df$Strain9[i] <- ifelse(df$gene_id[i] %in% rownames(genes_9), 
                          genes_9[which(rownames(genes_9) %in% df$gene_id[i]),]$Nine, 0)
  df$Strain160[i] <- ifelse(df$gene_id[i] %in% rownames(genes_160), 
                            genes_160[which(rownames(genes_160) %in% df$gene_id[i]),]$P_160, 0)
  print(i)
}
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
df$dmel_ort <- ort[match(df$gene_id, ort$V6),]$V1
df$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(df$dmel_ort), 'SYMBOL','FLYBASE'))

GO <- df
library(clusterProfiler)
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

dotplot(ego2)

extract = as.data.frame(ego2)
genes = unlist(strsplit(extract$geneID[1], '/'))
genes <- c(genes)
df = data.frame(gene = genes)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
df$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(df$gene), 'SYMBOL','FLYBASE'))
df$dmel.sym






#heatmap
df$Gene <- NA
df$Gene <- ifelse(df$dmel.sym != 'NULL', df$dmel.sym, as.character(df$gene_id))
df$gene_id <- NULL
df$dmel_ort <- NULL
df$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(df)
library('gplots')
ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 4, space = "Lab", 
                       name="log10(TPM reads)") + ggtitle('70 genes, log10(TPM), piRNA reads') +
  theme_bw() + xlab('Drosophila virilis strain')





#another genes, 9, but not in intersection
df = data.frame(gene_id = rownames(genes_9[!rownames(genes_9) %in% rownames(genes_160),]), Strain9 = NA, Strain160 = NA)
for (i in 1:nrow(df)){
  df$Strain9[i] <- ifelse(df$gene_id[i] %in% rownames(genes_9), 
                          genes_9[which(rownames(genes_9) %in% df$gene_id[i]),]$Nine, 0)
  df$Strain160[i] <- ifelse(df$gene_id[i] %in% rownames(genes_160), 
                            genes_160[which(rownames(genes_160) %in% df$gene_id[i]),]$P_160, 0)
  print(i)
}
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
df$dmel_ort <- ort[match(df$gene_id, ort$V6),]$V1
df$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(df$dmel_ort), 'SYMBOL','FLYBASE'))


GO <- df
library(clusterProfiler)
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

dotplot(ego2)

#heatmap
df$Gene <- NA
df$Gene <- ifelse(df$dmel.sym != 'NULL', df$dmel.sym, as.character(df$gene_id))
df$gene_id <- NULL
df$dmel_ort <- NULL
df$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(df)
library('gplots')

ggplot(data = df_melted[128:188,], aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "orange", 
                       midpoint = 2, space = "Lab", limits = c(2.2, 4.2),
                       name="log10(TPM reads)") + ggtitle('128-188 genes log10(TPM)\npiRNA only in 9 strain (>200TPM)') +
  theme_bw() + xlab('Drosophila virilis strain')



#На какие гены смотрим дальше?
GENES = intersect(RPM, TPM)
save = genes_160[rownames(genes_160) %in% GENES,]
write.csv(save, 'Important_piRNA_genes_37.csv', quote = F, row.names = T)




#F1, NR data: F1 
setwd('/mnt/raid/pro_milka/2021_year/GENES_PIRNA/raw_data_NR_andF2_JB/trimmed/small_RNA_after_filter')

library(Biostrings)
lib_81 <- readDNAStringSet('SRR2096081_trimmed.fq.piRNA.fasta') 
lib_81 <- lib_81[lib_81@ranges@width > 22,]
writeXStringSet(lib_81, '/mnt/raid/pro_milka/2021_year/GENES_PIRNA/raw_data_NR_andF2_JB/trimmed/small_RNA_after_filter/R_filtered/piRNA81_JB_F2.fasta')





#ALL
setwd('~/2021_year/GENES_PIRNA/raw_data_NR_andF2_JB/trimmed/small_RNA_after_filter/R_filtered/')
featurecounts = read.table('all_counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts$Length <- featurecounts$Length / 1000


featurecounts$pirna26.fasta.sam <- featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna26.fasta.sam 
featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna26.fasta.sam <- NULL
featurecounts$pirna38.fasta.sam <- featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna38.fasta.sam
featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna38.fasta.sam <- NULL
featurecounts$pirna50.fasta.sam <- featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna50.fasta.sam
featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna50.fasta.sam <- NULL
featurecounts$pirna61.fasta.sam <- featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna61.fasta.sam
featurecounts$X.mnt.raid.pro_milka.2021_year.GENES_PIRNA.raw_data_F0.trimmed.small_rna_fasta_after_filtration_script.R_filtered_piRNA.pirna61.fasta.sam <- NULL


featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam/sum(featurecounts$pirna26.fasta.sam)
featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam/sum(featurecounts$pirna38.fasta.sam)
featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam / sum(featurecounts$pirna50.fasta.sam)
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam * 1000000 / featurecounts$Length
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam / sum(featurecounts$pirna61.fasta.sam)
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F0_JB_9_Emb <- (featurecounts$pirna26.fasta.sam + featurecounts$pirna38.fasta.sam)/2
featurecounts$F0_JB_160_Emb <- (featurecounts$pirna50.fasta.sam + featurecounts$pirna61.fasta.sam)/2
featurecounts$Chr <- NULL
featurecounts$Start <- NULL
featurecounts$End <- NULL
featurecounts$Strand <- NULL
featurecounts$pirna26.fasta.sam <- NULL
featurecounts$pirna38.fasta.sam <- NULL
featurecounts$pirna50.fasta.sam <- NULL
featurecounts$pirna61.fasta.sam <- NULL


featurecounts$piRNA87.fasta.sam <- featurecounts$piRNA87.fasta.sam/sum(featurecounts$piRNA87.fasta.sam)
featurecounts$piRNA87.fasta.sam <- featurecounts$piRNA87.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F0_NR_9_Emb <- featurecounts$piRNA87.fasta.sam
featurecounts$piRNA87.fasta.sam <- NULL

featurecounts$piRNA63.fasta.sam <- featurecounts$piRNA63.fasta.sam/sum(featurecounts$piRNA63.fasta.sam)
featurecounts$piRNA63.fasta.sam <- featurecounts$piRNA63.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F0_NR_160_Emb <- featurecounts$piRNA63.fasta.sam
featurecounts$piRNA63.fasta.sam <- NULL


featurecounts$piRNA61.fasta.sam <- featurecounts$piRNA61.fasta.sam/sum(featurecounts$piRNA61.fasta.sam)
featurecounts$piRNA61.fasta.sam <- featurecounts$piRNA61.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F1_NR_160_9_Emb <- featurecounts$piRNA61.fasta.sam
featurecounts$piRNA61.fasta.sam <- NULL

featurecounts$piRNA62.fasta.sam <- featurecounts$piRNA62.fasta.sam/sum(featurecounts$piRNA62.fasta.sam)
featurecounts$piRNA62.fasta.sam <- featurecounts$piRNA62.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F1_NR_9_160_Emb <- featurecounts$piRNA62.fasta.sam
featurecounts$piRNA62.fasta.sam <- NULL


featurecounts$piRNA54.fasta.sam <- featurecounts$piRNA54.fasta.sam/sum(featurecounts$piRNA54.fasta.sam)
featurecounts$piRNA54.fasta.sam <- featurecounts$piRNA54.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F1_NR_160_9_Ovaries <- featurecounts$piRNA54.fasta.sam
featurecounts$piRNA54.fasta.sam <- NULL

featurecounts$piRNA56.fasta.sam <- featurecounts$piRNA56.fasta.sam/sum(featurecounts$piRNA56.fasta.sam)
featurecounts$piRNA56.fasta.sam <- featurecounts$piRNA56.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F1_NR_9_160_Ovaries <- featurecounts$piRNA56.fasta.sam
featurecounts$piRNA56.fasta.sam <- NULL


featurecounts$piRNA72_JB_F2.fasta.sam <- featurecounts$piRNA72_JB_F2.fasta.sam/sum(featurecounts$piRNA72_JB_F2.fasta.sam)
featurecounts$piRNA72_JB_F2.fasta.sam <- featurecounts$piRNA72_JB_F2.fasta.sam * 1000000 / featurecounts$Length
featurecounts$piRNA79_JB_F2.fasta.sam <- featurecounts$piRNA79_JB_F2.fasta.sam/sum(featurecounts$piRNA79_JB_F2.fasta.sam)
featurecounts$piRNA79_JB_F2.fasta.sam <- featurecounts$piRNA79_JB_F2.fasta.sam * 1000000 / featurecounts$Length

featurecounts$piRNA80_JB_F2.fasta.sam <- featurecounts$piRNA80_JB_F2.fasta.sam / sum(featurecounts$piRNA80_JB_F2.fasta.sam)
featurecounts$piRNA80_JB_F2.fasta.sam <- featurecounts$piRNA80_JB_F2.fasta.sam * 1000000 / featurecounts$Length
featurecounts$piRNA81_JB_F2.fasta.sam <- featurecounts$piRNA81_JB_F2.fasta.sam/ sum(featurecounts$piRNA81_JB_F2.fasta.sam)
featurecounts$piRNA81_JB_F2.fasta.sam <- featurecounts$piRNA81_JB_F2.fasta.sam * 1000000 / featurecounts$Length
featurecounts$F2_JB_rec_Emb <- (featurecounts$piRNA80_JB_F2.fasta.sam + featurecounts$piRNA81_JB_F2.fasta.sam)/2
featurecounts$F2_JB_dys_Emb <- (featurecounts$piRNA72_JB_F2.fasta.sam + featurecounts$piRNA79_JB_F2.fasta.sam)/2

featurecounts$piRNA72_JB_F2.fasta.sam <- NULL
featurecounts$piRNA79_JB_F2.fasta.sam <- NULL
featurecounts$piRNA80_JB_F2.fasta.sam <- NULL
featurecounts$piRNA81_JB_F2.fasta.sam <- NULL
featurecounts$Length <- NULL

#filter genes
setwd('~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/')
genes <- read.csv('R_filtered_piRNA/Important_piRNA_genes_37.csv')
genes <- genes$X

genes <- c('FBgn0204789', 'FBgn0203577', 'FBgn0205401', 'FBgn0282812', 'FBgn0208364')
genes <- c('FBgn0208016', 'FBgn0210197', 'FBgn0283301', 'FBgn0256614', 'FBgn0199144', 'FBgn0283080')



Intersect <- featurecounts[rownames(featurecounts) %in% genes,]
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
Intersect$dmel_ort <- ort[match(rownames(Intersect), ort$V6),]$V1
Intersect$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(Intersect$dmel_ort), 'SYMBOL','FLYBASE'))
Intersect$Gene <- NA
Intersect$Gene <- ifelse(Intersect$dmel.sym != 'NULL', Intersect$dmel.sym, rownames(Intersect))
Intersect$dmel_ort <- NULL
Intersect$dmel.sym <- NULL
library(reshape2)
df_melted <- melt(Intersect)
library('gplots')

ggplot(data = df_melted, aes(x=variable, y=Gene, fill=log10(value))) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 2.5, space = "Lab",
                       name="log10\n(TPM reads)") + ggtitle('log10(TPM) for genes with piRNA only in 160 strain') +
  theme_bw() + xlab('Progeny')

















#RNAseq F0, JB

setwd('~/2021_year/Kelleher/F0_RNA_seq_ref/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts_parents_RNA.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11)]
cond_1 = rep("9", 3)
cond_2 = rep("160", 3)
samples = names(featurecounts)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)
colData
dds = DESeqDataSetFromMatrix(countData=featurecounts, colData=colData, design = ~condition)
dds$condition = relevel(dds$condition, ref = "9")
dds = DESeq(dds)

#scale factors:
dds$sizeFactor
head(counts(dds))
head(counts(dds, normalized = TRUE))
head(featurecounts)
normalized = as.data.frame(counts(dds, normalized = TRUE))
normalized$Nine <- (normalized$SRR2096015_trimmed.fq.bam + normalized$SRR2096016_trimmed.fq.bam + normalized$SRR2096017_trimmed.fq.bam)/3
normalized$P <- (normalized$SRR2096018_trimmed.fq.bam + normalized$SRR2096019_trimmed.fq.bam + normalized$SRR2096020_trimmed.fq.bam)/3

genes_with_piRNA_in_160 <- data.frame(gene = rownames(genes_160), TPM = genes_160$P_160, basemean = NA)
for (i in 1:nrow(genes_with_piRNA_in_160)){
  genes_with_piRNA_in_160$basemean[i] <- normalized[which(rownames(normalized) %in% genes_with_piRNA_in_160$gene[i]),]$P
  print(i)
}

ggplot(data = genes_with_piRNA_in_160, aes(x = log10(TPM), y=log10(basemean))) + 
  geom_point() + theme_bw() + geom_smooth(method='lm') + 
  ggtitle('160 strain, 83 genes with high piRNA level') +
  xlab('log10(TPM), piRNA') + ylab('log10(basemean), mRNA')



cor.test(genes_with_piRNA_in_160$TPM, genes_with_piRNA_in_160$basemean)



genes_with_piRNA_in_9 <- data.frame(gene = rownames(genes_9), TPM = genes_9$Nine, basemean = NA)
for (i in 1:nrow(genes_with_piRNA_in_9)){
  genes_with_piRNA_in_9$basemean[i] <- normalized[which(rownames(normalized) %in% genes_with_piRNA_in_9$gene[i]),]$Nine
  print(i)
}

ggplot(data = genes_with_piRNA_in_9, aes(x = log10(TPM), y=log10(basemean))) + 
  geom_point() + theme_bw() + geom_smooth(method='lm') + 
  ggtitle('9 strain, 83 genes with high piRNA level') +
  xlab('log10(TPM), piRNA') + ylab('log10(basemean), mRNA')



cor.test(genes_with_piRNA_in_160$TPM, genes_with_piRNA_in_160$basemean)




before = readDNAStringSet('~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/SRR2096026_trimmed.fq.piRNA.fasta')
table(before@ranges@width)
after = readDNAStringSet('~/filtered2.fasta')
table(after@ranges@width)
after
