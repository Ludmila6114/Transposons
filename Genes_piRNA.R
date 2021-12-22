#F0, library length filtration

setwd('~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/')
library(Biostrings)
lib_26 <- readDNAStringSet('SRR2096026_trimmed.fq.piRNA.fasta') 
lib_26 <- lib_26[lib_26@ranges@width > 22,]
writeXStringSet(lib_26, '~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/R_filtered_piRNA/pirna26.fasta')
lib_38 <- readDNAStringSet('SRR2096038_trimmed.fq.piRNA.fasta') 
lib_38 <- lib_38[lib_38@ranges@width > 22,]
writeXStringSet(lib_38, '~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/R_filtered_piRNA/pirna38.fasta')
lib_50 <- readDNAStringSet('SRR2096050_trimmed.fq.piRNA.fasta') 
lib_50 <- lib_50[lib_50@ranges@width > 22,]
writeXStringSet(lib_50, '~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/R_filtered_piRNA/pirna50.fasta')
lib_61 <- readDNAStringSet('SRR2096061_trimmed.fq.piRNA.fasta') 
lib_61 <- lib_61[lib_61@ranges@width > 22,]
writeXStringSet(lib_61, '~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/R_filtered_piRNA/pirna61.fasta')


setwd('~/2021_year/GENES_PIRNA/raw_data_F0/trimmed/small_rna_fasta_after_filtration_script/R_filtered_piRNA/')
#RPM normalization
featurecounts = read.table('piRNA_counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9)]
featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam/sum(featurecounts$pirna26.fasta.sam)
featurecounts$pirna26.fasta.sam <- featurecounts$pirna26.fasta.sam * 1000000 
featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam/sum(featurecounts$pirna38.fasta.sam)
featurecounts$pirna38.fasta.sam <- featurecounts$pirna38.fasta.sam * 1000000
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam / sum(featurecounts$pirna50.fasta.sam)
featurecounts$pirna50.fasta.sam <- featurecounts$pirna50.fasta.sam * 1000000
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam / sum(featurecounts$pirna61.fasta.sam)
featurecounts$pirna61.fasta.sam <- featurecounts$pirna61.fasta.sam * 1000000
featurecounts$Nine <- (featurecounts$pirna26.fasta.sam + featurecounts$pirna38.fasta.sam)/2
featurecounts$P_160 <- (featurecounts$pirna50.fasta.sam + featurecounts$pirna61.fasta.sam)/2
ggplot(data = featurecounts, aes(x = log10(Nine), y = log10(P_160))) +
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




