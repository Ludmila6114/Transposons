#report
result = data.frame(genome = NA, TE = NA, Sum_length = NA)

#border
chinaN = read.table('~/compare_assemblies_plotly/November2021/report/masker/chinaN/chinaN_start.out', fill = T)
chinaN = chinaN[3:26,]
chinaN = na.omit(chinaN)
chinaN$TE = as.numeric(as.character(chinaN$V7)) - as.numeric(as.character(chinaN$V6))
#Просуммировать в каждой группе
one = as.data.frame(aggregate(chinaN$TE, by=list(chinaN$V10), FUN=sum))
chinaN_df = data.frame(genome = 'chinaN', TE = one$Group.1, Sum_length = one$x)
result = rbind(result, chinaN_df)
result = na.omit(result)


chinaM = read.table('~/compare_assemblies_plotly/November2021/report/masker/chinaM/chinaM_start.out', fill = T)
chinaM = chinaM[3:18,]
chinaM = na.omit(chinaM)
chinaM$TE = as.numeric(as.character(chinaM$V7)) - as.numeric(as.character(chinaM$V6))
two = as.data.frame(aggregate(chinaM$TE, by=list(chinaM$V10), FUN=sum))
chinaM_df = data.frame(genome = 'chinaM', TE = two$Group.1, Sum_length = two$x)
result = rbind(result, chinaM_df)

N_101 = read.table('~/compare_assemblies_plotly/November2021/report/masker/101N/101N_start.out', fill = T)
N_101 = N_101[3:11,]
N_101 = na.omit(N_101)
N_101$TE = as.numeric(as.character(N_101$V7)) - as.numeric(as.character(N_101$V6))
three = as.data.frame(aggregate(N_101$TE, by=list(N_101$V10), FUN=sum))
N_df = data.frame(genome = 'JapanN', TE = three$Group.1, Sum_length = three$x)
result = rbind(result, N_df)

M_101 = read.table('~/compare_assemblies_plotly/November2021/report/masker/101M/101M_start.out', fill = T)
M_101 = M_101[3:11,]
M_101 = na.omit(M_101)
M_101$TE = as.numeric(as.character(M_101$V7)) - as.numeric(as.character(M_101$V6))
fo = as.data.frame(aggregate(M_101$TE, by=list(M_101$V10), FUN=sum))
fo_df = data.frame(genome = 'JapanM', TE = fo$Group.1, Sum_length = fo$x)
result = rbind(result, fo_df)


P = read.table('~/compare_assemblies_plotly/November2021/report/masker/160/160_start.out', fill = T)
P = P[3:17,]
P = na.omit(P)
P$TE = as.numeric(as.character(P$V7)) - as.numeric(as.character(P$V6))
P_df = as.data.frame(aggregate(P$TE, by=list(P$V10), FUN=sum))
five = data.frame(genome = '160', TE = P_df$Group.1, Sum_length = P_df$x)
result = rbind(result, five)

ggplot(data = result, aes(y = Sum_length, fill = TE, x = genome)) + 
  geom_bar(stat = 'identity') + theme_bw() + xlab('Genome') + ylab('Sum TE length') + 
  ggtitle('Sum TE length within +- 100Kb windows\nfrom an inversia border')



####################
cont = read.table('~/compare_assemblies_plotly/November2021/query160/target_101M_query_160.paf', fill = T)
cont = cont[cont$V11 > 50000,]
cont$rez2 = paste(cont$V6,'_',cont$V5,'_',cont$V1)
cont$rez1 = paste(cont$V1,'_',cont$V5,'_',cont$V6)
#9_101N in table
View(table(cont$rez1))


open = read.table('~/compare_assemblies_plotly/November2021/report/101M_endX', fill = T)
open = open[open$V11 > 50000,]

########################################################
#12730300 to 17297960, CM017607.2
#Annotation

annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']
Gene <- Gene[seqnames(Gene) %in% 'CM017605.2',]
region = GRanges(seqnames = 'CM017607.2', ranges = IRanges(start = 12730300, end = 17297960))
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
Gene$dmel_ort <- ort[match(Gene$gene_id, ort$V6),]$V1
Gene$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(Gene$dmel_ort), 'SYMBOL','FLYBASE'))

Gene = Gene[start(Gene) >= 12730300 & end(Gene) <= 17297960,]
View(as.data.frame(Gene))

library(clusterProfiler)
Gene$dmel.sym
Gene.ENSEMBL <- bitr(Gene$dmel.sym, fromType = "SYMBOL",
                     toType = c("ENSEMBL"),
                     OrgDb = org.Dm.eg.db)

Gene.ENSEMBL
ego2 <- enrichGO(gene         = Gene.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ego2, showCategory=10)
Gene = as.data.frame(Gene)
Gene$X. <- NULL
write.csv(Gene, '~/Inversia_analysis/Genes_433_CM017607.2_inversia.csv')


##################################synteny in annotation:
#BEST_genome assembly fasta

china = rtracklayer::import('/mnt/raid/alexr/genome.annotation/china.neutral/china_neutral.gtf')
china = china[seqnames(china) == 'contig_36']
china = china[china$type == 'gene']
View(as.data.frame(china))
china = china[order(start(china)),]


annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']
Gene = Gene[Gene$gene_id %in% china$gene_id]
Gene <- Gene[seqnames(Gene) %in% 'CM017607.2',]
Gene <- Gene[order(start(Gene)),]
View(as.data.frame(Gene))


china = china[china$gene_id %in% Gene$gene_id]
which(Gene$gene_id %in% 'FBgn0205280')
df = data.frame(china = china$gene_id, p = rev(Gene$gene_id), number = NA)
#write.csv(df, '~/Inversia_analysis/DATA_FRAME_ANNOTATIONS')

for (i in 1:nrow(df)){
  if (as.character(df$china[i]) %in% df$p){
    df$number[i] = which(df$p %in% df$china[i])
  }
}

df$rez <- df$number == rownames(df)

#primers, what works? what didnt
###################




open = read.table('~/compare_assemblies_plotly/November2021/report/101N_end2', fill = T)
open = open[open$V11 > 50000,]
open = na.omit(open)










################################################################
#translocations
#annotation <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
#annotation <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/china.neutral/china_neutral.gtf')
annotation <- rtracklayer::import('/mnt/raid/sergey/bio-first/genomes_virilis/101N_annotation.gtf')

annotation <- annotation[annotation$type == 'gene']
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
annotation$dmel_ort <- ort[match(annotation$gene_id, ort$V6),]$V1
annotation$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(annotation$dmel_ort), 'SYMBOL','FLYBASE'))
table(is.na(annotation$dmel_ort))
#11795 with orts, 3242 without
annotation = annotation[!is.na(annotation$dmel_ort),]
dmel = rtracklayer::import('~/DMEL_ANNOTATION/dmel-all-r6.42.gtf')
dmel = dmel[dmel$type == 'gene']
table(seqnames(dmel))
dmel = dmel[seqnames(dmel) %in% c('2L', '2R', '3L', '3R', '4', 'X', 'Y', 'mitochondrion_genome')]
#dmel = dmel[seqnames(dmel) %in% c('X')]

#25 dmel was not in annotation
annotation = annotation[annotation$dmel_ort %in% dmel$gene_id,]
annotation$dmel_chr = seqnames(dmel[match(annotation$dmel_ort, dmel$gene_id),])

#View(as.data.frame(annotation))
#annotation = annotation[seqnames(annotation) %in% c('CM017604.2', 'CM017605.2', 'CM017606.2',
#                                                    'CM017607.2', 'CM017608.2', 'CM017609.2')]

#annotation = annotation[seqnames(annotation) %in% c('contig_50', 'scaffold_82', 'contig_84', 'contig_51',
#                                                       'contig_36', 'contig_49', 'contig_118')]

annotation = annotation[seqnames(annotation) %in% c('contig_13', 'scaffold_33', 'contig_7', 'contig_35',
                                                    'contig_32', 'contig_6', 'scaffold_22', 'contig_34')]


annotation$intrigue = paste(seqnames(annotation),'_',annotation$dmel_chr)
ggplot(data = as.data.frame(annotation), aes(x = intrigue)) + 
  geom_bar(stat = 'count', fill = 'indianred2') +
  coord_flip() +
  theme_bw() +
  ggtitle('Chromosomes in 160 Drosophila virilis\ncompare to Drosophila melanogaster\nfor 11756 genes') + 
  ylab('Chromosomes') +
  xlab('Count') +
  theme(text = element_text(size = 13)) +
  scale_y_continuous(breaks = c(0, 50, 100, 500, 1000, 1500, 2000), guide = guide_axis(angle = 45))


View(as.data.frame(table(as.character(annotation$intrigue))))

table(annotation$intrigue)
write.csv(annotation, '~/compare_assemblies_plotly/November2021/TRANSLOCATIONS/Genome_101N.csv', quote = F, row.names = F)


####################


rez_P = read.csv('~/compare_assemblies_plotly/November2021/TRANSLOCATIONS/Genome_160.csv')
rez_chinaN = read.csv('~/compare_assemblies_plotly/November2021/TRANSLOCATIONS/Genome_chinaN.csv')
rez_101N = read.csv('~/compare_assemblies_plotly/November2021/TRANSLOCATIONS/Genome_101N.csv')
table(rez_P$intrigue)

a = as.character(rez_P[rez_P$intrigue == 'CM017604.2 _ 2L',]$dmel.sym)
b = as.character(rez_chinaN[rez_chinaN$intrigue %in% c('contig_50 _ 2L'),]$dmel.sym)
c = as.character(rez_101N[rez_101N$intrigue %in% c('contig_13 _ 2L'),]$dmel.sym)

a = unique(a)
c = unique(c)
b = unique(b)

d = intersect(a, b)
m = intersect(d, c)
m = unique(m)

dmel = rtracklayer::import('~/DMEL_ANNOTATION/dmel-all-r6.42.gtf')
dmel = dmel[dmel$type == 'gene']

annotation <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
annotation <- annotation[annotation$type == 'gene']
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
annotation$dmel_ort <- ort[match(annotation$gene_id, ort$V6),]$V1
annotation$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(annotation$dmel_ort), 'SYMBOL','FLYBASE'))
annotation = annotation[!is.na(annotation$dmel_ort),]
annotation = annotation[annotation$dmel_ort %in% dmel$gene_id,]

annotation$is_TR = annotation$dmel.sym %in% m

annotation = as.data.frame(annotation)
annotation = annotation[annotation$is_TR == 'TRUE',]
table(annotation$seqnames)

an_delete = annotation[annotation$seqnames == 'CM017607.2',]
annotation = annotation[annotation$seqnames == 'CM017604.2' & ! annotation$dmel_ort %in% an_delete$dmel_ort,]

ggplot(data = annotation, 
       aes(x = start,
           y = end)) + geom_point()
annotation$dmel_ort
annotation$dmel.sym


View(as.data.frame(annotation))
View(as.data.frame(dmel))


###########dvir translocations
annotation1 <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
annotation2 <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/china.neutral/china_neutral.gtf')
annotation3 <- rtracklayer::import('/mnt/raid/sergey/bio-first/genomes_virilis/101N_annotation.gtf')
annotation1 <- annotation1[annotation1$type == 'gene']
annotation2 <- annotation2[annotation2$type == 'gene']
annotation3 <- annotation3[annotation3$type == 'gene']

library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
annotation1$dmel_ort <- ort[match(annotation1$gene_id, ort$V6),]$V1
annotation1$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(annotation1$dmel_ort), 'SYMBOL','FLYBASE'))

annotation2$dmel_ort <- ort[match(annotation2$gene_id, ort$V6),]$V1
annotation2$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(annotation2$dmel_ort), 'SYMBOL','FLYBASE'))

annotation3$dmel_ort <- ort[match(annotation3$gene_id, ort$V6),]$V1
annotation3$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(annotation3$dmel_ort), 'SYMBOL','FLYBASE'))



annotation1 = annotation1[seqnames(annotation1) %in% c('CM017604.2', 'CM017605.2', 'CM017606.2',
                                                    'CM017607.2', 'CM017608.2', 'CM017609.2')]

annotation2 = annotation2[seqnames(annotation2) %in% c('contig_50', 'scaffold_82', 'contig_84', 'contig_51',
                                                       'contig_36', 'contig_49', 'contig_118')]

annotation3 = annotation3[seqnames(annotation3) %in% c('contig_13', 'scaffold_33', 'contig_7', 'contig_35',
                                                    'contig_32', 'contig_6', 'scaffold_22', 'contig_34')]



annotation3[annotation3$gene_id %in% annotation1[seqnames(annotation1) == 'CM017609.2',]$gene_id &
              seqnames(annotation3) == 'contig_13',]$dmel.sym

seqnames(annotation3[annotation3$gene_id %in% N_P,])




View(as.data.frame(annotation1))
View(as.data.frame(annotation2))
View(as.data.frame(annotation3))





file = read.table('~/compare_assemblies_plotly/November2021/query160/target_chinaN_query_160.paf', fill = T)



##DE
#~/Inversia_analysis/Genes_433_CM017607.2_inversia.csv

#delele low-express
data = read.table('~/RNA_SEQ_LAB_GENOMES/for_pro_milka/alignements/Counts', header = TRUE)
data$ChinaN_rep1.bam = data$ChinaN_rep1.bam/data$Length
data$ChinaN_rep2.bam = data$ChinaN_rep2.bam/data$Length
data$ChinaM_rep1.bam = data$ChinaM_rep1.bam/data$Length
data$ChinaM_rep2.bam = data$ChinaM_rep2.bam/data$Length
data$X101N_rep1.bam = data$X101N_rep1.bam/data$Length
data$X101N_rep2.bam = data$X101N_rep2.bam/data$Length
data$X101M_rep1.bam = data$X101M_rep1.bam/data$Length
data$X101M_rep2.bam = data$X101M_rep2.bam/data$Length

data$ChinaN_rep1.bam = 1000000 * data$ChinaN_rep1.bam/(sum(data$ChinaN_rep1.bam))
data$ChinaN_rep2.bam = 1000000 * data$ChinaN_rep2.bam/(sum(data$ChinaN_rep2.bam))
data$ChinaM_rep1.bam = 1000000 * data$ChinaM_rep1.bam/(sum(data$ChinaM_rep1.bam))
data$ChinaM_rep2.bam = 1000000 * data$ChinaM_rep2.bam/(sum(data$ChinaM_rep2.bam))
data$X101N_rep1.bam = 1000000 * data$X101N_rep1.bam/(sum(data$X101N_rep1.bam))
data$X101N_rep2.bam = 1000000 * data$X101N_rep2.bam/(sum(data$X101N_rep2.bam))
data$X101M_rep1.bam = 1000000 * data$X101M_rep1.bam/(sum(data$X101M_rep1.bam))
data$X101M_rep2.bam = 1000000 * data$X101M_rep2.bam/(sum(data$X101M_rep2.bam))


df = data.frame(gene_id = data$Geneid, 
                ChinaN = 0.5*(data$ChinaN_rep1.bam + data$ChinaN_rep2.bam),
                ChinaM = 0.5*(data$ChinaM_rep1.bam + data$ChinaM_rep2.bam),
                N_101 = 0.5*(data$X101N_rep1.bam + data$X101N_rep2.bam),
                M_101 = 0.5*(data$X101M_rep1.bam + data$X101M_rep2.bam)
)
                


#only genes
genes = read.csv('~/Inversia_analysis/Genes_433_CM017607.2_inversia.csv')

df = df[df$gene_id %in% genes$gene_id,]

rownames(df) = df$gene_id
df$gene_id = NULL
df_melted = melt(df, variable.name = 'genome', value.name = gene_id)
df_melted = df_melted[df_melted$value > 10,]
df_melted = df_melted[df_melted$value < 1000,]


ggplot(data = df_melted, aes(x = df_melted$variable, y = df_melted$value)) + geom_boxplot()






library('DESeq2')
featurecounts = read.table('~/RNA_SEQ_LAB_GENOMES/for_pro_milka/alignements/Counts', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9)]
colData <- data.frame("replicate" = c("R1","R2", "R1","R2"), "treatment" = c("ChinaN","ChinaN", "ChinaM","ChinaM"))
attributes(colData)$row.names <- c("R1.N","R2.N","R1.M", "R2.M")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "ChinaM")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
resF2 <- results(dds50,contrast=c("treatment","ChinaN","ChinaM"))
resF2 <- as.data.frame(resF2)
expressed_and_passed_F2 <- resF2[resF2$padj < 0.05,]  



DE = expressed_and_passed_F2


library('DESeq2')
featurecounts = read.table('~/RNA_SEQ_LAB_GENOMES/for_pro_milka/alignements/Counts', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(10, 11, 12, 13)]
colData <- data.frame("replicate" = c("R1","R2", "R1","R2"), "treatment" = c("101N","101N", "101M","101M"))
attributes(colData)$row.names <- c("R1.N","R2.N","R1.M", "R2.M")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "101M")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
resF2 <- results(dds50,contrast=c("treatment","101N","101M"))
resF2 <- as.data.frame(resF2)
expressed_and_passed_F2 <- resF2[resF2$padj < 0.05,]  



length(intersect(rownames(DE), rownames(expressed_and_passed_F2)))


china_DE = DE[rownames(DE) %in% rownames(expressed_and_passed_F2),]
DE_101 = expressed_and_passed_F2[rownames(expressed_and_passed_F2) %in% rownames(DE),]
china_DE = china_DE[order(rownames(china_DE)),]
DE_101 = DE_101[order(rownames(DE_101)),]

rez = data.frame(gene_id = rownames(DE_101),
                 chinaN_DE = china_DE$log2FoldChange,
                 JapanN_DE = DE_101$log2FoldChange
)

ggplot(data = rez, aes(x = rez$chinaN_DE, y = rez$JapanN_DE)) + geom_point() + theme_bw() +
  geom_vline(xintercept = 0, col = 'red') +
  geom_hline(yintercept = 0, col = 'red') +
  xlab('log2 FC, chinaN vs chinaM') +
  ylab('log2 FC, 101N vs 101M')


library('DESeq2')
featurecounts = read.table('~/RNA_SEQ_LAB_GENOMES/for_pro_milka/alignements/Counts', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6:13)]
colData <- data.frame("replicate" = c("R1","R2", "R1","R2", "R1","R2", "R1","R2"), 
                      "treatment" = c("ChinaN","ChinaN", "ChinaM","ChinaM", "101N","101N", "101M","101M"))
attributes(colData)$row.names <- c("R1.cN","R2.cN","R1.cM", "R2.cM", "R1.JN","R2.JN","R1.JM", "R2.JM")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "101M")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
#resF2 <- results(dds50,contrast=c("treatment","101N","101M"))
#resF2 <- as.data.frame(resF2)
#expressed_and_passed_F2 <- resF2[resF2$padj < 0.05,]  

normalized = as.data.frame(counts(dds50, normalized = TRUE))
normalized = as.data.frame(t(normalized))
normalized$sample = colnames(normalized)

pca_res <- prcomp(normalized, scale. = TRUE)
autoplot(pca_res, colour = 'sample', data = normalized, size = 5)





