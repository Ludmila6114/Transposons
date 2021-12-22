setwd('~/2021_year/Yeast_RNAseq/')
library('DESeq2')

length = read.table('SRR941816.out', header=TRUE, sep="\t", row.names=1 )
length$gene <- rownames(length)

countData = read.table('SRR941816.out', header=TRUE, sep="\t", row.names=1 )
countData <- countData['SRR941816.bam']
countData$before_1_rep <- countData$SRR941816.bam
countData$SRR941816.bam <- NULL

before <- read.table('SRR941817.out', header = TRUE, sep = '\t', row.names = 1)
countData$before_2_rep <- before$SRR941817.bam

after_one <- read.table('SRR941818.out', header = TRUE, sep = '\t', row.names = 1)
countData$after_1_rep <- after_one$SRR941818.bam

after_two <- read.table('SRR941819.out', header = TRUE, sep = '\t', row.names = 1)
countData$after_2_rep <- after_two$SRR941819.bam


cond_1 = rep("0 min", 2)
cond_2 = rep("30 min", 2)

samples = names(countData)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)
colData

dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
dds$condition = relevel(dds$condition,"0 min")

dds = DESeq(dds)
#scale factors:
dds$sizeFactor

res = results(dds)
sorted = res[with(res, order(padj, -log2FoldChange)), ]
sorted.df = data.frame("id"=rownames(sorted),sorted)
write.table(sorted.df, file="result.txt", sep="\t", col.names=NA, quote=FALSE)


#Только при TRUE будут посчитаны нормализованные значения
nc = counts(dds,normalized=TRUE)
dt = data.frame("id"=rownames(nc),nc)
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

plotMA(res)
View(as.data.frame(res))

#rld <- rlog(dds)
#plotPCA(rld) + theme_bw() + ggtitle('PCA plot with rlog')
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),colData=colData(dds))
plotPCA( DESeqTransform( se ) ) + theme_bw() + ggtitle("PCA plot with DESeq Transform")

View(sorted.df)
sorted.df <- na.omit(sorted.df)
sorted.df$Expression_shift <- 'Genes without significant changes'
sorted.df[sorted.df$log2FoldChange > 1 & sorted.df$padj < 0.05,]$Expression_shift <- 'Up, log(FC) > 1, adj p < 0.05'
sorted.df[sorted.df$log2FoldChange < -1 & sorted.df$padj < 0.05,]$Expression_shift <- 'Down, log(FC) < -1, adj p < 0.05'
sorted.df <- sorted.df[log(sorted.df$baseMean) > 0,]

nrow(sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange >= 1,])
nrow(sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange <= -1,])
nrow(sorted.df[sorted.df$padj < 0.05 & abs(sorted.df$log2FoldChange) >= 3,])

library(ggthemes)
library(RColorBrewer)
library(ggplot2)
ggplot(data = sorted.df, aes(x = log2FoldChange, y = -log10(padj), col = Expression_shift)) + geom_point() + theme_bw() + 
  ggtitle("Plot which represents log2(FC) and -log10(adj. p-value) for each gene") +
  ylab('-log10 (Adjusted p-value)')  +
  xlab('log2 (FoldChange)') +
  scale_color_calc()

??scale_color_calc
#ggplot(data = sorted.df, aes(x = log(baseMean), y = log2FoldChange)) + geom_point() + theme_bw()

nrow(sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange >= 2,])
nrow(sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange <= -2,])
nrow(sorted.df[sorted.df$padj < 0.05 & abs(sorted.df$log2FoldChange) >= 2,])

sorted.df[sorted.df$padj < 0.05 & abs(sorted.df$log2FoldChange) >= 2,]

up <- data.frame('up' = sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange >= 1,]$id)
up$up <- unlist(lapply(strsplit(as.character(up$up), '-'), function(x)x[2]))
up[50,]
write.csv(up, 'up_1.csv', quote = F, row.names = F)
down <- data.frame('down' = sorted.df[sorted.df$padj < 0.05 & sorted.df$log2FoldChange <= -1,])
down$down <-  unlist(lapply(strsplit(as.character(down$down), '-'), function(x)x[2]))
down[51,]
write.csv(down, 'down_1.csv', quote = F, row.names = F)



#heatmap
gene = unlist(lapply(strsplit(as.character(dt[,1]), '-'), function(x)x[2]))
vals = as.matrix(dt[,2:ncol(dt)])
vals = jitter(vals, factor = 1, amount=0.00001)
score = NULL
for (i in 1:nrow(vals)) {
  row=vals[i,]
  zscore=(row-mean(row))/sd(row)
  score =rbind(score,zscore)
}
row.names(score) = gene
zscore=score
mat = as.matrix(zscore)

#pdf('heatmap.pdf')
colors = colorRampPalette(c("blue","white","red"))
heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(14,14),lhei=c(1,5), theme(axis.text.x=element_text(angle = 17 , hjust = 1, size = 9)))
#invisible(dev.off())


#corr between FC and gene_length?

match[as.character(length$gene), sorted.df]
match[sorted.df$id, rownames(length)]
