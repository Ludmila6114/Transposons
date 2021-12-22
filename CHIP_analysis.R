#160
#CHIP
library('edgeR')
library('csaw')
FDR_T <- 0.05
logFC_T <- log2(2)
setwd('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/Diff_Enrich_160_genome/')
#rep1: 1609, 9160 than rep2 1609 9160 only chip
bamFiles <- c('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_160_genome_rep1/SRR1533732_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_160_genome_rep1/SRR1533742_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_160_rep2/SRR1533733_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_160_rep2/SRR1533743_trimmed.fq.sorted.rmdup.UNIQ.bam'
              )
regions = read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/Diff_Enrich_160_genome/master.bed', sep = '\t', header = FALSE)


meta <- c()
meta[1] <- "master_list"
meta[2] <- "F1_160_9"
meta[3] <- "F1_9_160"
group <- rep(meta[2:3], 2)
group
group[group == meta[2]] <- 1
group[group == meta[3]] <- 0
group <- factor(group)
design <- model.matrix(~group)
colnames(design) <- c("intercept", "Type")

#"F1_160_9" = 1, TREATMENT?

ids <- paste(regions[,1], ":", regions[,2], "-",
             
             regions[,3], sep="")

d <- GRanges(regions[,1], IRanges(regions[,2],
                                  
                                  regions[,3]), id=ids)
dCounts <- regionCounts(bamFiles, regions=d)
data <- normFactors(dCounts, se.out=TRUE)
#QC plots
cpms_norm_N <- cpm(asDGEList(dCounts))
cpms_norm_Y <- cpm(asDGEList(data), normalized.lib.size=TRUE)
cpms_norm_N[cpms_norm_N == 0] <- min(cpms_norm_N[cpms_norm_N > 0])
cpms_norm_Y[cpms_norm_Y == 0] <- min(cpms_norm_Y[cpms_norm_Y > 0])
ylab <- "Expression Values"
nms <- rep("", 4)
boxplot(log2(assay(dCounts)+1),
        
        xlab="", ylab=ylab, las=2, names=nms, main="Raw")

boxplot(log2(cpms_norm_N),
        
        xlab="", ylab=ylab, las=2, names=nms, main="CPM")

boxplot(log2(cpms_norm_Y),
        
        xlab="", ylab=ylab, las=2, names=nms, main="CPM + TMM")

y <- asDGEList(data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)
FDR <- p.adjust(results$table$PValue, method = "BH")
out <- data.frame(chrom=regions[,1],
                  start=regions[,2],
                  end=regions[,3],
                  logCPM=results$table$logCPM,
                  logFC=results$table$logFC,
                  PValue=results$table$PValue,
                  FDR)
which_up <- out$logFC >= logFC_T & out$FDR <= FDR_T
which_down <- out$logFC <= -logFC_T & out$FDR <= FDR_T
table(which_up)
table(which_down)

View(as.data.frame(results))
table(which_up == which_down)

regions <- cbind(regions, as.data.frame(results))
regions <- cbind(regions, which_up)
regions <- cbind(regions, which_down)



############9
setwd('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/DE_9/')
#rep1: 1609, 9160 than rep2 1609 9160 only chip
bamFiles <- c('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_9_genome/SRR1533732_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_9_genome/SRR1533742_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_9_rep2/SRR1533733_trimmed.fq.sorted.rmdup.UNIQ.bam',
              '~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_9_rep2/SRR1533743_trimmed.fq.sorted.rmdup.UNIQ.bam'
)
regions = read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/DE_9/master_list.bed', sep = '\t', header = FALSE)
meta <- c()
meta[1] <- "master_list"
meta[2] <- "F1_160_9"
meta[3] <- "F1_9_160"
group <- rep(meta[2:3], 2)
group
group[group == meta[2]] <- 1
group[group == meta[3]] <- 0
group <- factor(group)
design <- model.matrix(~group)
colnames(design) <- c("intercept", "Type")

ids <- paste(regions[,1], ":", regions[,2], "-",
             
             regions[,3], sep="")

d <- GRanges(regions[,1], IRanges(regions[,2],
                                  
                                  regions[,3]), id=ids)
dCounts <- regionCounts(bamFiles, regions=d)
data <- normFactors(dCounts, se.out=TRUE)
#QC plots
cpms_norm_N <- cpm(asDGEList(dCounts))
cpms_norm_Y <- cpm(asDGEList(data), normalized.lib.size=TRUE)
cpms_norm_N[cpms_norm_N == 0] <- min(cpms_norm_N[cpms_norm_N > 0])
cpms_norm_Y[cpms_norm_Y == 0] <- min(cpms_norm_Y[cpms_norm_Y > 0])
ylab <- "Expression Values"
nms <- rep("", 4)
boxplot(log2(assay(dCounts)+1),
        
        xlab="", ylab=ylab, las=2, names=nms, main="Raw")

boxplot(log2(cpms_norm_N),
        
        xlab="", ylab=ylab, las=2, names=nms, main="CPM")

boxplot(log2(cpms_norm_Y),
        
        xlab="", ylab=ylab, las=2, names=nms, main="CPM + TMM")

y <- asDGEList(data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)
FDR <- p.adjust(results$table$PValue, method = "BH")
out <- data.frame(chrom=regions[,1],
                  start=regions[,2],
                  end=regions[,3],
                  logCPM=results$table$logCPM,
                  logFC=results$table$logFC,
                  PValue=results$table$PValue,
                  FDR)
which_up <- out$logFC >= logFC_T & out$FDR <= FDR_T
which_down <- out$logFC <= -logFC_T & out$FDR <= FDR_T
table(which_up)
table(which_down)


######################################################################

setwd('~/2021_year/Kelleher/F2_RNA_seq_ref/aln_to_ref/')
expressed_and_passed <- read.csv('DE_genes_3609.csv')
annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']
library(GenomicRanges)
library(RColorBrewer)


#chi-square, DE genes within 1Kb influence, inside also:
GENES_INSIDE <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
length(GENES_INSIDE)
table(GENES_INSIDE$gene_id %in% expressed_and_passed$X)
GENES_Outside <- Gene$gene_id[!Gene$gene_id %in% GENES_INSIDE$gene_id]
table(GENES_Outside %in% expressed_and_passed$X)
chi_inside = matrix(c(54, 67, 3555, 11361), nrow = 2, ncol = 2, byrow = TRUE)
chisq.test(chi_inside)


start(Gene) <- start(Gene) - 10000
end(Gene) <- end(Gene) + 10000

GENES_1k <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
GENES_1k <- GENES_1k$gene_id[!GENES_1k$gene_id %in% GENES_INSIDE$gene_id]
table(GENES_1k %in% expressed_and_passed$X)

GENES_Outside <- Gene$gene_id[(!Gene$gene_id %in% GENES_INSIDE$gene_id) & (!Gene$gene_id %in% GENES_1k)]
table(GENES_Outside %in% expressed_and_passed$X)
chi_inside = matrix(c(94, 269, 3641, 11092), nrow = 2, ncol = 2, byrow = TRUE)
fisher.test(chi_inside)



Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
Gene <- Gene[!is.na(Nearest_ins),]
Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
nrow(as.data.frame(Gene))


Gene <- Gene[Gene$gene_id %in% expressed_and_passed$X,]

GENES_INSIDE <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
GENES_INSIDE <- GENES_INSIDE$gene_id
GENES_WITH_STRAND <- subsetByOverlaps(Gene, TE, ignore.strand = FALSE)
GENES_WITH_STRAND <- GENES_WITH_STRAND$gene_id


Gene$closest_TE <- TE[Nearest_ins]$Element
Gene$closest_TE_chr <- TE[Nearest_ins]$chromatin
Gene$closest_TE_per <- TE[Nearest_ins]$percent
Gene$closest_TE_full <- TE[Nearest_ins]$full
Gene$TE_id <- TE[Nearest_ins]$id
Gene$closest_distance <- NA
Gene$closest_distance <- ifelse(Gene$gene_id %in% GENES_INSIDE, 0, 
                                ifelse(
                                  start(TE[Nearest_ins]) > end(Gene) &
                                    start(TE[Nearest_ins]) > start(Gene), 
                                  start(TE[Nearest_ins]) - end(Gene),
                                  ifelse(
                                    end(TE[Nearest_ins]) < end(Gene) &
                                      end(TE[Nearest_ins]) < start(Gene), 
                                    start(Gene) - end(TE[Nearest_ins]),
                                    'Uncorrect')))    

table(Gene[Gene$closest_distance < 1000,]$gene_id %in% expressed_and_passed$X)
table(Gene$closest_distance)

Gene[Gene$closest_distance < 1000,]$gene_id %in% GENES_INSIDE

Gene$closest_TE_piRNA_status <- TE[Nearest_ins,]$piRNA_summary
View(as.data.frame(Gene))

Gene$closest_distance <- as.numeric(as.character(Gene$closest_distance))
Gene <- Gene[Gene$closest_distance < 5000,]
length(Gene)

Gene$log2FC <- NA
for (i in 1:length(Gene)){
  Gene$log2FC[i] <- ifelse(Gene$gene_id[i] %in% expressed_and_passed$X,
                           expressed_and_passed[which(expressed_and_passed$id %in% Gene$gene_id[i]),]$log2FoldChange,'Nothing')
}

View(as.data.frame(Gene))

ggplot(data = as.data.frame(Gene), aes(x = closest_distance, fill = closest_TE)) + 
  geom_histogram() + theme_bw() + scale_fill_brewer(palette = "Paired") + facet_wrap(.~closest_TE_full) +
  scale_y_discrete(limits = c(0, 1, 2, 5, 10, 15, 20, 25))

ggplot(data = as.data.frame(Gene), aes(x = closest_distance, y = log2FC, col = closest_TE)) + 
  geom_point(size = 6) + theme_bw() + 
  xlim(c(0, 5000)) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle('Distance to the closest TE (any length)\nfrom DE gene with its log2FC\n103 Genes; y-axis limits: -1.5, +1.5') +
  ylim(c(-1.5, +1.5))

ggplot(data = as.data.frame(Gene), aes(x = closest_distance, y = log2FC, col = closest_TE)) + 
  geom_point(size = 4) + theme_bw() + 
  xlim(c(0, 5000)) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle('Distance to the closest TE (any length)\nfrom DE gene with its log2FC\n103 Genes: Helena(10), Paris(36), Penelope(23), Polyphemus(34)') + facet_wrap(.~closest_TE)

table(Gene$closest_TE)

Gene$group <- NA
Gene$group <- ifelse(Gene$closest_distance == 0, '0',
                     ifelse(Gene$closest_distance < 1000, '0-1Kb',
                            ifelse(Gene$closest_distance < 2000, '1-2Kb',
                                   ifelse(Gene$closest_distance < 3000, '2-3Kb',
                                          ifelse(Gene$closest_distance < 4000, '3-4Kb',
                                                 '4-5Kb')))))

table(Gene$group)
ggplot(data = as.data.frame(Gene), aes(x = group, y = log2FC, fill = group)) + geom_boxplot() + geom_jitter() + 
  theme_bw() +
  ggtitle('Each dot represent a gene, \ngroups show the distance to the closest TE\nlimits: -1.5, +1.5 log2FC') +
  ylim(c(-1.5, 1.5)) +  facet_wrap(.~closest_TE)



length(Gene)


Permtutation_test <- function(d, P){
  levels(d$label)
  print('Obserbed test-statistics:')
  Test_stat1 <- abs(mean(d[d$label == levels(d$label)[1],]$log2FC - mean(d[d$label == levels(d$label)[2],]$log2FC)))
  print(Test_stat1)
  Test_stat2 <- abs(median(d[d$label == levels(d$label)[1],]$log2FC - median(d[d$label == levels(d$label)[2],]$log2FC)))
  print(Test_stat2)
  
  set.seed(1999)
  n = nrow(d)
  var = d$log2FC
  Perm_matrix <- matrix(0, nrow = n, ncol = P)
  for (i in 1:P){
    Perm_matrix[,i] <- sample(var, size = n, replace = FALSE)
  }
  
  
  Perm_test_stat1 <- rep(0, P)
  Perm_test_stat2 <- rep(0, P)
  
  for (i in 1:P){
    Perm_test_stat1[i] <- abs(mean(Perm_matrix[d$label ==levels(d$label)[1], i]) - mean(Perm_matrix[d$label == levels(d$label)[2], i]) )
    Perm_test_stat2[i] <- abs(median(Perm_matrix[d$label == levels(d$label)[1], i]) - median(Perm_matrix[d$label == levels(d$label)[2], i]) )
  }
  
  stat1 = as.data.frame(table(Perm_test_stat1 > Test_stat1))
  stat2 = as.data.frame(table(Perm_test_stat2 > Test_stat2))
  
  print('p-value1:')
  print(stat1[stat1$Var == 'TRUE',]$Freq/P)
  print('p-value2:')
  print(stat2[stat2$Var == 'TRUE',]$Freq/P)
}

annotation <- annotation_Sasha[annotation_Sasha$type == 'gene']
annotation <- annotation[annotation$gene_id %in% expressed_and_passed$X,]
annotation$logFC <- NA
for (i in 1:length(annotation)){
  print(i)
  annotation$logFC[i] <- expressed_and_passed[which(expressed_and_passed$id %in% annotation$gene_id[i]),]$log2FoldChange
}



d <- data.frame(gene_id = Gene[Gene$group == '0' & Gene$closest_TE == 'Penelope' & Gene$gene_id %in% GENES_WITH_STRAND,]$gene_id, log2FC = Gene[Gene$group == '0' & Gene$closest_TE == 'Penelope' & Gene$gene_id %in% GENES_WITH_STRAND,]$log2FC, label = 'INSIDE')
d <- unique(d)
others = annotation[!annotation$gene_id %in% d$gene_id,]$logFC
ot = data.frame(log2FC = others, label = 'NOT INSIDE')
d$gene_id <- NULL
d = rbind(d, ot)
d$label <- as.factor(d$label)

Permtutation_test(d, 100000)
table(Gene$closest_TE_full)

View(as.data.frame(Gene))
ggplot(data = as.data.frame(Gene[Gene$group == '0']), aes(x = closest_TE_piRNA_status, y = log2FC, fill = closest_TE_piRNA_status)) + 
  geom_boxplot() + 
  geom_jitter() + theme_bw() +
  ylim(c(-1.5, 1.5)) +
  ggtitle('The level of gene expression with respect of\npiRNA status of the closest TE\nylim: -1.5, +1.5 (54 genes)')

length(Gene)
View(as.data.frame(Gene[Gene$group == '0' & Gene$closest_TE == 'Penelope' & Gene$gene_id %in% GENES_WITH_STRAND,]))



