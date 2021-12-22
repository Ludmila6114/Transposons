#This script correspond to Kelleheer call, may 2021.
#Insertions, their chipseq and piRNA status:

data <- readxl::read_xlsx('/mnt/raid/pro_milka/2021_year/Kel_redo/assTEv1_160_insertions_final.xlsx')
TE <- GRanges(data$Chromosome, IRanges(data$Begin, data$End), strand = data$Strand, Element = data$Element, id = data$id, percent = data$`% of canonical copy`, chromatin = data$Chromatin)
TE$full <- ifelse(TE$percent > 0.8, 'FULL', 'DIV, <90% length')
TE_1000_longer = TE
start(TE_1000_longer) = start(TE) - 1000
end(TE_1000_longer) = end(TE) + 1000

#de

regions <- regions[regions$which_up != FALSE | regions$which_down != FALSE,]
regions <- GRanges(seqnames = regions$V1, ranges = IRanges(start = regions$V2, end = regions$V3), log2FC = regions$logFC, up = regions$which_up, down = regions$which_down)
table(countOverlaps(TE_1000_longer, regions))

#У 19 инсерций что-то статистически значимо поменялось по метке.

#ужно записать эти результаты в инсерции TE, соотнести с геном и его экспрессией.




#END DE

chip_160_dysgenic_rep1 <- read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_160_genome_rep1/F1_9_160_H3K9me3_vs_Control_peaks.narrowPeak', sep = '\t', header = FALSE)
chip_160_dysgenic_rep2 <- read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_160_rep2/F1_9_160_H3K9me3_vs_Control_REP2_peaks.narrowPeak', sep = '\t', header = FALSE)
chip_160_dysgenic_rep1 <- GRanges(seqnames = chip_160_dysgenic_rep1$V1, ranges = IRanges(start = chip_160_dysgenic_rep1$V2, end = chip_160_dysgenic_rep1$V3), enrichment = chip_160_dysgenic_rep1$V7)
chip_160_dysgenic_rep2 <- GRanges(seqnames = chip_160_dysgenic_rep2$V1, ranges = IRanges(start = chip_160_dysgenic_rep2$V2, end = chip_160_dysgenic_rep2$V3), enrichment = chip_160_dysgenic_rep2$V7)

overlaps = findOverlaps(TE_1000_longer, chip_160_dysgenic_rep1)
overlaps2 = findOverlaps(TE_1000_longer, chip_160_dysgenic_rep2)
TE$Chip_dys <- rep(0, length(TE))

for (i in 1:length(overlaps)){
  print(i)
  TE[queryHits(overlaps[i]),]$Chip_dys <- chip_160_dysgenic_rep1[subjectHits(overlaps[i]),]$enrichment + TE[queryHits(overlaps[i]),]$Chip_dys
}
for (i in 1:length(overlaps2)){
  print(i)
  TE[queryHits(overlaps2[i]),]$Chip_dys <- chip_160_dysgenic_rep2[subjectHits(overlaps2[i]),]$enrichment + TE[queryHits(overlaps2[i]),]$Chip_dys
}

TE$Chip_dys <- TE$Chip_dys/2

chip_160_rec_rep1 <- read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_and_align_to_160_genome_rep1/F1_160_9_H3K9me3_vs_Control_peaks.narrowPeak', sep = '\t', header = FALSE)
chip_160_rec_rep2 <- read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/trimmed_160_rep2/F1_160_9_H3K9me3_vs_Control_REP2_peaks.narrowPeak', sep = '\t', header = FALSE)
chip_160_rec_rep1 <- GRanges(seqnames = chip_160_rec_rep1$V1, ranges = IRanges(start = chip_160_rec_rep1$V2, end = chip_160_rec_rep1$V3), enrichment = chip_160_rec_rep1$V7)
chip_160_rec_rep2 <- GRanges(seqnames = chip_160_rec_rep2$V1, ranges = IRanges(start = chip_160_rec_rep2$V2, end = chip_160_rec_rep2$V3), enrichment = chip_160_rec_rep2$V7)
over_rec = findOverlaps(TE_1000_longer, chip_160_rec_rep1)
over_r2 = findOverlaps(TE_1000_longer, chip_160_rec_rep2)
TE$Chip_rec <- rep(0, length(TE))

for (i in 1:length(over_rec)){
  print(i)
  TE[queryHits(over_rec[i]),]$Chip_rec <- chip_160_rec_rep1[subjectHits(over_rec[i]),]$enrichment + TE[queryHits(over_rec[i]),]$Chip_rec
}
for (i in 1:length(over_r2)){
  print(i)
  TE[queryHits(over_r2[i]),]$Chip_rec <- chip_160_rec_rep2[subjectHits(over_r2[i]),]$enrichment + TE[queryHits(over_r2[i]),]$Chip_rec
}

TE$Chip_rec <- TE$Chip_rec/2
View(as.data.frame(TE))

ggplot(data = as.data.frame(TE), aes(x = Chip_rec, y = Chip_dys)) + geom_point() + theme_bw() + 
  geom_smooth(method = 'lm') + ggtitle('Total H3K9me3 enrichment within 1Kb flanks\nfor each insertion between rec and dys samples')



.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

setwd('~/2021_year/TE_in_genomes_analysis/piRNA_insertions/')
bam_to_df <- function(path){
  bam <- scanBam(path)
  bam_field <- names(bam[[1]])
  list_bam <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list_bam)
  names(bam_df) <- bam_field
  bam_df <- data.frame(bam_df)
  bam_df$type <- ifelse(bam_df$qwidth == 21 | bam_df$qwidth == 22, 'siRNA', 
                        ifelse((bam_df$qwidth >= 23 & bam_df$qwidth <= 30), 'piRNA', 'smth'))
  bam_df <- bam_df[bam_df$flag != 4,]
  return(bam_df)
}

#160

bam_72 <- bam_to_df('SRR2096072.piRNA.fasta.__160.bam')
bam_79 <- bam_to_df('SRR2096079.piRNA.fasta.__160.bam')
bam_80 <- bam_to_df('SRR2096080.piRNA.fasta.__160.bam')
bam_81 <- bam_to_df('SRR2096081.piRNA.fasta.__160.bam')

piRNA_rec = GRanges(seqnames = bam_80$rname, strand = bam_80$strand, ranges = IRanges(start=bam_80$pos, width=bam_80$qwidth), type = bam_80$type)
piRNA_dis = GRanges(seqnames = bam_72$rname, strand = bam_72$strand, ranges = IRanges(start=bam_72$pos, width=bam_72$qwidth), type = bam_72$type)
piRNA_rec_2 = GRanges(seqnames = bam_81$rname, strand = bam_81$strand, ranges = IRanges(start=bam_81$pos, width=bam_81$qwidth), type = bam_81$type)
piRNA_dis_2 = GRanges(seqnames = bam_79$rname, strand = bam_79$strand, ranges = IRanges(start=bam_79$pos, width=bam_79$qwidth), type = bam_79$type)

piRNA_rec <- piRNA_rec[piRNA_rec$type == 'piRNA',]
piRNA_dis <- piRNA_dis[piRNA_dis$type == 'piRNA',]
piRNA_rec_2 <- piRNA_rec_2[piRNA_rec_2$type == 'piRNA',]
piRNA_dis_2 <- piRNA_dis_2[piRNA_dis_2$type == 'piRNA',]

TE$rec_piRNA_BODY <- (countOverlaps(TE, piRNA_rec, ignore.strand = TRUE)/226506 + countOverlaps(TE, piRNA_rec_2, ignore.strand = TRUE)/286181)/2  * 1000000
TE$dis_piRNA_BODY <- (countOverlaps(TE, piRNA_dis, ignore.strand = TRUE)/133438 + countOverlaps(TE, piRNA_dis_2, ignore.strand = TRUE)/119349)/2 * 1000000

TE$rec_piRNA_ONLY_FLANKS <- (countOverlaps(TE_1000_longer, piRNA_rec, ignore.strand = TRUE)/226506 + countOverlaps(TE_1000_longer, piRNA_rec_2, ignore.strand = TRUE)/286181)/2 *1000000 - TE$rec_piRNA_BODY
TE$dis_piRNA_ONLY_FLANKS <- (countOverlaps(TE_1000_longer, piRNA_dis, ignore.strand = TRUE)/133438 + countOverlaps(TE_1000_longer, piRNA_dis_2, ignore.strand = TRUE)/119349)/2 *1000000- TE$dis_piRNA_BODY

View(as.data.frame(TE))
#write.csv(as.data.frame(TE), '~/2021_year/TE_in_genomes_analysis/FINAL_TABLE_315_INS_160_CHIP_and_PIRNA.csv', quote = F, row.names = F)
#diff = as.data.frame(TE)
#diff = diff[diff$rec_piRNA_ONLY_FLANKS != diff$dis_piRNA_ONLY_FLANKS,]
#diff <- diff[,c('id', "rec_piRNA_ONLY_FLANKS", "dis_piRNA_ONLY_FLANKS")]
#diff$Dysgenic <- diff$dis_piRNA_ONLY_FLANKS
#diff$Reciprocal <- diff$rec_piRNA_ONLY_FLANKS
#diff$rec_piRNA_ONLY_FLANKS <- NULL
#diff$dis_piRNA_ONLY_FLANKS <- NULL
#df = melt(diff[61:119,])
#library('gplots')
#ggplot(data = df, aes(x=variable, y=id, fill=log2(value))) + 
#  geom_tile(color = "white") +
#  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    #                   midpoint = 5, space = "Lab",
    #                   name="log2(RPM)") + ggtitle('piRNA in flanking sequences for 160 strain\nlist2, filter: dysgenic value != reciprocal') +
  #theme_bw() + xlab('Drosophila virilis strain') + theme(axis.text.x = element_text(angle = 30 , hjust = 1, size = 10))

table(TE$rec_piRNA_ONLY_FLANKS)
table(TE$dis_piRNA_ONLY_FLANKS)
TE$piRNA_status <- NA

TE$piRNA_status <- ifelse(TE$rec_piRNA_ONLY_FLANKS == 0 & TE$dis_piRNA_ONLY_FLANKS == 0, 'both 0 piRNA',
                        ifelse(TE$rec_piRNA_ONLY_FLANKS == 0 & TE$dis_piRNA_ONLY_FLANKS != 0, 'Mostly dysgenic',
                          ifelse(TE$dis_piRNA_ONLY_FLANKS == 0 & TE$rec_piRNA_ONLY_FLANKS !=0, 'Mostly reciprocal',
                                 ifelse(TE$rec_piRNA_ONLY_FLANKS/TE$dis_piRNA_ONLY_FLANKS >=2, 'Mostly reciprocal',
                                      ifelse(TE$dis_piRNA_ONLY_FLANKS/TE$rec_piRNA_ONLY_FLANKS >=2, 'Mostly dysgenic',
                                               'Equal')))))


setwd('~/2021_year/Kelleher/F2_RNA_seq_ref/aln_to_ref/')
expressed_and_passed <- read.csv('DE_genes_3609.csv')
annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']
library(GenomicRanges)
library(RColorBrewer)


#chi-square, DE genes within 1Kb influence, inside also:
#GENES_INSIDE <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
#length(GENES_INSIDE)
#table(GENES_INSIDE$gene_id %in% expressed_and_passed$X)
#GENES_Outside <- Gene$gene_id[!Gene$gene_id %in% GENES_INSIDE$gene_id]
#table(GENES_Outside %in% expressed_and_passed$X)
#chi_inside = matrix(c(54, 67, 3555, 11361), nrow = 2, ncol = 2, byrow = TRUE)
#chisq.test(chi_inside)
#fisher.test(chi_inside)

#start(Gene) <- start(Gene) - 1000
#end(Gene) <- end(Gene) + 1000
#GENES_1k <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
#GENES_1k <- GENES_1k$gene_id[!GENES_1k$gene_id %in% GENES_INSIDE$gene_id]
#length(GENES_1k)
#table(GENES_1k %in% expressed_and_passed$X)
#GENES_Outside <- Gene$gene_id[(!Gene$gene_id %in% GENES_INSIDE$gene_id) & (!Gene$gene_id %in% GENES_1k)]
#table(GENES_Outside %in% expressed_and_passed$X)
#chi_inside = matrix(c(15, 35, 3540, 11326), nrow = 2, ncol = 2, byrow = TRUE)
#fisher.test(chi_inside)


Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
Gene <- Gene[!is.na(Nearest_ins),]
Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
nrow(as.data.frame(Gene))

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
Gene$closest_TE_strand <- strand(TE[Nearest_ins])
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

Gene$closest_TE_Chip_dys <- TE[Nearest_ins]$Chip_dys
Gene$closest_TE_Chip_rec <- TE[Nearest_ins]$Chip_rec
Gene$IS_DE <- ifelse(Gene$gene_id %in% expressed_and_passed$X, 'DE', 'not DE')
Gene$TE_start <- start(TE[Nearest_ins])
Gene$TE_end <- end(TE[Nearest_ins])
Gene$piRNA_status <- TE[Nearest_ins]$piRNA_status

Gene <- Gene[Gene$gene_id %in% GENES_INSIDE,]
Gene$exon_distance <- NA
for (i in 1:length(Gene)){
  print(i)
  Gene$exon_distance[i] <- ifelse(strand(Gene[i]) == '+' | strand(Gene[i]) != '-',
                                  min(abs(Gene[i]$TE_start - start(Gene[i])), abs(Gene[i]$TE_end - start(Gene[i]))),
                                  min(abs(end(Gene[i]) - Gene[i]$TE_end), abs(end(Gene[i]) - Gene[i]$TE_start)))
}

Gene$log2FC <- NA
for (i in 1:length(Gene)){
  print(i)
  Gene$log2FC[i] <- ifelse(Gene$gene_id[i] %in% expressed_and_passed$X,
                           expressed_and_passed[which(expressed_and_passed$id %in% Gene$gene_id[i]),]$log2FoldChange,'Nothing')
}
ggplot(data = as.data.frame(Gene[Gene$gene_id %in% GENES_INSIDE,]), aes(x = closest_TE_Chip_rec, y = closest_TE_Chip_dys, col = IS_DE)) +
  labs(col = 'Is a gene DE?') +
  geom_point(size = 4)  + ggtitle('Total enrichment of H3K9me3 tag in insertions\nand flanking sequences for insertions inside of a gene') +
  xlab('H3K9me3 tag in reciprocal sample, total enrichment') + 
  ylab('H3K9me3 tag in dysgenic sample, total enrichment') +
  theme_bw(base_size = 13) 

table(Gene[Gene$gene_id %in% GENES_INSIDE,]$IS_DE)
Gene$log2FC
ggplot(data = as.data.frame(Gene[Gene$gene_id %in% expressed_and_passed$X,]), 
       aes(x = closest_TE_Chip_rec, y = as.numeric(log2FC), col = closest_TE)) +
  geom_point(size = 3) + theme_bw(base_size = 13) +
  labs(col = 'TE inside a gene') + ggtitle('Each dot represents a DE gene with its log2FC (ref=reciprocal)\nand total enrichment of H3K9me3 tag for insertions inside a gene') +
  xlab('H3K9me3 tag in reciprocal sample, total enrichment') + ylab('log2FC')

ggplot(data = as.data.frame(Gene[Gene$gene_id %in% expressed_and_passed$X,]), 
       aes(x = closest_TE_Chip_dys, y = as.numeric(log2FC), col = closest_TE)) +
  geom_point(size = 3) + theme_bw(base_size = 13) +
  labs(col = 'TE inside a gene') + ggtitle('Each dot represents a DE gene with its log2FC (ref=reciprocal)\nand total enrichment of H3K9me3 tag for insertions inside a gene') +
  xlab('H3K9me3 tag in dysgenic sample, total enrichment') + ylab('log2FC')


length(Gene[Gene$gene_id %in% expressed_and_passed$X & Gene$log2FC > 0 & Gene$closest_TE_Chip_rec != 0,])
chi_inside = matrix(c(25, 14, 6, 9), nrow = 2, ncol = 2, byrow = TRUE)
fisher.test(chi_inside)

inside <- Gene[Gene$gene_id %in% GENES_INSIDE,]
View(as.data.frame(inside))

length(inside[strand(inside) != inside$closest_TE_strand & inside$gene_id %in% expressed_and_passed$X,])
length(inside[strand(inside) != inside$closest_TE_strand,])
chi_inside = matrix(c(29, 32, 25, 35), nrow = 2, ncol = 2, byrow = TRUE)
fisher.test(chi_inside)

View(as.data.frame(Gene))
length(Gene)
length(inside$closest_TE_strand)

table(Gene$piRNA_status, Gene$IS_DE)

Gene <- Gene[Gene$gene_id %in% expressed_and_passed$X]
Gene
d <- data.frame(gene_id = Gene[Gene$piRNA_status == 'Mostly reciprocal',]$gene_id, log2FC = Gene[Gene$piRNA_status == 'Mostly reciprocal',]$log2FC, label = 'both0')
d <- unique(d)
d
others = Gene[!Gene$gene_id %in% d$gene_id,]$log2FC
ot = data.frame(log2FC = others, label = 'NOT INSIDE')
ot
d$gene_id <- NULL
d = rbind(d, ot)
d
d$label <- as.factor(d$label)
d$log2FC <- as.numeric(as.character(d$log2FC))

Permtutation_test(d, 100000)


ggplot(data = as.data.frame(Gene[Gene$gene_id %in% expressed_and_passed$X,]), 
       aes(x = piRNA_status, y = as.numeric(log2FC), fill = piRNA_status)) + geom_boxplot() + theme_bw() + ylim(c(-1.5, +1.5)) +
  geom_jitter() + labs(fill = 'piRNA in insertion inside a gene') +
  ggtitle('Each dot represent a DE gene with an insertion inside and its piRNA "status"') +
  xlab('piRNA in flanking sequences for insertions') + ylab('log2FC')

#Genes_piRNA_plus <- Gene[Gene$closest_TE_PIRBA %in% c('Mostly dysgenic','Mostly reciprocal'),]$gene_id
#Genes_piRNA_minus <- Gene[Gene$closest_TE_PIRBA %in% c('both 0 piRNA','Equal'),]$gene_id
#table(Genes_piRNA_minus %in% expressed_and_passed$id)
#table(Genes_piRNA_plus %in% expressed_and_passed$id)
#chi_inside = matrix(c(1008, 3143, 2593, 8259), nrow = 2, ncol = 2, byrow = TRUE)
#fisher.test(chi_inside)
#Gene$closest_distance <- as.numeric(as.character(Gene$closest_distance))
#Gene <- Gene[Gene$closest_distance < 5000,]
#Gene <- Gene[Gene$gene_id %in% expressed_and_passed$X,]
#length(Gene)
#Gene <- Gene[Gene$gene_id %in% GENES_INSIDE,]

Gene <- Gene[Gene$gene_id %in% GENES_INSIDE,]
Gene$exon_distance <- NA
for (i in 1:length(Gene)){
  print(i)
  Gene$exon_distance[i] <- ifelse(strand(Gene[i]) == '+' | strand(Gene[i]) != '-',
                             min(abs(Gene[i]$TE_start - start(Gene[i])), abs(Gene[i]$TE_end - start(Gene[i]))),
                             min(abs(end(Gene[i]) - Gene[i]$TE_end), abs(end(Gene[i]) - Gene[i]$TE_start)))
}

Gene$log2FC <- NA
for (i in 1:length(Gene)){
  print(i)
  Gene$log2FC[i] <- ifelse(Gene$gene_id[i] %in% expressed_and_passed$X,
                           expressed_and_passed[which(expressed_and_passed$id %in% Gene$gene_id[i]),]$log2FoldChange,'Nothing')
}

ggplot(data = as.data.frame(Gene[Gene$gene_id %in% expressed_and_passed$X,]),
       aes(x = exon_distance, y = as.numeric(as.character(log2FC)), col = closest_TE_full)) + 
  geom_point(size = 4) + xlab('Is a gene DE?') + 
  ylab('distance to a gene TSS') + 
  ggtitle('Distance to TSS for all DE/not DE genes\nwith TE inside') + theme_bw() +
  facet_wrap(.~closest_TE, scales = 'free')

t.test(Gene$exon_distance[Gene$IS_DE == 'DE'], Gene$exon_distance[Gene$IS_DE == 'not DE'],)


#length(Gene[Gene$closest_TE_PIRBA == 'Mostly dysgenic' & Gene$log2FC < 0,])
#chi_inside = matrix(c(12, 21, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
#fisher.test(chi_inside)
#length(Gene)
#length(Gene[Gene$closest_TE_CHIP == 'Mostly reciprocal' & Gene$log2FC < 0,])
#chi_inside = matrix(c(9,2, 7, 1), nrow = 2, ncol = 2, byrow = TRUE)
#fisher.test(chi_inside)

ggplot(data = as.data.frame(Gene), aes(x = closest_distance, y = log2FC, col = closest_TE)) + 
  geom_point(size = 6) + theme_bw() + 
  xlim(c(0, 5000)) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle('Distance to the closest TE (lim: 5Kb)\nfrom DE gene with its log2FC\n103 Genes; y-axis limits: -1.5, +1.5') +
  ylim(c(-1.5, +1.5))


View(as.data.frame(Gene))
ggplot(data = as.data.frame(Gene), aes(x = closest_distance, y = log2FC, col = closest_TE, shape = closest_TE_CHIP)) + 
  geom_point(size = 4) + theme_bw() + 
  xlim(c(0, 5000)) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle('Distance to the closest TE (any length)\nfrom DE gene with its log2FC\n103 Genes: Helena(10), Paris(36), Penelope(23), Polyphemus(34)') +
  facet_wrap(.~closest_TE_PIRBA)

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
  ggtitle('Each dot represent a gene, \ngroups show the distance to the closest TE\nlimits: -1.5, +1.5 log2FC (ref = reciprocal)')+
  ylim(c(-1.5, 1.5))


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



d <- data.frame(gene_id = Gene[Gene$group == '0',]$gene_id, log2FC = Gene[Gene$group == '0',]$log2FC, label = 'INSIDE')
d <- unique(d)
others = annotation[!annotation$gene_id %in% d$gene_id,]$logFC
ot = data.frame(log2FC = others, label = 'NOT INSIDE')
d$gene_id <- NULL
d = rbind(d, ot)
d$label <- as.factor(d$label)

Permtutation_test(d, 100000)
length(GENES_INSIDE)






#######################################################################

ins9 <- readxl::read_xlsx('~/2021_year/assTEv1_9_insertions_final.xlsx')
ins9$Chromatin <- ifelse(ins9$Chromatin == 'Heterochromatic?', 'Heterochromatic', 
                         ifelse(ins9$Chromatin == 'Euchromatic',
                                'Euchromatic', 'Heterochromatic'))
ins9$Chromosome <- paste(ins9$Chromosome, '_', sep = '')
M9 <- GRanges(seqnames = ins9$Chromosome, strand= ins9$Strand, ranges = IRanges(start = ins9$Begin, end = ins9$End),
              Element = ins9$Element, id = ins9$id, percent = ins9$`% of canonical copy`,
              chromatin = ins9$Chromatin)
M9$full <-  ifelse(M9$percent > 0.9, 'FULL', 'DIV, <90% length')
#chip_9 <- read.table('~/2021_year/TE_in_genomes_analysis/Aravin_chip_analysis/DE_9/master_list.bed', sep = '\t', header = FALSE)
#chip_9 <- GRanges(seqnames = chip_9$V1, ranges = IRanges(start = chip_9$V2, end = chip_9$V3))
M9_1000_longer = M9
start(M9_1000_longer) = start(M9) - 1000
end(M9_1000_longer) = end(M9) + 1000
#M9$Chip <- countOverlaps(M9_1000_longer, chip_9, ignore.strand = TRUE)

length(M9)
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

setwd('~/2021_year/TE_in_genomes_analysis/piRNA_insertions/')
bam_to_df <- function(path){
  bam <- scanBam(path)
  bam_field <- names(bam[[1]])
  list_bam <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list_bam)
  names(bam_df) <- bam_field
  bam_df <- data.frame(bam_df)
  bam_df$type <- ifelse(bam_df$qwidth == 21 | bam_df$qwidth == 22, 'siRNA', 
                        ifelse((bam_df$qwidth >= 23 & bam_df$qwidth <= 30), 'piRNA', 'smth'))
  bam_df <- bam_df[bam_df$flag != 4,]
  return(bam_df)
}

#9

bam_72 <- bam_to_df('SRR2096072.piRNA.fasta.__9.bam')
bam_79 <- bam_to_df('SRR2096079.piRNA.fasta.__9.bam')
bam_80 <- bam_to_df('SRR2096080.piRNA.fasta.__9.bam')
bam_81 <- bam_to_df('SRR2096081.piRNA.fasta.__9.bam')

piRNA_rec = GRanges(seqnames = bam_80$rname, strand = bam_80$strand, ranges = IRanges(start=bam_80$pos, width=bam_80$qwidth), type = bam_80$type)
piRNA_dis = GRanges(seqnames = bam_72$rname, strand = bam_72$strand, ranges = IRanges(start=bam_72$pos, width=bam_72$qwidth), type = bam_72$type)
piRNA_rec_2 = GRanges(seqnames = bam_81$rname, strand = bam_81$strand, ranges = IRanges(start=bam_81$pos, width=bam_81$qwidth), type = bam_81$type)
piRNA_dis_2 = GRanges(seqnames = bam_79$rname, strand = bam_79$strand, ranges = IRanges(start=bam_79$pos, width=bam_79$qwidth), type = bam_79$type)

piRNA_rec <- piRNA_rec[piRNA_rec$type == 'piRNA',]
piRNA_dis <- piRNA_dis[piRNA_dis$type == 'piRNA',]
piRNA_rec_2 <- piRNA_rec_2[piRNA_rec_2$type == 'piRNA',]
piRNA_dis_2 <- piRNA_dis_2[piRNA_dis_2$type == 'piRNA',]

M9$rec_piRNA_BODY <- (countOverlaps(M9, piRNA_rec, ignore.strand = TRUE)/444618 + countOverlaps(M9, piRNA_rec_2, ignore.strand = TRUE)/564392)/2  * 1000000
M9$dis_piRNA_BODY <- (countOverlaps(M9, piRNA_dis, ignore.strand = TRUE)/240512 + countOverlaps(M9, piRNA_dis_2, ignore.strand = TRUE)/222382)/2 * 1000000

M9$rec_piRNA_ONLY_FLANKS <- (countOverlaps(M9_1000_longer, piRNA_rec, ignore.strand = TRUE)/444618 + countOverlaps(M9_1000_longer, piRNA_rec_2, ignore.strand = TRUE)/564392)/2 *1000000 - M9$rec_piRNA_BODY
M9$dis_piRNA_ONLY_FLANKS <- (countOverlaps(M9_1000_longer, piRNA_dis, ignore.strand = TRUE)/240512 + countOverlaps(M9_1000_longer, piRNA_dis_2, ignore.strand = TRUE)/222382)/2 *1000000- M9$dis_piRNA_BODY

View(as.data.frame(M9))
write.csv(as.data.frame(M9), '~/2021_year/TE_in_genomes_analysis/FINAL_TABLE_122_INS_9_CHIP_and_PIRNA.csv', quote = F, row.names = F)


















######################################DE
setwd('~/2021_year/Kelleher/F2_RNA_seq_ref/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11, 12, 13)]
cond_1 = rep("disgenic_young", 4)
cond_2 = rep("nondisgenic_young", 4)
samples = names(featurecounts)
condition = factor(c(cond_1, cond_2))
colData = data.frame(samples=samples, condition=condition)
colData
dds = DESeqDataSetFromMatrix(countData=featurecounts, colData=colData, design = ~condition)
dds$condition = relevel(dds$condition, ref = "nondisgenic_young")
dds = DESeq(dds)

#scale factors:
dds$sizeFactor
head(counts(dds))
head(counts(dds, normalized = TRUE))
head(featurecounts)
normalized = as.data.frame(counts(dds, normalized = TRUE))
normalized$gene_id <- rownames(normalized)
normalized$width <- NA
ref <- rtracklayer::import('~/virilis_fly_base_1.07/dvir-all-r1.07.gtf')
ref <- ref[ref$type == 'gene',]
for (i in 1:nrow(normalized)){
  normalized$width[i] <- ref[which(ref$gene_id %in% normalized$gene_id[i]),]@ranges@width
  print(i)
}

sorted = res[with(res, order(padj, -log2FoldChange)), ]
sorted.df = data.frame("id"=rownames(sorted),sorted)
sorted.df <- na.omit(sorted.df)
sorted.df <- sorted.df[sorted.df$padj < 0.05,]
#PCA
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
plotPCA( DESeqTransform( se ) )
#write.csv(sorted.df, 'DE_genes_3609.csv')


GO <- sorted.df
#Посмотрим, кто это, по ортологам dmel.
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(GO$id, ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym
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

ego <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(ego, showCategory=30)
nrow(sorted.df[abs(sorted.df$log2FoldChange) > 0.5,])


#piRNA to 44 TE


setwd('~/2021_year/Kelleher/F2_RNA_seq_ref/aln_to_ref/')
expressed_and_passed <- read.csv('DE_genes_3609.csv')
annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']


length(Result$gene_id)
length(unique(Result$gene_id))

A = as.data.frame(table(Result$group))
ggplot(data = A, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity", fill = 'royalblue2') + theme_bw() + coord_flip() + xlab('Distance') + ylab('Frequency')


