data <- readxl::read_xlsx('/mnt/raid/pro_milka/2021_year/Teflon_insertions/assTEv1_160_insertions_final.xlsx')
TE <- GRanges(data$Chromosome, IRanges(data$Begin, data$End), strand = data$Strand, Element = data$Element, 
              id = data$id, percent = data$`% of canonical copy`, chromatin = data$Chromatin)
TE
genome_9 <- readDNAStringSet('~/9_POLISHED/9_correct_names.fasta')
genome_160 <- readDNAStringSet('~/160_chromosome_assembly/160_chromosome.fasta')
names(genome_160) <- unlist(lapply(strsplit(names(genome_160), ' '), function(x)x[1]))

TE$flanc_left <- IRanges( start = start(ranges(TE))-3000, end=start(ranges(TE)))
TE$flanc_rigth <- IRanges( start = end(ranges(TE)), end=end(ranges(TE))+3000)
TE_left <- GRanges(seqnames = seqnames(TE), strand = strand(TE), ranges=TE$flanc_left)
TE_left
TE$left_fl_seq <- getSeq(genome_160, TE_left)
TE_rigth <- GRanges(seqnames = seqnames(TE), strand = strand(TE), ranges=TE$flanc_rigth)
TE$rigth_fl_seq <- getSeq(genome_160, TE_rigth)
TE$dist <- NA
count = 0
count_0 = vector()
count_4 = 0
for (N in 1:length(TE)){
  left_flanc <- TE$left_fl_seq[N]
  names(left_flanc) <- paste(N, ' left_flanc_seq')
  writeXStringSet(left_flanc, '~/left_flanc.fasta')
  left_blast <- system('blastn -task megablast -query ~/left_flanc.fasta -db ~/9_POLISHED/9_correct_names.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  left_df <- as.data.frame(do.call('rbind', strsplit(left_blast, '\t')))
  
  right_flanc <- TE$rigth_fl_seq[N]
  names(right_flanc) <- paste(N, ' rigth_flanc_seq')
  writeXStringSet(right_flanc, '~/right_flanc.fasta')
  right_blast <- system('blastn -task megablast -query ~/right_flanc.fasta -db ~/9_POLISHED/9_correct_names.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  right_df <- as.data.frame(do.call('rbind', strsplit(right_blast, '\t')))
  #################
  left_df <- left_df[as.numeric(as.character(left_df$V4)) > 2000,]
  right_df <- right_df[as.numeric(as.character(right_df$V4)) > 2000,]
  #print(paste(N, ' ', nrow(left_df), ' ', nrow(right_df)))
  if (N %% 10 == 0){
    print(paste(N, ' insertion done'))
  }
  if(nrow(left_df) == 1 & nrow(right_df) == 1){
    if (as.character(left_df$V2) == as.character(right_df$V2)){
      right_df$V9 <- as.numeric(as.character(right_df$V9))
      right_df$V10 <- as.numeric(as.character(right_df$V10))
      left_df$V9 <- as.numeric(as.character(left_df$V9))
      left_df$V10 <- as.numeric(as.character(left_df$V10))
      left = IRanges(min(left_df$V9, left_df$V10), max(left_df$V9, left_df$V10))
      right = IRanges(min(right_df$V9, right_df$V10), max(right_df$V9, right_df$V10))
      if (length(findOverlaps(left, right)) != 0){
        TE$dist[N] <- 0
      }
      else{
        TE$dist[N] = min(abs(start(right) - end(left)), abs(start(left) - end(right)))
        new_blast = GRanges(seqnames = left_df$V2, ranges = IRanges(start = min(end(right), end(left)), end = max(start(left), start(right))))
        strange = getSeq(genome_9, new_blast)
        names(strange) <- paste(N, 'strange')
        writeXStringSet(strange, '~/strange.fasta')
        blast_strange <- system('blastn -task megablast -query ~/strange.fasta -db ~/TE_bases/44_TE/dvir_full-size_TEs.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
        blast_df <- as.data.frame(do.call('rbind', strsplit(blast_strange, '\t')))
        if (nrow(blast_df) == 0){
          print(paste('nothing in 44 TE with length:', width(strange)))
          count_0 = c(count_0, width(strange))
          if (width(strange) > 1000){
            setwd('~/')
            writeXStringSet(strange, paste(width(strange), '_TE_STRANGE.fasta'))
          }
        }
        else{
        count_4 = count_4 + 1
        blast_df$TE <- TE[N]$Element
        print(paste(blast_df$V2, ' IN 160:', TE[N]$Element))
        }
      }
    }
  }
}


table(is.na(TE$dist))
#151 всего
#46 запапились на наши 4 (30%)
#37 полный 0 (24%)
#68 незамапились на базу с разной длиной (45%)
#Из них 6 имеют длину больше 1000

TE <- TE[!is.na(TE$dist),]
TE <- TE[TE$dist < 20,]
length(TE)
table(TE$Element)

###################################slide2: piRNA and 44 TE F0
bam_analysis <- function(path, sample_id, RPM){
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam <- scanBam(path)
  print('bam file was downloaded')
  bam_field <- names(bam[[1]])
  list_bam <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list_bam)
  names(bam_df) <- bam_field
  bam_df <- data.frame(bam_df)
  bam_df$type <- ifelse(bam_df$qwidth == 21 | bam_df$qwidth == 22, 'siRNA', 
                        ifelse((bam_df$qwidth >= 23 & bam_df$qwidth <= 30), 'piRNA', 
                               ifelse(bam_df$qwidth >= 31, 'Ultra_long', 'Ultra_short')))
  print(table(bam_df$type))
  
  result <- data.frame(TE = unique(bam_df$rname))
  
  arr <- function(my_data, name){
    result$piRNA <- NA
    result$siRNA <- NA
    result$Ultra_long <- NA
    result$Ultra_short <- NA
    
    
    for (i in 1:nrow(result)){
      if (result$TE[i] %in% my_data$Var1){
      result$piRNA[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'piRNA')),]$Freq
      #result$siRNA[i] <-  my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'siRNA')),]$Freq 
      #result$Ultra_long[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'Ultra_long')),]$Freq
      #result$Ultra_short[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'Ultra_short')),]$Freq
      }
    }
    
    return(result)
  }
  
  data <- as.data.frame(table(bam_df$rname, bam_df$type))
  result <- arr(data, sample_id)
  result$strain <- sample_id
  result$piRNA <- result$piRNA * 1000000/RPM
  #result$siRNA <- result$siRNA * 1000000/RPM
  #result$Ultra_long <- result$Ultra_long * 1000000/RPM
  #result$Ultra_short <- result$Ultra_short * 1000000/RPM
  print('ALL')
  return(result)
}

setwd('~/2021_year/GENES_piRNA_automatic/trimmed/Unmapped_fasta/align_piRNA_to_44_TE/')
base44_26 <- bam_analysis('SRR2096026_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096026', 3274511)
base44_38 <- bam_analysis('SRR2096038_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096038', 3456429)
base44_50 <- bam_analysis('SRR2096050_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096050', 1935610)
base44_61 <- bam_analysis('SRR2096061_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096061', 1542255)

base44_26 <- base44_26[order(base44_26$TE),]
base44_38 <- base44_38[order(base44_38$TE),]
base44_50 <- base44_50[order(base44_50$TE),]
base44_61 <- base44_61[order(base44_61$TE),]

F0_piRNA = data.frame(TE = base44_26$TE, P_160 = (base44_50$piRNA + base44_61$piRNA)/2, M_9 = (base44_26$piRNA + base44_38$piRNA)/2)
ggplot(data = F0_piRNA, aes(x = log2(P_160), y = log2(M_9))) + geom_point() + theme_bw() +
  geom_text(aes(label=unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),hjust=0, vjust=0) +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('Total amount of piRNA between F0 Embryos 0-2h, RPM') +
  xlab('log2(piRNA from 160 Embryo, RPM)') +
  ylab('log2(piRNA from 9 Embryo, RPM)')


#############################################mRNA F0 to 44TE

setwd('~/2021_year/raw_RNA_F0_JB/trimmed/aln_to_44/')





###########################################################slide3

TE_1000_longer = TE
start(TE_1000_longer) = start(TE) - 1000
end(TE_1000_longer) = end(TE) + 1000
#Правда ли, что суммарная пиРНК по фланкам уникальных инсерций в 160:
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


bam_72 <- bam_to_df('SRR2096072.piRNA.fasta.__160.bam')
bam_79 <- bam_to_df('SRR2096079.piRNA.fasta.__160.bam')
bam_80 <- bam_to_df('SRR2096080.piRNA.fasta.__160.bam')
bam_81 <- bam_to_df('SRR2096081.piRNA.fasta.__160.bam')
bam_26 <- bam_to_df('check_in_parents_data_from_GENESPIRNA_auto/SRR2096026_piRNA_to_160.bam')
bam_38 <- bam_to_df('check_in_parents_data_from_GENESPIRNA_auto/SRR2096038_piRNA_to_160.bam')
bam_50 <- bam_to_df('check_in_parents_data_from_GENESPIRNA_auto/SRR2096050_piRNA_to_160.bam')
bam_61 <- bam_to_df('check_in_parents_data_from_GENESPIRNA_auto/SRR2096061_piRNA_to_160.bam')


piRNA_9 <- GRanges(seqnames = bam_26$rname, strand = bam_26$strand, ranges = IRanges(start = bam_26$pos, width = bam_26$qwidth))
piRNA_9_2 <- GRanges(seqnames = bam_38$rname, strand = bam_38$strand, ranges = IRanges(start = bam_38$pos, width = bam_38$qwidth))
piRNA_160 <- GRanges(seqnames = bam_50$rname, strand = bam_50$strand, ranges = IRanges(start = bam_50$pos, width = bam_50$qwidth))
piRNA_160_2 <- GRanges(seqnames = bam_61$rname, strand = bam_61$strand, ranges = IRanges(start = bam_61$pos, width = bam_61$qwidth))

piRNA_rec = GRanges(seqnames = bam_80$rname, strand = bam_80$strand, ranges = IRanges(start=bam_80$pos, width=bam_80$qwidth), type = bam_80$type)
piRNA_dis = GRanges(seqnames = bam_72$rname, strand = bam_72$strand, ranges = IRanges(start=bam_72$pos, width=bam_72$qwidth), type = bam_72$type)
piRNA_rec_2 = GRanges(seqnames = bam_81$rname, strand = bam_81$strand, ranges = IRanges(start=bam_81$pos, width=bam_81$qwidth), type = bam_81$type)
piRNA_dis_2 = GRanges(seqnames = bam_79$rname, strand = bam_79$strand, ranges = IRanges(start=bam_79$pos, width=bam_79$qwidth), type = bam_79$type)

TE$rec_piRNA_BODY <- (countOverlaps(TE, piRNA_rec, ignore.strand = TRUE)/226506 + countOverlaps(TE, piRNA_rec_2, ignore.strand = TRUE)/286181)/2  * 1000000
TE$dis_piRNA_BODY <- (countOverlaps(TE, piRNA_dis, ignore.strand = TRUE)/133438 + countOverlaps(TE, piRNA_dis_2, ignore.strand = TRUE)/119349)/2 * 1000000
TE$rec_piRNA_ONLY_FLANKS <- (countOverlaps(TE_1000_longer, piRNA_rec, ignore.strand = TRUE)/226506 + countOverlaps(TE_1000_longer, piRNA_rec_2, ignore.strand = TRUE)/286181)/2 *1000000 - TE$rec_piRNA_BODY
TE$dys_piRNA_ONLY_FLANKS <- (countOverlaps(TE_1000_longer, piRNA_dis, ignore.strand = TRUE)/133438 + countOverlaps(TE_1000_longer, piRNA_dis_2, ignore.strand = TRUE)/119349)/2 *1000000- TE$dis_piRNA_BODY

TE$parent_9_BODY <- (countOverlaps(TE, piRNA_9, ignore.strand = TRUE)/257327 + countOverlaps(TE, piRNA_9_2, ignore.strand = TRUE)/269824)/2 *1000000
TE$parent_160_BODY <- (countOverlaps(TE, piRNA_160, ignore.strand = TRUE)/201072 + countOverlaps(TE, piRNA_160_2, ignore.strand = TRUE)/151954)/2*1000000
TE$parent_9_FLANK <- (countOverlaps(TE_1000_longer, piRNA_9, ignore.strand = TRUE)/257327 + countOverlaps(TE_1000_longer, piRNA_9_2, ignore.strand = TRUE)/269824)/2 *1000000 - TE$parent_9_BODY
TE$parent_160_FLANK <- (countOverlaps(TE_1000_longer, piRNA_160, ignore.strand = TRUE)/201072 + countOverlaps(TE_1000_longer, piRNA_160_2, ignore.strand = TRUE)/151954)/2*1000000 - TE$parent_160_BODY

View(as.data.frame(TE))

table(TE$parent_9_FLANK == 0 & TE$parent_160_FLANK == 0)
TE$TE_length = ifelse(TE$percent > 0.9, 'Full', '<90%')


####################
ggplot(data = as.data.frame(TE), aes(x = parent_160_FLANK, y = parent_9_FLANK, col = TE_length)) +
  geom_point(size = 4) + theme_bw() + ggtitle('Each dot represent piRNA enrichment\nof 60 asymmetric insertions in 160 genome\n36 of them have 0 piRNA in 9 and 160 strain, 0-2h Embryos') +
  xlab('piRNA enrichment in 1Kb flanking sequences, piRNA from 160 Embryo, RPM') +
  ylab('piRNA enrichment in 1Kb flanking sequences\npiRNA from 9 Embryo, RPM')
library(reshape)
TE_df = as.data.frame(TE)
#table(TE_df$rec_piRNA_ONLY_FLANKS == 0 & TE_df$dys_piRNA_ONLY_FLANKS == 0)
#ggplot(data = as.data.frame(TE_df), aes(x = dys_piRNA_ONLY_FLANKS, y = rec_piRNA_ONLY_FLANKS, col = TE_length)) +
#  geom_point(size = 4) + theme_bw() + ggtitle('Each dot represent piRNA enrichment\nof 60 asymmetric insertions in 160 genome\n37 of them have 0 piRNA in Dys and Rec samples, F2 0-2h Embryos') +
#  xlab('piRNA enrichment in 1Kb flanking sequences, piRNA from Dysgenic Embryo, RPM') +
#  ylab('piRNA enrichment in 1Kb flanking sequences\npiRNA from Reciprocal Embryo, RPM') + facet_wrap(.~Element, scales = 'free')

#TE_df <- TE_df[,c('rec_piRNA_ONLY_FLANKS', 'dys_piRNA_ONLY_FLANKS')]
TE_df <- TE_df[,c('parent_9_FLANK', 'parent_160_FLANK')]
melted = melt(TE_df)
ggplot(data = melted, aes(x = variable, y =value)) + geom_boxplot() + 
  theme_bw() + ggtitle('Total amount of piRNA in 60 asymmetric ins\nin160 genome between 9 and 160 embryos F0') +
  ylab('RPM, piRNA') + xlab('F0 0-2h embryos')
melted
wilcox.test(melted$value ~ melted$variable) 


##############################################slide4
setwd('~/2021_year/GENES_piRNA_automatic/trimmed/Unmapped_fasta/align_piRNA_to_44_TE/')
base44_72 <- bam_analysis('SRR2096072_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096072', 1523998)
base44_79 <- bam_analysis('SRR2096079_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096079', 1502275)
base44_80 <- bam_analysis('SRR2096080_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096080', 2582650)
base44_81 <- bam_analysis('SRR2096081_trimmed.fq.smallRNA_filtered.fasta.piRNA.fasta.bam', 'SRR2096081', 3215711)


base44_72 <- base44_72[order(base44_72$TE),]
base44_79 <- base44_79[order(base44_79$TE),]
base44_80 <- base44_80[order(base44_80$TE),]
base44_81 <- base44_81[order(base44_81$TE),]

F2_piRNA = data.frame(TE = base44_72$TE, Dysgenic = (base44_72$piRNA + base44_79$piRNA)/2, Reciprocal = (base44_80$piRNA + base44_81$piRNA)/2)
F2_piRNA$TE_label <- NA
F2_piRNA$TE_label <- ifelse(as.character(F2_piRNA$TE) %in% c('Polyphemus_full#DNA/TcMar-Tc1', 'Paris_full#DNA/TcMar-Tc1', 'Helena-2-2014_full#LINE/Jockey', 'Penelope_full#LINE/Penelope'), as.character(F2_piRNA$TE), NA)

ggplot(data = F2_piRNA, aes(x = log2(Dysgenic), y = log2(Reciprocal))) + geom_point(size = 2) + theme_bw() +
  geom_text(aes(label=unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),hjust=0, vjust=0) +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('Total amount of piRNA between F2 Embryos 0-2h, RPM') +
  xlab('log2(piRNA from Dysgenic Embryo, RPM)') +
  ylab('log2(piRNA from Reciprocal Embryo, RPM)') + theme(text = element_text(size = 14))


####################################################slide 7
setwd('~/2021_year/Kelleher/F2_RNA_seq_ref/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11, 12, 13)]
colData <- data.frame("replicate" = c("R1","R2","R3","R4", "R1","R2","R3","R4"), "treatment" = c("dys","dys", "dys","dys", "rec","rec", "rec","rec"))
attributes(colData)$row.names <- c("R1.dys","R2.dys","R3.dys","R4.dys", "R1.rec","R2.rec", "R3.rec","R4.rec")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "rec")
dds$replicate <- relevel(dds$replicate, ref = "R1")

keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","dys","rec"))
res2 <- as.data.frame(res2)

expressed_and_passed <- res2[res2$padj < 0.05,]
annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'gene']
library(GenomicRanges)
library(RColorBrewer)

Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
Gene <- Gene[!is.na(Nearest_ins),]
Nearest_ins <- nearest(Gene, TE, ignore.strand = T)
nrow(as.data.frame(Gene))

GENES_INSIDE <- subsetByOverlaps(Gene, TE, ignore.strand = TRUE)
GENES_INSIDE <- GENES_INSIDE$gene_id


Gene$closest_TE <- TE[Nearest_ins]$Element
Gene$closest_TE_per <- TE[Nearest_ins]$percent
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

Gene$IS_DE <- ifelse(Gene$gene_id %in% rownames(expressed_and_passed), 'DE', 'not DE')
Gene$TE_start <- start(TE[Nearest_ins])
Gene$TE_end <- end(TE[Nearest_ins])
Gene <- Gene[as.numeric(as.character(Gene$closest_distance)) < 5000,]
Gene
table(Gene$gene_id %in% rownames(expressed_and_passed))
table(Gene$gene_id %in% rownames(res2))
Gene <- Gene[Gene$gene_id %in% rownames(res2),]

Gene$log2FC <- NA
for (i in 1:length(Gene)){
  print(i)
  Gene$log2FC[i] <- ifelse(Gene$gene_id[i] %in% rownames(res2),
                           res2[which(rownames(res2) %in% Gene$gene_id[i]),]$log2FoldChange,'Nothing')
}

Gene$log2FC
table(Gene$gene_id %in% rownames(expressed_and_passed))

ggplot(data = as.data.frame(Gene), aes(x = as.numeric(as.character(closest_distance)), y = as.numeric(as.character(log2FC)), col = IS_DE)) +
  geom_point(size = 3) + theme_bw()+
  ggtitle('Distance to the closest TE (60 ins)\nfor all genes with limits 5000') + 
  labs(col = 'Is a gene DE?') +
  xlab('Distance from gene to TE, bp') + ylab('log2FC, ref = reciprocal') +
  theme(text = element_text(size = 14))
  
Gene$group <- NA
Gene$closest_distance <- as.numeric(Gene$closest_distance)
Gene$group <- ifelse(Gene$closest_distance == 0, '0',
                     ifelse(Gene$closest_distance < 1000, '0-1Kb',
                            ifelse(Gene$closest_distance < 2000, '1-2Kb',
                                   ifelse(Gene$closest_distance < 3000, '2-3Kb',
                                          ifelse(Gene$closest_distance < 4000, '3-4Kb',
                                                 '4-5Kb')))))

sum(table(Gene$group))
table(Gene$group)
ggplot(data = as.data.frame(Gene[Gene$IS_DE == 'DE',]), aes(x = group, y = log2FC, fill = group)) + geom_boxplot() + geom_jitter() + 
  theme_bw() +
  ggtitle('Each dot represent a gene, \ngroups show the distance to the closest TE\nlimits: -1.5, +1.5 log2FC (ref = reciprocal)')+
  ylim(c(-1.5, 1.5))



table(Gene$IS_DE)




Permtutation_test <- function(d, P){
  levels(d$label)
  print('Obserbed test-statistics:')
  Test_stat1 <- abs(abs(mean(d[d$label == levels(d$label)[1],]$log2FC)) - abs(mean(d[d$label == levels(d$label)[2],]$log2FC)))
  print(Test_stat1)
  Test_stat2 <- abs(abs(median(d[d$label == levels(d$label)[1],]$log2FC)) - abs(median(d[d$label == levels(d$label)[2],]$log2FC)))
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
    Perm_test_stat1[i] <- abs(abs(mean(Perm_matrix[d$label ==levels(d$label)[1], i])) - abs(mean(Perm_matrix[d$label == levels(d$label)[2], i]) ))
    Perm_test_stat2[i] <- abs(abs(median(Perm_matrix[d$label == levels(d$label)[1], i])) - abs(median(Perm_matrix[d$label == levels(d$label)[2], i]) ))
  }
  
  stat1 = as.data.frame(table(Perm_test_stat1 > Test_stat1))
  stat2 = as.data.frame(table(Perm_test_stat2 > Test_stat2))
  
  print('p-value1:')
  print(stat1[stat1$Var == 'TRUE',]$Freq/P)
  print('p-value2:')
  print(stat2[stat2$Var == 'TRUE',]$Freq/P)
}

res2$log2FoldChange
nrow(res2)

table(Gene$group)

d <- data.frame(gene_id = Gene[Gene$group == '4-5Kb',]$gene_id, log2FC = Gene[Gene$group == '4-5Kb',]$log2FC, label = 'INSIDE')
d <- unique(d)
others = res2[!rownames(res2) %in% d$gene_id,]$log2FoldChange
ot = data.frame(log2FC = others, label = 'NOT INSIDE')
d$gene_id <- NULL
d = rbind(d, ot)
d$label <- as.factor(d$label)

Permtutation_test(d, 100000)


View(as.data.frame(Gene))

