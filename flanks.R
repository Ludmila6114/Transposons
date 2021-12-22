library(Biostrings)
library(ORFik)

##check for each genome, use only repeat masker
base = readDNAStringSet('~/TE_bases/44_TE/dvir_full-size_TEs.fasta')
TE = data.frame(TE = unlist(lapply(strsplit(base@ranges@NAMES, '#'), function(x)x[1])), width = base@ranges@width)

#160
Neutral = read.table('~/RepeatMasker_UPD_NEWBASE_FINAL/RepeatMasker_160/160_chromosome.fasta.out', fill = T)
Neutral = Neutral[c(3:22845),]
Neutral = na.omit(Neutral)
Neutral = read.table('~/RepeatMasker_UPD_NEWBASE_FINAL/RepeatMasker_9/9_correct_names.fasta.out', fill = T)
#Neutral = Neutral[c(3:25909),]
#Neutral = na.omit(Neutral)

#Neutral = read.table('~/check_china_homology/chinaN/BEST_chinaN_assembly.fasta.out', fill = T)
#Neutral = Neutral[c(3:25088),]
#Neutral = na.omit(Neutral)

#Neutral = read.table('~/check_china_homology/chinaM/assembly.fasta.out', fill = T)
#Neutral = Neutral[c(3:27806),]
#Neutral = na.omit(Neutral)

#genome = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/1051_china_M_assembly/assembly.fasta')
#genome = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/Genome_assembly_chinaN/assembly_all_reads_without_short_mt/BEST_chinaN_assembly.fasta')

genome = readDNAStringSet('~/160_chromosome_assembly/160_chromosome.fasta')
genome@ranges@NAMES = unlist(lapply(strsplit(genome@ranges@NAMES, ' '), function(x)x[1]))

#genome = readDNAStringSet('~/9_POLISHED/9_BEST_NAMES.fasta')
#genome@ranges@NAMES = paste(genome@ranges@NAMES,'_', sep = '')

Neutral$TE_width = TE[match(as.character(Neutral$V10), TE$TE),]$width
Neutral$range = as.numeric(as.character(Neutral$V7)) - as.numeric(as.character(Neutral$V6))
Neutral$percent = Neutral$range/Neutral$TE_width * 100
Neutral <- Neutral[Neutral$percent >= 90,]
TE_table = data.frame(table(as.character(Neutral$V10)))
TE_table$class = Neutral[match(TE_table$Var1, as.character(Neutral$V10)),]$V11

ggplot(data = TE_table, aes(x = reorder(Var1, Freq), y = Freq, fill = class)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  theme_bw() + ggtitle('Number of copies at least 90%, RepeatMasker\ngenome: 160 strain') +
  xlab('TE copy number') +
  ylab('TE name') +
  labs(fill = 'TE class') + 
  theme(text = element_text(size = 15))

Neutral_GR = GRanges(seqnames = as.character(Neutral$V5), ranges = IRanges(start = as.numeric(as.character(Neutral$V6)),
                                                                           end = as.numeric(as.character(Neutral$V7))),
                     strand = ifelse(as.character(Neutral$V9) == '+', '+', '-'),
                     Element = as.character(Neutral$V10),
                     class = as.character(Neutral$V11),
                     percent = Neutral$percent, sequence = NA)

Neutral_GR$sequence = getSeq(genome, Neutral_GR)
View(as.data.frame(Neutral_GR))

for (i in 1:length(Neutral_GR)){
  Neutral_GR$PLUS_MAX[i] = max(unlist(findORFs(as.character(Neutral_GR$sequence[i]), startCodon = "ATG", minimumLength = 0))@width)
  Neutral_GR$MINUS_MAX[i] = max(unlist(findORFs(as.character(reverseComplement(Neutral_GR$sequence[i])), startCodon = "ATG", minimumLength = 0))@width)
  if ((i %% 100) == 0){
    print(i)
  }
}

length(Neutral_GR[Neutral_GR$Element == 'Nausicaa_full' & Neutral_GR$PLUS_MAX == 4821])/length(Neutral_GR[Neutral_GR$Element == 'Nausicaa_full',])

Neutral_GR$Number = 1:length(Neutral_GR) 
TE_table$best_ORF_number = NA
TE_table$ORF_length = NA
TE_table$base_ORF_length = NA
TE_table$strand = NA

base@ranges@NAMES = unlist(lapply(strsplit(base@ranges@NAMES, '#'), function(x)x[1]))
for (i in 1:nrow(TE_table)){
  TE_table$base_ORF_length[i] = max(unlist(findORFs(as.character(base[base@ranges@NAMES == TE_table$Var1[i],]), startCodon = "ATG", minimumLength = 0))@width)
}

for (i in 1:nrow(TE_table)){
  TE_table$ORF_length[i] = max(Neutral_GR[Neutral_GR$Element == TE_table$Var1[i],]$PLUS_MAX)
  TE_table$strand[i] = as.character(strand(Neutral_GR[(Neutral_GR$PLUS_MAX == max(Neutral_GR[Neutral_GR$Element == TE_table$Var1[i],]$PLUS_MAX)) &
                                           Neutral_GR$Element  == TE_table$Var1[i],]))[1]
  
  
  TE_table$best_ORF_number[i] = Neutral_GR[(Neutral_GR$PLUS_MAX == max(Neutral_GR[Neutral_GR$Element == TE_table$Var1[i],]$PLUS_MAX)) &
                                             Neutral_GR$Element  == TE_table$Var1[i],]$Number[1]
}


ggplot(data = as.data.frame(Neutral_GR[Neutral_GR$Element == 'dvir-786_full',]), aes(x = PLUS_MAX)) +
  geom_histogram(fill = 'lightblue3') + 
  theme_bw() +
  xlab('ORF length, bp') + 
  ylab('Count') +
  ggtitle('dvir-786_full in 160')




TE_table = na.omit(TE_table)
ggplot(data = TE_table, aes(x = ORF_length, y = ORF_length_2, col = class)) + geom_point() + 
  xlab('Max ORF length in 160') +
  ylab('Max ORF length in ChinaM') + theme_bw() +  geom_label_repel(aes(label = TE_table$Var1),
                                                               box.padding   = 0.2, 
                                                               point.padding = 0.5,
                                                               segment.color = 'grey50') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', col = 'red')



ggplot(data = Neutral, aes(x = V10, y = percent, fill = V11)) + geom_boxplot() + 
  facet_wrap(~V11, scales = 'free') + coord_flip() + theme_bw() +
  geom_jitter() +
  ggtitle('ChinaN full-length TEs, at least 90% of canon') +
  xlab('TE name & class') +
  labs(fill = 'TE class') +
  ylab('TE in genome, percent from canonical copy')


######################################





setwd('~/check_china_homology/')

base = readDNAStringSet('~/TE_bases/44_TE/dvir_full-size_TEs.fasta')
TE = data.frame(TE = unlist(lapply(strsplit(base@ranges@NAMES, '#'), function(x)x[1])), width = base@ranges@width)

#genome <- readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/ASSEMBLY_101M/assembly.fasta')
#genome <- readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/101N_genome_assembly/assembly.fasta')
genome = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/Genome_assembly_chinaN/assembly_all_reads_without_short_mt/BEST_chinaN_assembly.fasta')
#genome = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/1051_china_M_assembly/assembly.fasta')
contigs = data.frame(contig = genome@ranges@NAMES, width = genome@ranges@width)

#Neutral = read.table('M_101/assembly.fasta.out', fill = T)
#Neutral = read.table('101N/assembly.fasta.out', fill = T)
#Neutral = read.table('chinaN/BEST_chinaN_assembly.fasta.out', fill = T)
#Neutral = read.table('chinaM/assembly.fasta.out', fill = T)
#Neutral = Neutral[Neutral$V10 == 'DROGYPSYZ_full',]
Neutral = read.table('~/RepeatMasker_UPD_NEWBASE_FINAL/RepeatMasker_160/160_chromosome.fasta.out', fill = T)

Neutral = Neutral[c(3:22845),]
Neutral = na.omit(Neutral)

Neutral$TE_width = TE[match(as.character(Neutral$V10), TE$TE),]$width
Neutral$contig_length = contigs[match(as.character(Neutral$V5), contigs$contig),]$width
Neutral$range = as.numeric(as.character(Neutral$V7)) - as.numeric(as.character(Neutral$V6))
Neutral$percent = Neutral$range/Neutral$TE_width * 100
Neutral <- Neutral[Neutral$percent >= 90,]
Neutral = Neutral[Neutral$contig_length - as.numeric(as.character(Neutral$V7)) >= 3000 & as.numeric(as.character(Neutral$V6)) >= 3000,]


TE_table = data.frame(table(as.character(Neutral$V10)))
TE_table$class = Neutral[match(TE_table$Var1, as.character(Neutral$V10)),]$V11


ggplot(data = Neutral, aes(x = V10, y = percent, fill = V11)) + geom_boxplot() + 
  facet_wrap(~V11, scales = 'free') + coord_flip() + theme_bw() +
  geom_jitter() +
  ggtitle('ChinaN full-length TEs, at least 90% of canon') +
  xlab('TE name & class') +
  labs(fill = 'TE class') +
  ylab('TE in genome, percent from canonical copy')



ggplot(data = TE_table, aes(x = V1, fill = V11)) + geom_bar() + coord_flip() + theme_bw()


#ggplot(data = Neutral, aes(x = reorder(V10,percent), y = percent,fill = V11)) + 
#  geom_boxplot() + 
#  coord_flip() + 
#  theme_bw() +
#  geom_jitter() +
#  ggtitle('ChinaN full-length TEs, at least 90%')


ggplot(data = Neutral, aes(x = forcats::fct_infreq(unlist(lapply(strsplit(as.character(V10), '_'), function(x)x[1]))), fill = V11)) + geom_bar() + 
  coord_flip() + theme_bw() + scale_x_discrete() + facet_wrap(~V11, scales = 'free') + 
  ggtitle('ChinaN TEs with at least 90% length') + ylab('Count') 
+ xlab('TE name & class') + labs(fill = 'TE class')



Neutral_GR = GRanges(seqnames = as.character(Neutral$V5), ranges = IRanges(start = as.numeric(as.character(Neutral$V6)),
                                                                           end = as.numeric(as.character(Neutral$V7))),
                     strand = ifelse(as.character(Neutral$V9) == '+', '+', '-'),
                     Element = as.character(Neutral$V10),
                     class = as.character(Neutral$V11),
                     percent = Neutral$percent, status = NA)



Neutral_GR$flanc_left <- IRanges( start = start(ranges(Neutral_GR))-3000, end=start(ranges(Neutral_GR)))
Neutral_GR$flanc_right <- IRanges( start = end(ranges(Neutral_GR)), end=end(ranges(Neutral_GR))+3000)
TE_left <- GRanges(seqnames = seqnames(Neutral_GR), strand = strand(Neutral_GR), ranges=Neutral_GR$flanc_left)
Neutral_GR$left_fl_seq <- getSeq(genome, TE_left)

TE_right <- GRanges(seqnames = seqnames(Neutral_GR), strand = strand(Neutral_GR), ranges=Neutral_GR$flanc_right)
Neutral_GR$rigth_fl_seq <- getSeq(genome, TE_right)

#View(as.data.frame(Neutral_GR))
Neutral_GR$dist <- NA
count = 0
count_0 = vector()
count_4 = 0

#genome_101M = readDNAStringSet('/mnt/raid/pro_milka/GENOME_ASSEMBLIES_ONT/ASSEMBLY_101M/assembly.fasta')
genome_101M = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/1051_china_M_assembly/assembly.fasta')
#genome_101M = readDNAStringSet('/mnt/raid/pro_milka/GENOME_ASSEMBLIES_ONT/101N_genome_assembly/assembly.fasta')
#genome_101M = readDNAStringSet('~/GENOME_ASSEMBLIES_ONT/Genome_assembly_chinaN/assembly_all_reads_without_short_mt/BEST_chinaN_assembly.fasta')
#genome_101M = readDNAStringSet('~/9_POLISHED/9_BEST_NAMES.fasta')
#genome_101M = readDNAStringSet('~/160_chromosome_assembly/160_chromosome.fasta')
#names(genome_101M) <- unlist(lapply(strsplit(names(genome_101M), ' '), function(x)x[1]))


for (N in 1:length(Neutral_GR)){
  left_flanc <- Neutral_GR$left_fl_seq[N]
  names(left_flanc) <- paste(N, ' left_flanc_seq')
  writeXStringSet(left_flanc, '~/left_flanc.fasta')
  left_blast <- system('blastn -task megablast -query ~/left_flanc.fasta -db ~/GENOME_ASSEMBLIES_ONT/1051_china_M_assembly/assembly.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  left_df <- as.data.frame(do.call('rbind', strsplit(left_blast, '\t')))
  
  right_flanc <- Neutral_GR$rigth_fl_seq[N]
  names(right_flanc) <- paste(N, ' rigth_flanc_seq')
  writeXStringSet(right_flanc, '~/right_flanc.fasta')
  right_blast <- system('blastn -task megablast -query ~/right_flanc.fasta -db ~/GENOME_ASSEMBLIES_ONT/1051_china_M_assembly/assembly.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
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
        Neutral_GR$dist[N] <- 0
      }
      else{
        Neutral_GR$dist[N] = min(abs(start(right) - end(left)), abs(start(left) - end(right)))
        new_blast = GRanges(seqnames = left_df$V2, ranges = IRanges(start = min(end(right), end(left)), end = max(start(left), start(right))))
        strange = getSeq(genome_101M, new_blast)
        names(strange) <- paste(N, 'strange')
        writeXStringSet(strange, '~/strange.fasta')
        blast_strange <- system('blastn -task megablast -query ~/strange.fasta -db ~/TE_bases/44_TE/dvir_full-size_TEs.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
        blast_df <- as.data.frame(do.call('rbind', strsplit(blast_strange, '\t')))
        if (nrow(blast_df) == 0){
          print(paste('nothing in 44 TE with length:', width(strange), 'number', N))
          count_0 = c(count_0, width(strange))
          if (width(strange) > 1000){
            setwd('~/')
            print(paste('number:', N, 'first genome TE'))
            print('In the term of second genome:')
            print(new_blast)
            #writeXStringSet(strange, paste(width(strange),'_TE_STRANGE_',N,'.fasta'))
          }
          else{
            print(paste('smth unmapped with length:', width(strange), 'number', N))
          }
          
          
        }
        else{
          count_4 = count_4 + 1
          blast_df$TE <- Neutral_GR[N]$Element
          blast_df$V2 = unlist(lapply(strsplit(as.character(blast_df$V2), '#'), function(x)x[1]))
          if ((Neutral_GR$dist[N]/as.numeric(as.character(Neutral_GR[N]@ranges@width))*100) > 105){
            print(paste('smth with my TE, number', N, 'Dist:', Neutral_GR$dist[N], 'TE canon length:', Neutral_GR[N]@ranges@width))
          } 
          
          if (!(Neutral_GR[N]$Element %in% blast_df$V2)){
            print(paste(blast_df$V2, ' IN chinaN it was:', Neutral_GR[N]$Element, 'Initial number:', N))
          }
        }
      }
    }
  }
}

Neutral[159,]
Neutral[654,]
table(!is.na(Neutral_GR$dist))
Neutral_GR = Neutral_GR[!is.na(Neutral_GR$dist),]
length(Neutral_GR[Neutral_GR$dist < 1000,])
Neutral_GR$p = Neutral_GR$dist/Neutral_GR@ranges@width
length(Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,])


save_160 = Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,]
save_9 = Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,]
save_101N = Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,]
save_chinaM = Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,]

save_160$ID = paste(seqnames(save_160), ranges(save_160))
save_9$ID = paste(seqnames(save_9), ranges(save_9))
save_101N$ID = paste(seqnames(save_101N), ranges(save_101N))
save_chinaM$ID = paste(seqnames(save_chinaM), ranges(save_chinaM))

View(as.data.frame(save_9[save_9$ID %in% intersect(intersect(intersect(save_101N$ID, save_160$ID), save_9$ID), save_chinaM$ID),]))

intersect(save_101N$ID, save_160$ID)
##########Distance and class
ggplot(data = as.data.frame(Neutral_GR[Neutral_GR$dist < 10000,]), aes(x = p, y = reorder(Element, p), fill = class)) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  ylab('Element') +
  xlab('Percent of canon between flanks in 9 genome')

View(as.data.frame(Neutral_GR))
########################3



ggplot(data = as.data.frame(Neutral_GR[Neutral_GR$p >= 0.9 & Neutral_GR$p < 2,]), aes(x = unlist(lapply(strsplit(as.character(Element), '_'), function(x)x[1])), fill = class)) +
  geom_bar() + 
  facet_wrap(~class, scales = 'free') + 
  coord_flip() + 
  theme_bw() +
  ggtitle('Транспозоны, присутствуюшие в chinaN (полноразмерно),\nи сохранившие локализацию в аналогичных фланках в 160') +
  labs(fill = 'TE class') +
  xlab('TE name & class') + ylab('Count')


View(as.data.frame(Neutral_GR))
table(Neutral$V10, Neutral$V11)

Neutral_GR$p = Neutral_GR$dist/Neutral_GR@ranges@width
length(Neutral_GR[Neutral_GR$dist < 20,])
min(Neutral_GR[Neutral_GR$dist < 950,]$p)

length(Neutral_GR[Neutral_GR$Element == 'NewTE4_full'
                  & Neutral_GR$dist > 20,])
