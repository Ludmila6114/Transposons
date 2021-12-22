#Подготовка
setwd('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11)]
colData <- data.frame("replicate" = c("R1","R2", "R3", "R1","R2", "R3"), 
                      "treatment" = c("9","9","9", "160", "160", "160"))
attributes(colData)$row.names <- c("9_YOUNG_R1","9_YOUNG_R2", "9_YOUNG_R3", 
                                   "160_YOUNG_R1","160_YOUNG_R2", "160_YOUNG_R3")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "9")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","160","9"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  

normalized_counts <- counts(dds50, normalized=TRUE)
normalized_counts$name = rownames(normalized_counts)
featurecounts$name = rownames(featurecounts)
'FBgn0205574' %in% u 


write.csv(as.data.frame(res2), 'expressed_F0_YOUNG.csv', quote = F)
write.csv(as.data.frame(expressed_and_passed), 'DE_F0_YOUNG.csv', quote = F)

GO <- expressed_and_passed
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(rownames(GO), ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym

#Down:
GO.ENSEMBL <- bitr(GO[GO$log2FoldChange >0,]$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)
ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ego, showCategory=30)
length(GO.ENSEMBL$ENSEMBL)
goplot(ego2)


expressed_and_passed$label <- ort[match(rownames(expressed_and_passed), ort$V6),]$V1
expressed_and_passed$label <- as.character(mapIds(org.Dm.eg.db, as.character(expressed_and_passed$label), 'SYMBOL','FLYBASE'))
expressed_and_passed$label = ifelse(expressed_and_passed$label == 'NULL', rownames(expressed_and_passed), expressed_and_passed$label)
expressed_and_passed$label = ifelse(expressed_and_passed$pvalue < 0.001 & abs(expressed_and_passed$log2FoldChange) > 3, expressed_and_passed$label, NA)
expressed_and_passed$label[!is.na(expressed_and_passed$label)]

ggplot(data = expressed_and_passed, aes(x = log2FoldChange, y = -log10(pvalue), col = log10(baseMean))) +
  theme_bw() + geom_point() + 
  geom_label_repel(aes(label = expressed_and_passed$label),
                   box.padding   = 0.2, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  xlab('log2FC, ref = 9')



#d <- as.character(mapIds(org.Dm.eg.db, as.character(GO[GO$log2FoldChange > 0,]$dmel_ort), 'UNIPROT', 'FLYBASE'))
#k <- as.data.frame(enrichKEGG(gene = d,organism = 'dme',keyType = 'uniprot'))

up = bitr(GO[GO$log2FoldChange >0,]$dmel.sym, fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = org.Dm.eg.db)

up = up$ENTREZID
length(up)

down = bitr(GO[GO$log2FoldChange <0,]$dmel.sym, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Dm.eg.db)

down = down$ENTREZID
length(down)

up_kegga = limma::topKEGG(limma::kegga(de = up, 
                                 species = 'Dm'))

down_kegga = limma::topKEGG(limma::kegga(de = down, 
                                      species = 'Dm'))



#Сравнение закладки транскриптов транспозонов у молодых родителей:

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
  bam_df$type <- 'mRNA'
  
  result <- data.frame(TE = unique(bam_df$rname))
  
  arr <- function(my_data, name){
    result$mRNA <- NA
    for (i in 1:nrow(result)){
      if (result$TE[i] %in% my_data$Var1){
        result$mRNA[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'mRNA')),]$Freq
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
  result$mRNA <- result$mRNA * 1000000/RPM
  #result$siRNA <- result$siRNA * 1000000/RPM
  #result$Ultra_long <- result$Ultra_long * 1000000/RPM
  #result$Ultra_short <- result$Ultra_short * 1000000/RPM
  print('ALL')
  return(result)
}

young_9_15 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096015_trimmed.fq.bam', '9_sample15', 29331)
young_9_16 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096016_trimmed.fq.bam', '9_sample16', 28285)
young_9_17 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096017_trimmed.fq.bam', '9_sample17', 23133)

TEs <- readDNAStringSet('~/TE_bases/44_TE/dvir_full-size_TEs.fasta')
vector = names(TEs)

vector[!vector %in% young_9_17$TE]
young_9_15[43,] <- c("dvir-1002_full#RC/Helitron", 0, '9_sample15')
young_9_15[44,] <- c("dvir-552_full#LINE/Jockey", 0, '9_sample15')

young_9_16[40,] <- c("dvir-1002_full#RC/Helitron", 0, '9_sample16')
young_9_16[41,] <- c("dvir-552_full#LINE/Jockey", 0, '9_sample16')
young_9_16[42,] <- c("dvir-786_full#RC/Helitron", 0, '9_sample16')
young_9_16[43,] <- c("Paris_full#DNA/TcMar-Tc1", 0, '9_sample16')
young_9_16[44,] <- c("Polyphemus_full#DNA/TcMar-Tc1", 0, '9_sample16')

young_9_17[40,] <- c("dvir-1002_full#RC/Helitron", 0, '9_sample17')
young_9_17[41,] <- c("dvir-390_full#LINE/CR1", 0, '9_sample17')
young_9_17[42,] <- c("dvir-593_full#DNA/TcMar-Tc1", 0, '9_sample17')
young_9_17[43,] <- c("Penelope_full#LINE/Penelope", 0, '9_sample17')
young_9_17[44,] <- c("Polyphemus_full#DNA/TcMar-Tc1", 0, '9_sample17')


young_9_15 <- young_9_15[order(young_9_15$TE),]
young_9_16 <- young_9_16[order(young_9_16$TE),]
young_9_17 <- young_9_17[order(young_9_17$TE),]

young_9 = young_9_15
young_9$strain <- NULL
young_9$mRNA_young_9 = (as.numeric(as.character(young_9_15$mRNA)) + 
                          as.numeric(as.character(young_9_16$mRNA)) +
                          as.numeric(as.character(young_9_17$mRNA))) /3 
young_9$mRNA <- NULL


young_160_18 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096018_trimmed.fq.bam', '160_sample18', 15079)
young_160_19 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096019_trimmed.fq.bam', '160_sample19', 14824)
young_160_20 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/SRR2096020_trimmed.fq.bam', '160_sample20', 12890)

TEs <- readDNAStringSet('~/TE_bases/44_TE/dvir_full-size_TEs.fasta')
vector = names(TEs)
vector[!vector %in% young_160_20$TE]

young_160_18[41,] <- c("dvir-1002_full#RC/Helitron", 0, '160_sample18')
young_160_18[42,] <- c("dvir-390_full#LINE/CR1", 0, '160_sample18')
young_160_18[43,] <- c("dvir-552_full#LINE/Jockey", 0, '160_sample18')
young_160_18[44,] <- c( "dvir-593_full#DNA/TcMar-Tc1", 0, '160_sample18')

young_160_19[41,] <- c("dvir-1002_full#RC/Helitron", 0, '160_sample19')
young_160_19[42,] <- c("dvir-390_full#LINE/CR1", 0, '160_sample19')
young_160_19[43,] <- c("dvir-552_full#LINE/Jockey", 0, '160_sample19')
young_160_19[44,] <- c("NewTE7_full#DNA/TcMar-Tc1", 0, '160_sample19')

young_160_20[42,] <- c("dvir-390_full#LINE/CR1", 0, '160_sample20')
young_160_20[43,] <- c("dvir-552_full#LINE/Jockey", 0, '160_sample20')
young_160_20[44,] <- c( "dvir-593_full#DNA/TcMar-Tc1", 0, '160_sample20')

young_160_18 <- young_160_18[order(young_160_18$TE),]
young_160_19 <- young_160_19[order(young_160_19$TE),]
young_160_20 <- young_160_20[order(young_160_20$TE),]

young_160 = young_160_18
young_160$strain <- NULL
young_160$mRNA_young_160 = (as.numeric(as.character(young_160_18$mRNA)) + 
                          as.numeric(as.character(young_160_19$mRNA)) +
                          as.numeric(as.character(young_160_20$mRNA))) /3 
young_160$mRNA <- NULL

young_df = young_160
young_df$mRNA_young_9 = young_9$mRNA_young_9
write.csv(young_df, '~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/TEs_young_RPM.csv')
young_df = read.csv('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_44/TEs_young_RPM.csv')

young_df$TE_length = TEs@ranges@width
names(TEs) == young_df$TE

young_df = young_df[young_df$mRNA_young_160 != 0 & young_df$mRNA_young_9 !=0,]
library(ggrepel)

ggplot(data = young_df, aes(x = log10(mRNA_young_9/TE_length), y = log10(mRNA_young_160/TE_length))) + 
  geom_point(color = ifelse(log10(young_df$mRNA_young_160)/log10(young_df$mRNA_young_9) > 2 | log10(young_df$mRNA_young_160)/log10(young_df$mRNA_young_9) < 0.5, 'red', 'black')) +
  xlab('Young 9, log10 TPM') + ylab('Young 160, log10 TPM')+
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  #geom_text_repel(aes(label=unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),hjust=0, vjust=0) +
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                  box.padding   = 0.2, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() + theme(text = element_text(size = 20))  + ggtitle('TEs expression level in F0 embryos\nfrom young (7-15 days) parents')   


#Есть ли корреляция между числом ольших копий транспозлона и его экскпрессией в эухроматине

P = read.csv('~/2021_year/TE_in_genomes_analysis/160_TE_90%_counts.csv')
M = read.csv('~/2021_year/TE_in_genomes_analysis/9_TE_90%_counts.csv')


cor_check = data.frame(TE = rep(unlist(lapply(strsplit(as.character(young_df$TE), '_'), function(x)x[1])), 2),
                       expression_9 = rep(young_df$mRNA_young_9/young_df$TE_length, 2),
                       expression_160 = rep(young_df$mRNA_young_160/young_df$TE_length, 2),
                       chromatin = c(rep('Euchromatin', 44), 
                                    rep('Heterochromatin', 44)),
                       copy_9 = NA, copy_160 = NA)
TEs$`Nausicaa_full#LTR/Gypsy`

for (i in 1:nrow(P)){
  if (P$Var2[i] == 'Euchromatin'){
    for (j in 1:44){
      if (cor_check$TE[j] == unlist(lapply(strsplit(as.character(P$Var1[i]), '_'), function(x)x[1]))){
        cor_check$copy_160[j] = P$Freq[i]
      }
    }
  }
  else if (P$Var2[i] == 'Heterochromatin'){
    for (j in 45:88){
      if (cor_check$TE[j] == unlist(lapply(strsplit(as.character(P$Var1[i]), '_'), function(x)x[1]))){
        cor_check$copy_160[j] = P$Freq[i]
      }
    }
  }
  print(i)
}

for (i in 1:nrow(M)){
  if (M$Var2[i] == 'Euchromatin'){
    for (j in 1:44){
      if (cor_check$TE[j] == unlist(lapply(strsplit(as.character(M$Var1[i]), '_'), function(x)x[1]))){
        cor_check$copy_9[j] = M$Freq[i]
      }
    }
  }
  else if (M$Var2[i] == 'Heterochromatin'){
    for (j in 45:88){
      if (cor_check$TE[j] == unlist(lapply(strsplit(as.character(M$Var1[i]), '_'), function(x)x[1]))){
        cor_check$copy_9[j] = M$Freq[i]
      }
    }
  }
  print(i)
}



ggplot(data = cor_check[cor_check$chromatin == 'Euchromatin',], aes(x = copy_9, y = expression_9)) + 
  geom_point() + 
  xlab('Number of at least 90%-length copies in 9 strain, euchr') +
  ylab('TEs expression in TPM') + 
  theme_bw() + 
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                                box.padding   = 0.2, 
                                point.padding = 0.5,
                                segment.color = 'grey50') +
  theme(text = element_text(size = 20))  + 
  ggtitle('TEs expression level and copy numbers\nin 9 embryos\nfrom young (7-15 days) parents\nEUCHROMATIN ') 
  

cor.test(cor_check[cor_check$chromatin == 'Euchromatin',]$expression_160, 
     cor_check[cor_check$chromatin == 'Euchromatin',]$copy_160)





#Сравнение пиРНК у родителей (возраст непонятен, поискать в статье)

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
    #result$Ultra_long <- NA
    result$Ultra_short <- NA
    
    
    for (i in 1:nrow(result)){
      if (result$TE[i] %in% my_data$Var1){
        result$piRNA[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'piRNA')),]$Freq
        result$siRNA[i] <-  my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'siRNA')),]$Freq 
        #result$Ultra_long[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'Ultra_long')),]$Freq
        result$Ultra_short[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'Ultra_short')),]$Freq
      }
    }
    return(result)
  }
  
  data <- as.data.frame(table(bam_df$rname, bam_df$type))
  result <- arr(data, sample_id)
  result$strain <- sample_id
  result$piRNA <- result$piRNA * 1000000/RPM
  result$siRNA <- result$siRNA * 1000000/RPM
  #result$Ultra_long <- result$Ultra_long * 1000000/RPM
  result$Ultra_short <- result$Ultra_short * 1000000/RPM
  print('ALL')
  return(result)
}

setwd('~/2021_year/GENES_piRNA_automatic/trimmed/Unmapped_fasta/filtered_but_all_length_smallRNA/')
base44_26 <- bam_analysis('SRR2096026_trimmed.fq.smallRNA_filtered.fasta.bam', 'SRR2096026', 3524966)
base44_38 <- bam_analysis('SRR2096038_trimmed.fq.smallRNA_filtered.fasta.bam', 'SRR2096038', 3730311)
base44_50 <- bam_analysis('SRR2096050_trimmed.fq.smallRNA_filtered.fasta.bam', 'SRR2096050', 2166446)
base44_61 <- bam_analysis('SRR2096061_trimmed.fq.smallRNA_filtered.fasta.bam', 'SRR2096061', 1735121)

base44_26 <- base44_26[order(base44_26$TE),]
base44_38 <- base44_38[order(base44_38$TE),]
base44_50 <- base44_50[order(base44_50$TE),]
base44_61 <- base44_61[order(base44_61$TE),]



F0_piRNA = data.frame(TE = base44_26$TE, 
                      P_160 = (base44_50$piRNA + base44_61$piRNA)/2, 
                      M_9 = (base44_26$piRNA + base44_38$piRNA)/2)

ggplot(data = F0_piRNA, aes(y = log2(P_160), x = log2(M_9))) + 
  geom_point(color = ifelse(log2(F0_piRNA$P_160)/log2(F0_piRNA$M_9) > 2 | log2(F0_piRNA$P_160)/log2(F0_piRNA$M_9) < 0.5, 'red', 'black')) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('Total amount of piRNA between F0 Embryos 0-2h, RPM') +
  ylab('log2(piRNA from 160 Embryo, RPM)') +
  xlab('log2(piRNA from 9 Embryo, RPM)') +  
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                                                              box.padding   = 0.3, 
                                                              point.padding = 0.5,
                                                              segment.color = 'grey50') +
  theme_classic() +
  theme(text = element_text(size = 20))  


F0_siRNA = data.frame(TE = base44_26$TE, 
                      P_160 = (base44_50$siRNA + base44_61$siRNA)/2, 
                      M_9 = (base44_26$siRNA + base44_38$siRNA)/2)

ggplot(data = F0_siRNA, aes(y = log2(P_160), x = log2(M_9))) + 
  geom_point(color = ifelse(log2(F0_siRNA$P_160)/log2(F0_siRNA$M_9) > 2 | log2(F0_siRNA$P_160)/log2(F0_siRNA$M_9) < 0.5, 'red', 'black')) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('Total amount of siRNA between F0 Embryos 0-2h, RPM') +
  ylab('log2(siRNA from 160 Embryo, RPM)') +
  xlab('log2(siRNA from 9 Embryo, RPM)') +  
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                   box.padding   = 0.3, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +
  theme(text = element_text(size = 20)) 


F0_smth = data.frame(TE = base44_26$TE, 
                      P_160 = (base44_50$Ultra_short + base44_61$Ultra_short)/2, 
                      M_9 = (base44_26$Ultra_short + base44_38$Ultra_short)/2)

ggplot(data = F0_smth, aes(y = log2(P_160), x = log2(M_9))) + 
  geom_point(color = ifelse(log2(F0_smth$P_160)/log2(F0_smth$M_9) > 2 | log2(F0_smth$P_160)/log2(F0_smth$M_9) < 0.5, 'red', 'black')) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('Total amount of small RNA < 21 between F0 Embryos 0-2h, RPM') +
  ylab('log2(small RNA < 21bp from 160 Embryo, RPM)') +
  xlab('log2(small RNA < 21bp from 9 Embryo, RPM)') +  
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                   box.padding   = 0.3, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +
  theme(text = element_text(size = 20)) 


#Кореляция уровня экспрессии транспозона и пиРНК к нему
table(young_df$TE == F0_piRNA$TE)
Summary = young_df
Summary$piRNA_9 = F0_piRNA$M_9
Summary$piRNA_160 = F0_piRNA$P_160
Summary$mRNA_young_9

Summary$siRNA_9 = F0_siRNA$M_9
Summary$siRNA_160 = F0_siRNA$P_160


ggplot(data = Summary, aes(y = log10(siRNA_9), x = log10(mRNA_young_9))) + 
  geom_point(color = ifelse(log10(Summary$siRNA_9)/log10(Summary$mRNA_young_9) > 2 | log10(Summary$siRNA_9)/log10(Summary$mRNA_young_9) < 0.5, 'red', 'black')) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('TE expression and siRNA expression in 9 strain') +
  xlab('log10(mRNA of TEs) from 9 Embryo, RPM') +
  ylab('log10(siRNA) from 9 Embryo, RPM') +  
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                   box.padding   = 0.3, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +
  theme(text = element_text(size = 20)) 


ggplot(data = Summary, aes(y = log10(siRNA_160), x = log10(piRNA_160))) + 
  geom_point(color = ifelse(log10(Summary$siRNA_160)/log10(Summary$piRNA_160) > 2 | log10(Summary$siRNA_160)/log10(Summary$piRNA_160) < 0.5, 'red', 'black')) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  ggtitle('piRNA expression and siRNA expression\nin 160 strain') +
  xlab('log10(piRNA) from 160 Embryo, RPM') +
  ylab('log10(siRNA) from 160 Embryo, RPM') +  
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                   box.padding   = 0.3, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() +
  theme(text = element_text(size = 20)) 




###########################СТАРЫЕ:


############## Сравнение генов старых мух друг с другом:
setwd('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11)]
colData <- data.frame("replicate" = c("R1","R2", "R3", "R1","R2", "R3"), "treatment" = c("9","9", "9", "160", "160", "160"))
attributes(colData)$row.names <- c("9_OLD_R1","9_OLD_R2", "9_OLD_R3", "160_OLD_R1","160_OLD_R2", "160_OLD_R3")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "9")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
resOLD <- results(dds50,contrast=c("treatment","160","9"))
resOLD <- as.data.frame(resOLD)
expressed_and_passedOLD <- resOLD[resOLD$padj < 0.05,]  

write.csv(as.data.frame(res2), 'expressed_F0_OLD.csv', quote = F)
write.csv(as.data.frame(expressed_and_passed), 'DE_F0_OLD.csv', quote = F)

GO <- expressed_and_passedOLD
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(rownames(GO), ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym
GO.ENSEMBL <- bitr(GO[GO$log2FoldChange <0,]$dmel.sym, fromType = "SYMBOL",
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
length(GO.ENSEMBL$ENSEMBL)




expressed_and_passed_SHARED = expressed_and_passed[rownames(expressed_and_passed) %in% rownames(expressed_and_passedOLD),]
expressed_and_passed_SHARED$OLD = NA
expressed_and_passed_SHARED$OLD = 
  expressed_and_passedOLD[match(rownames(expressed_and_passed_SHARED), rownames(expressed_and_passedOLD)),]$log2FoldChange

ggplot(data = expressed_and_passed_SHARED, aes(x = log2FoldChange, y = OLD)) + 
  geom_point() + 
  theme_bw() +
  ggtitle('Направление различий экспрессии\nу общих для молодых и старых\nдифф экспрессных генов (1850 шт)') +
  xlab('log2FC, Молодые') + ylab('log2FC, Старые') +
  theme(text = element_text(size = 20)) 



DE_in_old = expressed_and_passedOLD[!rownames(expressed_and_passedOLD) %in% rownames(expressed_and_passed),]
nrow(DE_in_old)
GO <- DE_in_old
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(rownames(GO), ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym
GO$dmel_ort
nrow(GO)

table(twenty %in% rownames(GO))

GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

twenty %in% GO.ENSEMBL$ENSEMBL

GO.ENSEMBL = GO.ENSEMBL[GO.ENSEMBL$ENSEMBL %in% twenty,]
GO.ENSEMBL

GO.Flybase = bitr(unique(GO.ENSEMBL$SYMBOL), fromType = "SYMBOL",
                  toType = c("FLYBASE"),
                  OrgDb = org.Dm.eg.db)
GO.Flybase




#Down:
GO.ENSEMBL <- bitr(GO$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

nrow(GO.ENSEMBL)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ego, showCategory=30)
length(GO.ENSEMBL$ENSEMBL)
goplot(ego2)


twenty = as.data.frame(ego)
twenty = strsplit(as.character(ego$geneID), '/')
twenty = unlist(twenty)
which(twenty %in% GO.ENSEMBL$ENSEMBL)

'FBgn0053892' %in% GO.ENSEMBL$ENSEMBL

twenty
FBGN = GO.ENSEMBL[GO.ENSEMBL$ENSEMBL %in% twenty,]
FBGN$SYMBOL

nrow(GO)

GO.ENSEMBL

#Уровень закладки транскриптов транспозонов у старых родителей:

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
  bam_df$type <- 'mRNA'
  
  result <- data.frame(TE = unique(bam_df$rname))
  
  arr <- function(my_data, name){
    result$mRNA <- NA
    for (i in 1:nrow(result)){
      if (result$TE[i] %in% my_data$Var1){
        result$mRNA[i] <- my_data[which((my_data$Var1 %in% result$TE[i]) & (my_data$Var2 == 'mRNA')),]$Freq
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
  result$mRNA <- result$mRNA * 1000000/RPM
  #result$siRNA <- result$siRNA * 1000000/RPM
  #result$Ultra_long <- result$Ultra_long * 1000000/RPM
  #result$Ultra_short <- result$Ultra_short * 1000000/RPM
  print('ALL')
  return(result)
}
old_9_21 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096021_trimmed.fq.bam', '9_sample21', 24769)
old_9_22 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096022_trimmed.fq.bam', '9_sample22', 24516)
old_9_23 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096023_trimmed.fq.bam', '9_sample23', 29600)

TEs <- readDNAStringSet('~/TE_bases/44_TE/dvir_full-size_TEs.fasta')
vector = names(TEs)
vector[!vector %in% old_9_23$TE]
old_9_21[41,] <- c("dvir-1002_full#RC/Helitron", 0, '9_sample21')
old_9_21[42,] <- c("dvir-552_full#LINE/Jockey", 0, '9_sample21')
old_9_21[43,] <- c("dvir-593_full#DNA/TcMar-Tc1", 0, '9_sample21')
old_9_21[44,] <- c("dvir-786_full#RC/Helitron", 0, '9_sample21')

old_9_22[39,] <- c("dvir-1002_full#RC/Helitron", 0, '9_sample22')
old_9_22[40,] <- c("dvir-390_full#LINE/CR1", 0, '9_sample22')
old_9_22[41,] <- c("dvir-593_full#DNA/TcMar-Tc1", 0, '9_sample22')
old_9_22[42,] <- c("dvir-786_full#RC/Helitron", 0, '9_sample22')
old_9_22[43,] <- c("NewTE7_full#DNA/TcMar-Tc1", 0, '9_sample22')
old_9_22[44,] <- c("Polyphemus_full#DNA/TcMar-Tc1", 0, '9_sample22')

old_9_23[43,] <- c("dvir-593_full#DNA/TcMar-Tc1", 0, '9_sample23')
old_9_23[44,] <- c("NewTE7_full#DNA/TcMar-Tc1", 0, '9_sample23')

old_9_21 <- old_9_21[order(old_9_21$TE),]
old_9_22 <- old_9_22[order(old_9_22$TE),]
old_9_23 <- old_9_23[order(old_9_23$TE),]

old_9 = old_9_21
old_9$strain <- NULL
old_9$mRNA_old_9 = (as.numeric(as.character(old_9_21$mRNA)) + 
                          as.numeric(as.character(old_9_22$mRNA)) +
                          as.numeric(as.character(old_9_23$mRNA))) /3 
old_9$mRNA <- NULL

old_160_25 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096025_trimmed.fq.bam', '160_25', 11698)
old_160_27 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096027_trimmed.fq.bam', '160_27', 11259)
old_160_28 <- bam_analysis('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/SRR2096028_trimmed.fq.bam', '160_28', 12572)

vector[!vector %in% old_160_28$TE]
old_160_25[41,] <- c("dvir-552_full#LINE/Jockey", 0, '160_25')
old_160_25[42,] <- c("dvir-690_full#LINE/Jockey", 0, '160_25')
old_160_25[43,] <- c("dvir-786_full#RC/Helitron", 0, '160_25')
old_160_25[44,] <- c("NewTE7_full#DNA/TcMar-Tc1", 0, '160_25')

old_160_27[40,] <- c("dvir-1002_full#RC/Helitron", 0, '160_27')
old_160_27[41,] <- c("dvir-390_full#LINE/CR1", 0, '160_27')
old_160_27[42,] <- c("dvir-593_full#DNA/TcMar-Tc1", 0, '160_27')
old_160_27[43,] <- c("dvir-690_full#LINE/Jockey", 0, '160_27')
old_160_27[44,] <- c("dvir-786_full#RC/Helitron", 0, '160_27')

old_160_28[41,] <- c("dvir-390_full#LINE/CR1", 0, '160_28')
old_160_28[42,] <- c("dvir-552_full#LINE/Jockey", 0, '160_28')
old_160_28[43,] <- c("dvir-942_full#LTR/Copia", 0, '160_28')
old_160_28[44,] <- c("NewTE7_full#DNA/TcMar-Tc1", 0, '160_28')

old_160_25 <- old_160_25[order(old_160_25$TE),]
old_160_27 <- old_160_27[order(old_160_27$TE),]
old_160_28 <- old_160_28[order(old_160_28$TE),]

old_160 = old_160_25
old_160$strain <- NULL
old_160$mRNA_old_160 = (as.numeric(as.character(old_160_25$mRNA)) + 
                      as.numeric(as.character(old_160_27$mRNA)) +
                      as.numeric(as.character(old_160_28$mRNA))) /3 
old_160$mRNA <- NULL


old = old_9
old$old_160_mRNA = old_160$mRNA_old_160
write.csv(old, '~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_44/TEs_old_RPM.csv')
old = old[old$mRNA_old_9 != 0 & old$old_160_mRNA !=0,]
library(ggrepel)

ggplot(data = old, aes(x = log10(mRNA_old_9), y = log10(old_160_mRNA))) + 
  geom_point(color = ifelse(log10(old$mRNA_old_9)/log10(old$old_160_mRNA) > 2 | log10(old$mRNA_old_9)/log10(old$old_160_mRNA) < 0.5, 'red', 'black')) +
  xlab('Old 9, log10 RPM') + ylab('Old 160, log10 RPM')+
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  #geom_text_repel(aes(label=unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),hjust=0, vjust=0) +
  geom_label_repel(aes(label = unlist(lapply(strsplit(as.character(TE), '_'), function(x)x[1]))),
                   box.padding   = 0.2, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic() + theme(text = element_text(size = 20))  + ggtitle('TEs expression level in F0 embryos\nfrom old (15-25 days) parents')   



#######################################################9 stareeet:
#DE genes between young and old 9 strain:

setwd('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8)]
setwd('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_ref/')
featurecounts2 = read.table('counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts2 = featurecounts2[,c(6, 7, 8)]
featurecounts <- featurecounts[order(rownames(featurecounts)),]
featurecounts2 <- featurecounts2[order(rownames(featurecounts2)),]

featurecounts <- cbind(featurecounts, featurecounts2)


colData <- data.frame("replicate" = c("R1","R2","R3", "R1","R2", "R3"), 
                      "treatment" = c("Young","Young", "Young", "Old", "Old", "Old"))
attributes(colData)$row.names <- c("YOUNG_R1","YOUNG_R2", "YOUNG_R3", "OLD_R1","OLD_R2", "OLD_R3")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Young")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Old","Young"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  
expressed_and_passed <- na.omit(expressed_and_passed)

write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_F0_JB/9strain_15_DE_genes_old_vs_young.csv', quote = F)
expressed_and_passed <- res2

GO <- expressed_and_passed
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(rownames(GO), ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym
GO.ENSEMBL <- bitr(GO[GO$log2FoldChange >0,]$dmel.sym, fromType = "SYMBOL",
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
View(as.data.frame(ego2))
length(GO.ENSEMBL$ENSEMBL)

View(as.data.frame(GO))

expressed_and_passed




ggplot(data = expressed_and_passed, aes(x = log2FoldChange, y = -log10(pvalue), col = log10(baseMean))) +
  theme_bw() + geom_point() + 
  geom_label_repel(aes(label = ifelse(GO$dmel.sym != 'NULL', GO$dmel.sym, rownames(expressed_and_passed))),
                   box.padding   = 0.2, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  xlab('log2FC, ref = 9')



table(twenty %in% rownames(expressed_and_passed))
twenty = c("FBgn0002781", "FBgn0053892", "FBgn0053886", "FBgn0053906",
           "FBgn0053908", "FBgn0053910", "FBgn0053880", "FBgn0053896",
           "FBgn0053902", "FBgn0053874", "FBgn0053894", "FBgn0053900",
           "FBgn0061209", "FBgn0053882", "FBgn0053888", "FBgn0053904",
           "FBgn0053884", "FBgn0053868", "FBgn0053876", "FBgn0053878",
           "FBgn0053872", "FBgn0053890", "FBgn0053870", "FBgn0053898")




###############Что меняется у 160й при старении?

setwd('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(9, 10, 11)]
setwd('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_ref/')
featurecounts2 = read.table('counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts2 = featurecounts2[,c(9, 10, 11)]
featurecounts <- featurecounts[order(rownames(featurecounts)),]
featurecounts2 <- featurecounts2[order(rownames(featurecounts2)),]
featurecounts <- cbind(featurecounts, featurecounts2)

colData <- data.frame("replicate" = c("R1","R2", "R3", "R1","R2", "R3"), 
                      "treatment" = c("Young","Young", "Young", "Old", "Old", "Old"))
attributes(colData)$row.names <- c("YOUNG_R1","YOUNG_R2","YOUNG_R3","OLD_R1","OLD_R2", "OLD_R3")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Young")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Old","Young"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  
expressed_and_passed <- na.omit(expressed_and_passed)
write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_F0_JB/160_strain_old_vs_young_DE_genes.csv')

GO <- expressed_and_passed
nrow(GO[GO$log2FoldChange < 0,])

library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(rownames(GO), ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym
GO[GO$dmel_ort == 'FBgn0053868',]
GO


GO.ENSEMBL <- bitr(GO[GO$log2FoldChange <0,]$dmel.sym, fromType = "SYMBOL",
                   toType = c("ENSEMBL"),
                   OrgDb = org.Dm.eg.db)

table(twenty %in% GO.ENSEMBL$ENSEMBL)

ego2 <- enrichGO(gene         = GO.ENSEMBL$ENSEMBL,
                 OrgDb         = org.Dm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

ego <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(ego2, showCategory=30)
goplot(ego2)
View(as.data.frame(ego2))
length(GO.ENSEMBL$ENSEMBL)

View(as.data.frame(GO))


table(twenty %in% rownames(res2))
table(twenty %in% GO.ENSEMBL$ENSEMBL)



twenty
GO.ENSEMBL$ENSEMBL




##########################################################################################################
########Посмотрим, увеличивается ли пиРНК в 9й линии со старением. Но данных по пиРНК нет, придется смотреть экспрессию TEs.

TEs = young_df
TEs$old_9 = old$mRNA_old_9
TEs$old_160 = old$old_160_mRNA
TEs$label = ifelse(as.numeric(as.character(TEs$mRNA_young_160)) > 1000, as.character(TEs$TE), NA)
TEs$label

ggplot(data = TEs, aes(x = log10(mRNA_young_160), y = log10(old_160))) +
  geom_point() + 
  xlab('Young, log10 RPM') + ylab('Old, log10 RPM') +
  geom_abline(slope = 1, intercept = 0, col = 'red', linetype = 'dashed') +
  geom_label_repel(aes(label=unlist(lapply(strsplit(as.character(TEs$TE), '_'), function(x)x[1]))),
                   box.padding   = 0.2, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  ggtitle('The level of TEs expression in 160 strain') + theme_bw() + theme(text = element_text(size = 20)) 


ggplot(data = TEs, aes(x = old_9/mRNA_young_9, y = old_160/mRNA_young_160, col = log10(mRNA_young_160))) +
  geom_point() + xlab('9 strain mRNA old/young') + ylab('160 strain mRNA, old/young')+
  geom_abline(slope = 0, intercept = 1, col = 'red', linetype = 'dashed') +
  geom_vline(xintercept = 1, linetype = 'dashed', col = 'red') +
  geom_label_repel(aes(label=unlist(lapply(strsplit(as.character(label), '_'), function(x)x[1]))),
                   hjust=0, vjust=0, angle = +40,
                   box.padding = 0.2,
                   point.padding = 0.5) + theme_bw() + ylim(c(0,2))


##################################






#F2

setwd('~/2021_year/raw_RNA_data_F2_JB/young_embryo/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11, 12, 13)]
colData <- data.frame("replicate" = c("R1","R2", "R3", "R4", "R1","R2", "R3", "R4"), 
                      "treatment" = c("Dys","Dys","Dys", "Dys", "Rec", "Rec", "Rec", "Rec"))
attributes(colData)$row.names <- c("Dys_YOUNG_R1","Dys_YOUNG_R2", "Dys_YOUNG_R3", "Dys_YOUNG_R4",
                                   "Rec_YOUNG_R1","Rec_YOUNG_R2", "Rec_YOUNG_R3", "Rec_YOUNG_R4")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Rec")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Dys","Rec"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  
write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_data_F2_JB/young_embryo/aln_to_ref/DE_genes_young_embryo_ref_reciprocal.csv', quote = F)

###########################################################



setwd('~/2021_year/raw_RNA_data_F2_JB/old_embryo/trimmed/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9, 10, 11, 12, 13)]
colData <- data.frame("replicate" = c("R1","R2", "R3", "R4", "R1","R2", "R3", "R4"), 
                      "treatment" = c("Dys","Dys","Dys", "Dys", "Rec", "Rec", "Rec", "Rec"))
attributes(colData)$row.names <- c("Dys_R1","Dys_R2", "Dys_R3", "Dys_R4",
                                   "Rec_R1","Rec_R2", "Rec_R3", "Rec_R4")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Rec")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Dys","Rec"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,] 

write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_data_F2_JB/old_embryo/trimmed/DE_genes_ref_reciprocal.csv', quote = F)

#Dysgenic embryo stareet

setwd('~/2021_year/raw_RNA_data_F2_JB/young_embryo/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(6, 7, 8, 9)]
setwd('~/2021_year/raw_RNA_data_F2_JB/old_embryo/trimmed/')
featurecounts2 = read.table('aln_to_ref/counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts2 = featurecounts2[,c(6, 7, 8, 9)]
featurecounts <- featurecounts[order(rownames(featurecounts)),]
featurecounts2 <- featurecounts2[order(rownames(featurecounts2)),]

featurecounts <- cbind(featurecounts, featurecounts2)


colData <- data.frame("replicate" = c("R1","R2","R3", "R4", "R1","R2", "R3", "R4"), 
                      "treatment" = c("Young","Young", "Young", "Young", "Old", "Old", "Old", "Old"))
attributes(colData)$row.names <- c("YOUNG_R1","YOUNG_R2", "YOUNG_R3", "YOUNG_R4", "OLD_R1","OLD_R2", "OLD_R3", "OLD_R4")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Young")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Old","Young"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  
expressed_and_passed <- na.omit(expressed_and_passed)
write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_data_F2_JB/Dysgenic_stareet_DE_genes.csv', quote = F)


#Recirrocal embryo stareet

setwd('~/2021_year/raw_RNA_data_F2_JB/young_embryo/aln_to_ref/')
library('DESeq2')
featurecounts = read.table('Counts.out', header=TRUE, sep="\t", row.names=1 )
featurecounts = featurecounts[,c(10,11, 12, 13)]
setwd('~/2021_year/raw_RNA_data_F2_JB/old_embryo/trimmed/')
featurecounts2 = read.table('aln_to_ref/counts_old.out', header=TRUE, sep="\t", row.names=1 )
featurecounts2 = featurecounts2[,c(10, 11, 12, 13)]
featurecounts <- featurecounts[order(rownames(featurecounts)),]
featurecounts2 <- featurecounts2[order(rownames(featurecounts2)),]

featurecounts <- cbind(featurecounts, featurecounts2)


colData <- data.frame("replicate" = c("R1","R2","R3", "R4", "R1","R2", "R3", "R4"), 
                      "treatment" = c("Young","Young", "Young", "Young", "Old", "Old", "Old", "Old"))
attributes(colData)$row.names <- c("YOUNG_R1","YOUNG_R2", "YOUNG_R3", "YOUNG_R4", "OLD_R1","OLD_R2", "OLD_R3", "OLD_R4")
#analyze
dds <- DESeqDataSetFromMatrix(countData = featurecounts, colData = colData, design = ~ replicate + treatment)
dds$treatment <- relevel(dds$treatment, ref = "Young")
dds$replicate <- relevel(dds$replicate, ref = "R1")
keep50 <- rowMeans(counts(dds)) >= 50
dds50 <- dds[keep50,]
dds50 <- DESeq(dds50)
res2 <- results(dds50,contrast=c("treatment","Old","Young"))
res2 <- as.data.frame(res2)
expressed_and_passed <- res2[res2$padj < 0.05,]  
expressed_and_passed <- na.omit(expressed_and_passed)

write.csv(as.data.frame(expressed_and_passed), '~/2021_year/raw_RNA_data_F2_JB/Reciprical_stareet.csv', quote = F)






table(rownames(expressed_and_passed_F2) %in% rownames(expressed_and_passed))
table(rownames(expressed_and_passed) %in% rownames(expressed_and_passed_F2))


genes = expressed_and_passed[rownames(expressed_and_passed) %in% rownames(expressed_and_passed_F2),]
genes_2 = expressed_and_passed_F2[rownames(expressed_and_passed_F2) %in% rownames(expressed_and_passed),]
genes = genes[order(rownames(genes)),]
genes_2 = genes_2[order(rownames(genes_2)),]
x = data.frame(parents = genes$log2FoldChange, children = genes_2$log2FoldChange) 
ggplot(data = x, aes(x = parents, y = children)) + geom_point()





#normalized = as.data.frame(counts(dds50, normalized = TRUE))
##normalized = normalized[rownames(normalized) %in% rownames(parents),]
#parents = parents[rownames(parents) %in% rownames(normalized),]
#parents = parents[order(rownames(parents)),]
#normalized = normalized[order(rownames(normalized)),]
#a = cbind(parents, normalized)
#normalized = as.data.frame(a)
#normalized = as.data.frame(t(normalized))
#normalized
#normalized$sample = as.factor(colnames(normalized))
#pca_res <- prcomp(normalized[1:6920,], scale. = TRUE)

#normalized = normalized[normalized$sample %in% c('9_YOUNG_R1', '160_YOUNG_R1', 'R1.dys', 'R2.dys', 'R3.dys', 'R4.dys', 'R1.rec', 'R2.rec', 'R3.rec', 'R4.rec'),]
#autoplot(pca_res, data = normalized) + theme_bw() +ggtitle('PCA plot, young and old flies')







#Results:

Young_F0 = read.csv('~/2021_year/raw_RNA_F0_JB/young_embr/trimmed/aln_to_ref/DE_F0_YOUNG.csv')
nrow(Young_F0[Young_F0$log2FoldChange < 0,])
Old_F0 = read.csv('~/2021_year/raw_RNA_F0_JB/old_embr/trimmed/aln_to_ref/DE_F0_OLD.csv')

Nine_stareet = read.csv('~/2021_year/raw_RNA_F0_JB/9strain_15_DE_genes_old_vs_young.csv')
P_stareet = read.csv('~/2021_year/raw_RNA_F0_JB/160_strain_old_vs_young_DE_genes.csv')

Young_embryo = read.csv('~/2021_year/raw_RNA_data_F2_JB/young_embryo/aln_to_ref/DE_genes_young_embryo_ref_reciprocal.csv')
Old_embryo = read.csv('~/2021_year/raw_RNA_data_F2_JB/old_embryo/trimmed/DE_genes_ref_reciprocal.csv')

Dysgenic_stareet = read.csv('~/2021_year/raw_RNA_data_F2_JB/Dysgenic_stareet_DE_genes.csv')
Reciprocal_stareet = read.csv('~/2021_year/raw_RNA_data_F2_JB/Reciprical_stareet.csv')


GO <- Reciprocal_stareet
library(org.Dm.eg.db)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
GO$dmel_ort <- ort[match(GO$X, ort$V6),]$V1
GO$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(GO$dmel_ort), 'SYMBOL','FLYBASE'))
library(clusterProfiler)
GO$dmel.sym

#Down:
GO.ENSEMBL <- bitr(GO[GO$log2FoldChange >0,]$dmel.sym, fromType = "SYMBOL",
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
goplot(ego2)


genes_that_goes_down_with_age_P = P_stareet[P_stareet$log2FoldChange < 0,]






