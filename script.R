#SNPs 9 and 160

P_strain = read.table('~/soft/snpEff/160JB_annotated.vcf', skip = 62)
P_strain = P_strain[P_strain$V7 == 'PASS',]
P_strain$AF = unlist(lapply(strsplit(as.character(P_strain$V10), ':'), function(x)x[5]))
P_strain$pos = paste(P_strain$V1,'_',P_strain$V2)

annotation_Sasha <- rtracklayer::import('~/soft/snpEff/data/Dvir_160/genes.gtf')
SNPs = GRanges(seqnames = P_strain$V1, 
               ranges = IRanges(start = P_strain$V2, 
                                end = P_strain$V2 + 1))

length(SNPs)
damage = subsetByOverlaps(annotation_Sasha, SNPs, ignore.strand = TRUE)

P_strain$ANN = unlist(lapply(strsplit(as.character(P_strain$V8), ';'), function(x)x[30]))
P_strain$GenomicLoci = unlist(lapply(strsplit(as.character(P_strain$ANN), '\\|'), function(x)x[2]))
P_strain$Effect = unlist(lapply(strsplit(as.character(P_strain$ANN), '\\|'), function(x)x[3]))
P_strain$Effectic = unlist(lapply(strsplit(as.character(P_strain$ANN), '\\|'), function(x)x[6]))

table(P_strain$GenomicLoci, P_strain$Effectic, P_strain$Effect)


genes = read.csv('~/soft/snpEff/160JB_statistics.genes.txt', sep = '\t', skip = 1, header = TRUE)
table(genes$GeneId %in% damage$gene_id)
length(unique(genes$GeneId))

table(unique(genes$GeneId) %in% damage$gene_id)


M_strain = read.table('~/soft/snpEff/9JB_annotated.vcf', skip = 62)
M_strain$AF = unlist(lapply(strsplit(as.character(M_strain$V10), ':'), function(x)x[5]))
M_strain = M_strain[M_strain$V7 == 'PASS',]
M_strain$pos = paste(M_strain$V1,'_',M_strain$V2)

homo_in_M = M_strain[M_strain$AF > 0.9,]
homo_in_M = homo_in_M[!homo_in_M$pos %in% P_strain$pos,]

homo_in_M$ANN = unlist(lapply(strsplit(as.character(homo_in_M$V8), ';'), function(x)x[30]))
homo_in_M$GenomicLoci = unlist(lapply(strsplit(as.character(homo_in_M$ANN), '\\|'), function(x)x[2]))
homo_in_M$Effect = unlist(lapply(strsplit(as.character(homo_in_M$ANN), '\\|'), function(x)x[3]))
homo_in_M$Effectic = unlist(lapply(strsplit(as.character(homo_in_M$ANN), '\\|'), function(x)x[6]))

table(homo_in_M$GenomicLoci, homo_in_M$Effectic, homo_in_M$Effect)


annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
annotation_Sasha = annotation_Sasha[annotation_Sasha$type == 'gene',]
SNPs = GRanges(seqnames = homo_in_M[homo_in_M$Effectic == 'transcript',]$V1, 
               ranges = IRanges(start = homo_in_M[homo_in_M$Effectic == 'transcript',]$V2, 
                                                                                      end = homo_in_M[homo_in_M$Effectic == 'transcript',]$V2 + 1))


Genes = subsetByOverlaps(annotation_Sasha, SNPs, ignore.strand = TRUE)
length(unique(Genes$gene_id))
Genes

table(unique(Genes$gene_id) %in% genes$X.GeneName)

genes = read.csv('~/soft/snpEff/9JB_LOF_statistics.genes.txt', sep = '\t', skip = 1, header = TRUE)
ort <- read.table('~/dmel_orthologs_in_drosophila_species_fb_2019_05.tsv', fill = T, header = F)
genes$dmel_ort <- ort[match(genes$X.GeneName, ort$V6),]$V1
genes$dmel.sym <- as.character(mapIds(org.Dm.eg.db, as.character(genes$dmel_ort), 'SYMBOL','FLYBASE'))


GO.ENSEMBL <- bitr(genes$dmel.sym, fromType = "SYMBOL",
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


unlist(lapply(strsplit(as.character(homo_in_M$ANN[1]), '\\|'), function(x)x[3]))
strsplit(as.character(homo_in_M$ANN[1]), '|')
strsplit(homo_in_M$ANN[1], '|')


SNPs = GRanges(seqnames = homo_in_M[homo_in_M$Effect == 'HIGH',]$V1, ranges = IRanges(start = homo_in_M[homo_in_M$Effect == 'HIGH',]$V2, 
                                                                                      end = homo_in_M[homo_in_M$Effect == 'HIGH',]$V2 + 1))
annotation_Sasha <- rtracklayer::import('/mnt/raid/alexr/genome.annotation/FINAL/Dvir_160_liftoff_recursive.gtf')
Gene <- annotation_Sasha[annotation_Sasha$type == 'exon']

genes_with_snps = subsetByOverlaps(annotation_Sasha, SNPs, ignore.strand = TRUE)
genes_with_snps

genes_with_snps[genes_with_snps$gene_id == 'FBgn0205684',]

View(as.data.frame(findOverlaps(SNPs, reduce(Gene[Gene$gene_id == 'FBgn0204789',]))))


length(reduce(Gene[Gene$gene_id == 'FBgn0201134',]))





genes = read.csv('~/soft/snpEff/9JB_LOF_statistics.genes.txt', sep = '\t', skip = 1, header = TRUE)
genes2 = read.csv('~/soft/snpEff/9JB_statistics.genes.txt', sep = '\t', skip = 1, header = TRUE)

table(genes == genes2)



















#####################ALLELIC specific transcription:

dys_160 = read.table('~/2021_year/AL_SPECIFIC/SRR2096029_160', header=TRUE, sep="\t", row.names=1 )
dys_9 = read.table('~/2021_year/AL_SPECIFIC/SRR2096029_9.counts', header=TRUE, sep="\t", row.names=1)

dys_160 = dys_160[order(rownames(dys_160)),]
dys_9 = dys_9[order(rownames(dys_9)),]

dys_160 = dys_160[rownames(dys_160) %in% rownames(dys_9),]
dys_9 = dys_9[rownames(dys_9) %in% rownames(dys_160),]

df = data.frame(dys_160 = dys_160$SRR2096029_trimmed_to_hybrid.bam, dys_9 = dys_9$SRR2096029_trimmed_to_hybrid.bam)
df$gene = rownames(dys_160)
df = df[df$dys_160 != 0 | df$dys_9 != 0,]
df$ratio = df$dys_160/(df$dys_9 + df$dys_160)

df = df[df$dys_160 + df$dys_9 > 100,]
#df$DE_parents = ifelse(df$gene %in% rownames(expressed_and_passed), 'DE', 'not')
#ggplot(data = df, aes(x = log(dys_160), y = log(dys_9), col = DE_parents)) + geom_point() + theme_bw() 
df[df$ratio == 1,]

rec_160 = read.table('~/2021_year/AL_SPECIFIC/SRR2096033_160', header=TRUE, sep="\t", row.names=1)
rec_9 = read.table('~/2021_year/AL_SPECIFIC/SRR2096033_9.counts', header=TRUE, sep="\t", row.names=1)

rec_160 = rec_160[order(rownames(rec_160)),]
rec_9 = rec_9[order(rownames(rec_9)),]
rec_160 = rec_160[rownames(rec_160) %in% rownames(rec_9),]
rec_9 = rec_9[rownames(rec_9) %in% rownames(rec_160),]

g = data.frame(rec_160 = rec_160$SRR2096033_trimmed_to_hybrid.bam, rec_9 = rec_9$SRR2096033_trimmed_to_hybrid.bam)
g$gene = rownames(rec_160)
g = g[g$rec_160 != 0 | g$rec_9 != 0,]

g = g[g$gene %in% df$gene,]
df = df[df$gene %in% g$gene,]
g$ratio = g$rec_160/(g$rec_9 + g$rec_160)


df$rec_160 = g$rec_160
df$rec_9 = g$rec_9
df$ratio_rec = g$ratio

#Imprinting
Imprinting = df[df$ratio > 0.6 & df$ratio_rec < 0.4 | df$ratio < 0.4 & df$ratio_rec > 0.6,]

#AS:
AS = df[df$ratio > 0.7 & df$ratio_rec > 0.7 | df$ratio < 0.3 & df$ratio_rec < 0.3,]















summ = data.frame('160/(160 + 9)' = c(df$ratio, g$ratio), 
                  type = c(rep('Dysgenic embryo', 8273), rep('Reciprocal embryo', 8273)))


ggplot(data = summ, aes(x = summ$X160..160...9., fill = type)) + geom_density(alpha = 0.3) + theme_bw() + 
  xlab('Unique transcripts for 160 divided by sum unique transcripts for 160 and 9 for each gene')


summ$gene = rep(df$gene, 2)
#summ$label = ifelse(summ$gene %in% rownames(expressed_and_passed_F2), '1','0')
summ$label_parents = ifelse(summ$gene %in% rownames(expressed_and_passed), '1', '0')

#ggplot(data = summ, aes(x = gene, fill = type, y = X160..160...9., col = label)) + geom_point() + theme_bw() + 
#  xlab('Unique transcripts for 160 divided by sum unique transcripts for 160 and 9 for each gene') + facet_wrap(.~label_parents)


diference = data.frame(gene = df$gene, ratio_dys = df$ratio)
rec = data.frame(gene = df$gene, ratio_rec = g$ratio)
diference$ratio_rec = rec$ratio_rec
diference$DE_parents = ifelse(diference$gene %in% rownames(expressed_and_passed), 'DE', 'not')
diference$DE_F2 = ifelse(diference$gene %in% rownames(expressed_and_passed_F2), 'DE_F2', 'not_DE_F2')

diference = diference[round(diference$ratio_dys*10) != round(diference$ratio_rec*10),]
ggplot(data = diference, aes(x = gene, y = ratio_rec-ratio_dys, col = DE_F2)) +
  geom_point() + theme_bw()



#table(rownames(expressed_and_passed) %in% diference$gene)

#nrow(diference[diference$ratio_dys < 0.75 & diference$ratio_rec > 0.75,])


#df[df$dys_160 == 0 & df$dys_9 > 100,]
#df[df$dys_160 > 100 & df$dys_9 == 0,]

ggplot(data = df, aes(x = dys_160, y = dys_9))
                      + geom_point()








