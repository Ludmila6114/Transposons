ref = rtracklayer::import('~/Anopheles/reference/GCF_000005575.2_AgamP3_genomic.gtf')
genome = readDNAStringSet('~/Anopheles/reference/GCF_000005575.2_AgamP3_genomic.fasta')
genome@ranges@NAMES = unlist(lapply(strsplit(genome@ranges@NAMES, ' '), function(x)x[1]))
ref$sequence = getSeq(genome, ref)
ref = ref[ref$type == 'gene',]

miRNA = ref[ref$gene_biotype == 'miRNA',]$sequence
names(miRNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'miRNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(miRNA, '~/Anopheles/reference/miRNA.fasta')

rRNA = ref[ref$gene_biotype == 'rRNA',]$sequence
names(rRNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'rRNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(rRNA, '~/Anopheles/reference/rRNA.fasta')

snoRNA = ref[ref$gene_biotype == 'snoRNA',]$sequence
names(snoRNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'snoRNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(snoRNA, '~/Anopheles/reference/snoRNA.fasta')

snRNA = ref[ref$gene_biotype == 'snRNA',]$sequence
names(snRNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'snRNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(snRNA, '~/Anopheles/reference/snRNA.fasta')

tRNA = ref[ref$gene_biotype == 'tRNA',]$sequence
names(tRNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'tRNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(tRNA, '~/Anopheles/reference/tRNA.fasta')

RNase_MRP_RNA = ref[ref$gene_biotype == 'RNase_MRP_RNA',]$sequence
names(RNase_MRP_RNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'RNase_MRP_RNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(RNase_MRP_RNA, '~/Anopheles/reference/RNase_MRP_RNA.fasta')

SRP_RNA = ref[ref$gene_biotype == 'SRP_RNA',]$sequence
names(SRP_RNA) = unlist(lapply(strsplit(ref[ref$gene_biotype == 'SRP_RNA',]$gene_id, '_'), function(x)x[2]))
writeXStringSet(SRP_RNA, '~/Anopheles/reference/SRP_RNA.fasta')

#extract CDS:
ref = ref[ref$type == 'CDS',]
ref_seq = ref$sequence
names(ref_seq) = paste(ref$gene_id,seqnames(ref),':',start(ref),'-',end(ref))
writeXStringSet(ref_seq, '~/Anopheles/reference/CDS_reference.fasta')



#extract UTRs, 300bp, roughly

genome = readDNAStringSet('~/Anopheles/reference/GCF_000005575.2_AgamP3_genomic.fasta')
genome@ranges@NAMES = unlist(lapply(strsplit(genome@ranges@NAMES, ' '), function(x)x[1]))

ref = rtracklayer::import('~/Anopheles/reference/GCF_000005575.2_AgamP3_genomic.gtf')
ref = ref[ref$type == 'CDS',]
ref$contig_length = genome[match(as.character(seqnames(ref)), genome@ranges@NAMES),]@ranges@width
ref = ref[ref$contig_length - as.numeric(as.character(end(ref))) >= 301 & as.numeric(as.character(start(ref))) >= 301,]

UTR1 = ref
UTR2 = ref
tmp = start(UTR1)
start(UTR1) = start(UTR1) - 300
end(UTR1) = tmp
UTR1$sequence = getSeq(genome, UTR1)
t = end(UTR2)
end(UTR2) = end(UTR2) + 300
start(UTR2) = t
UTR2$sequence = getSeq(genome, UTR2)

first_UTR = UTR1$sequence
names(first_UTR) = paste(seqnames(UTR1),':',start(UTR1),':',end(UTR1))
writeXStringSet(first_UTR, '~/Anopheles/reference/first_UTR_extracted.fasta')

second_UTR = UTR2$sequence
names(second_UTR) = paste(seqnames(UTR2),':',start(UTR2),':',end(UTR2))
writeXStringSet(second_UTR, '~/Anopheles/reference/second_UTR_extracted.fasta')


#initial length distribution

library('ShortRead')
male_head_1 = readFastq('~/Anopheles/trimmed_libraries/male_head_rep1_trimmed.fastq')
male_head_1 = male_head_1[male_head_1@sread@ranges@width < 35,]
male_head_2 = readFastq('~/Anopheles/trimmed_libraries/male_head_rep2_trimmed.fastq')
male_head_2 = male_head_2[male_head_2@sread@ranges@width < 35,]
male_head_3 = readFastq('~/Anopheles/trimmed_libraries/male_head_rep3_trimmed.fastq')
male_head_3 = male_head_3[male_head_3@sread@ranges@width < 35,]

male_tor_1 = readFastq('~/Anopheles/trimmed_libraries/male_thorax_rep1_trimmed.fastq')
male_tor_1 = male_tor_1[male_tor_1@sread@ranges@width < 35,]
male_tor_2 = readFastq('~/Anopheles/trimmed_libraries/male_thorax_rep2_trimmed.fastq')
male_tor_2 = male_tor_2[male_tor_2@sread@ranges@width < 35,]
male_tor_3 = readFastq('~/Anopheles/trimmed_libraries/male_thorax_rep3_trimmed.fastq')
male_tor_3 = male_tor_3[male_tor_3@sread@ranges@width < 35,]


female_head_1 = readFastq('~/Anopheles/trimmed_libraries/female_head_rep1_trimmed.fastq')
female_head_1 = female_head_1[female_head_1@sread@ranges@width < 35,]
female_head_2 = readFastq('~/Anopheles/trimmed_libraries/female_head_rep2_trimmed.fastq')
female_head_2 = female_head_2[female_head_2@sread@ranges@width < 35,]
female_head_3 = readFastq('~/Anopheles/trimmed_libraries/female_head_rep3_trimmed.fastq')
female_head_3 = female_head_3[female_head_3@sread@ranges@width < 35,]

female_tor_1 = readFastq('~/Anopheles/trimmed_libraries/female_thorax_rep1_trimmed.fastq')
female_tor_1 = female_tor_1[female_tor_1@sread@ranges@width < 35,]
female_tor_2 = readFastq('~/Anopheles/trimmed_libraries/female_thorax_rep2_trimmed.fastq')
female_tor_2 = female_tor_2[female_tor_2@sread@ranges@width < 35,]
female_tor_3 = readFastq('~/Anopheles/trimmed_libraries/female_thorax_rep3_trimmed.fastq')
female_tor_3 = female_tor_3[female_tor_3@sread@ranges@width < 35,]


library(ggplot2)
mh1 = ggplot(data = as.data.frame(male_head_1@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep1') + theme(text = element_text(size = 15))
mh2 = ggplot(data = as.data.frame(male_head_2@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep2') + theme(text = element_text(size = 15))
mh3 = ggplot(data = as.data.frame(male_head_3@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep3') + theme(text = element_text(size = 15))
mt1 = ggplot(data = as.data.frame(male_tor_1@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep1') + theme(text = element_text(size = 15))
mt2 = ggplot(data = as.data.frame(male_tor_2@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep2') + theme(text = element_text(size = 15))
mt3 = ggplot(data = as.data.frame(male_tor_3@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep3') + theme(text = element_text(size = 15))


fh1 = ggplot(data = as.data.frame(female_head_1@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep1') + theme(text = element_text(size = 15))
fh2 = ggplot(data = as.data.frame(female_head_2@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep2') + theme(text = element_text(size = 15))
fh3 = ggplot(data = as.data.frame(female_head_3@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep3') + theme(text = element_text(size = 15))
ft1 = ggplot(data = as.data.frame(female_tor_1@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep1') + theme(text = element_text(size = 15))
ft2 = ggplot(data = as.data.frame(female_tor_2@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep2') + theme(text = element_text(size = 15))
ft3 = ggplot(data = as.data.frame(female_tor_3@sread@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep3') + theme(text = element_text(size = 15))

mt3
library(gridExtra)
grid.arrange(mh1, mh2, mh3, mt1, mt2, mt3, fh1, fh2, fh3, ft1, ft2, ft3, nrow = 4, ncol = 3)


#aligned to base
mh1 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_head_rep1_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
mh2 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_head_rep2_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
mh3 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_head_rep3_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
mt1 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_thorax_rep1_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
mt2 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_thorax_rep2_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
mt3 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/male_thorax_rep3_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')

fh1 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_head_rep1_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
fh2 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_head_rep2_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
fh3 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_head_rep3_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
ft1 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_thorax_rep1_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
ft2 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_thorax_rep2_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')
ft3 = readDNAStringSet('~/Anopheles/bowtie_to_xRNA_base/aligned_to_base/female_thorax_rep3_trimmed.fastq.ALL_ALIGNED_AND_NO.sam.fasta')


pic_mh1 = ggplot(as.data.frame(mh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep1, xRNAs') + theme(text = element_text(size = 15))
pic_mh2 = ggplot(as.data.frame(mh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep2, xRNAs') + theme(text = element_text(size = 15))
pic_mh3 = ggplot(as.data.frame(mh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep3, xRNAs') + theme(text = element_text(size = 15))
pic_mt1 = ggplot(as.data.frame(mt1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep1, xRNAs') + theme(text = element_text(size = 15))
pic_mt2 = ggplot(as.data.frame(mt2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep2, xRNAs') + theme(text = element_text(size = 15))
pic_mt3 = ggplot(as.data.frame(mt3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep3, xRNAs') + theme(text = element_text(size = 15))

pic_fh1 = ggplot(as.data.frame(fh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep1, xRNAs') + theme(text = element_text(size = 15))
pic_fh2 = ggplot(as.data.frame(fh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep2, xRNAs') + theme(text = element_text(size = 15))
pic_fh3 = ggplot(as.data.frame(fh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep3, xRNAs') + theme(text = element_text(size = 15))
pic_ft1 = ggplot(as.data.frame(ft1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep1, xRNAs') + theme(text = element_text(size = 15))
pic_ft2 = ggplot(as.data.frame(ft2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep2, xRNAs') + theme(text = element_text(size = 15))
pic_ft3 = ggplot(as.data.frame(ft3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep3, xRNAs') + theme(text = element_text(size = 15))


grid.arrange(pic_mh1, pic_mh2, pic_mh3, pic_mt1, pic_mt2, pic_mt3, pic_fh1, pic_fh2, pic_fh3, pic_ft1, pic_ft2, pic_ft3, nrow = 4, ncol = 3)


#miRNA sequences:
mh1 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_head_rep1_filtered.fasta.sam.fasta')
mh2 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_head_rep2_filtered.fasta.sam.fasta')
mh3 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_head_rep3_filtered.fasta.sam.fasta')
mt1 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_thorax_rep1_filtered.fasta.sam.fasta')
mt2 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_thorax_rep2_filtered.fasta.sam.fasta')
mt3 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/male_thorax_rep3_filtered.fasta.sam.fasta')

fh1 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_head_rep1_filtered.fasta.sam.fasta')
fh2 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_head_rep2_filtered.fasta.sam.fasta')
fh3 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_head_rep3_filtered.fasta.sam.fasta')
ft1 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_thorax_rep1_filtered.fasta.sam.fasta')
ft2 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_thorax_rep2_filtered.fasta.sam.fasta')
ft3 = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/female_thorax_rep3_filtered.fasta.sam.fasta')


pic_mh1 = ggplot(as.data.frame(mh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male head rep1, miRNA') + theme(text = element_text(size = 15))
pic_mh2 = ggplot(as.data.frame(mh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male head rep2, miRNA') + theme(text = element_text(size = 15))
pic_mh3 = ggplot(as.data.frame(mh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male head rep3, miRNA') + theme(text = element_text(size = 15))
pic_mt1 = ggplot(as.data.frame(mt1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_mt2 = ggplot(as.data.frame(mt2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_mt3 = ggplot(as.data.frame(mt3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Male thorax rep3, miRNA') + theme(text = element_text(size = 15))

pic_fh1 = ggplot(as.data.frame(fh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female head rep1, miRNA') + theme(text = element_text(size = 15))
pic_fh2 = ggplot(as.data.frame(fh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female head rep2, miRNA') + theme(text = element_text(size = 15))
pic_fh3 = ggplot(as.data.frame(fh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female head rep3, miRNA') + theme(text = element_text(size = 15))
pic_ft1 = ggplot(as.data.frame(ft1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_ft2 = ggplot(as.data.frame(ft2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_ft3 = ggplot(as.data.frame(ft3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA length') + ylab('Density') + ggtitle('Female thorax rep3, miRNA') + theme(text = element_text(size = 15))

pic_ft1
grid.arrange(pic_mh1, pic_mh2, pic_mh3, pic_mt1, pic_mt2, pic_mt3, pic_fh1, pic_fh2, pic_fh3, pic_ft1, pic_ft2, pic_ft3, nrow = 4, ncol = 3)



#notmiRNA:
setwd('/mnt/raid/pro_milka/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/notmiRNA/fasta_notmirna')
mh1 = readDNAStringSet('male_head_rep1_filtered.fasta.sam.fasta')
mh1 = mh1[mh1@ranges@width < 35,]
mh2 = readDNAStringSet('male_head_rep2_filtered.fasta.sam.fasta')
mh2 = mh2[mh2@ranges@width < 35,]
mh3 = readDNAStringSet('male_head_rep3_filtered.fasta.sam.fasta')
mh3 = mh3[mh3@ranges@width < 35,]
mt1 = readDNAStringSet('male_thorax_rep1_filtered.fasta.sam.fasta')
mt1 = mt1[mt1@ranges@width < 35,]
mt2 = readDNAStringSet('male_thorax_rep2_filtered.fasta.sam.fasta')
mt2 = mt2[mt2@ranges@width < 35,]
mt3 = readDNAStringSet('male_thorax_rep3_filtered.fasta.sam.fasta')
mt3 = mt3[mt3@ranges@width < 35,]


fh1 = readDNAStringSet('female_head_rep1_filtered.fasta.sam.fasta')
fh1 = fh1[fh1@ranges@width < 35,]
fh2 = readDNAStringSet('female_head_rep2_filtered.fasta.sam.fasta')
fh2 = fh2[fh2@ranges@width < 35,]
fh3 = readDNAStringSet('female_head_rep3_filtered.fasta.sam.fasta')
fh3 = fh3[fh3@ranges@width < 35,]
ft1 = readDNAStringSet('female_thorax_rep1_filtered.fasta.sam.fasta')
ft1 = ft1[ft1@ranges@width < 35,]
ft2 = readDNAStringSet('female_thorax_rep2_filtered.fasta.sam.fasta')
ft2 = ft2[ft2@ranges@width < 35,]
ft3 = readDNAStringSet('female_thorax_rep3_filtered.fasta.sam.fasta')
ft3 = ft3[ft3@ranges@width < 35,]


pic_mh1 = ggplot(as.data.frame(mh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep1, small RNAs') + theme(text = element_text(size = 15))
pic_mh2 = ggplot(as.data.frame(mh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep2, small RNAs') + theme(text = element_text(size = 15))
pic_mh3 = ggplot(as.data.frame(mh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male head rep3, small RNAs') + theme(text = element_text(size = 15))
pic_mt1 = ggplot(as.data.frame(mt1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep1, small RNAs') + theme(text = element_text(size = 15))
pic_mt2 = ggplot(as.data.frame(mt2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep2, small RNAs') + theme(text = element_text(size = 15))
pic_mt3 = ggplot(as.data.frame(mt3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Male thorax rep3, small RNAs') + theme(text = element_text(size = 15))

pic_fh1 = ggplot(as.data.frame(fh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep1, small RNAs') + theme(text = element_text(size = 15))
pic_fh2 = ggplot(as.data.frame(fh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep2, small RNAs') + theme(text = element_text(size = 15))
pic_fh3 = ggplot(as.data.frame(fh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female head rep3, small RNAs') + theme(text = element_text(size = 15))
pic_ft1 = ggplot(as.data.frame(ft1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep1, small RNAs') + theme(text = element_text(size = 15))
pic_ft2 = ggplot(as.data.frame(ft2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep2, small RNAs') + theme(text = element_text(size = 15))
pic_ft3 = ggplot(as.data.frame(ft3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('small RNA length') + ylab('Density') + ggtitle('Female thorax rep3, small RNAs') + theme(text = element_text(size = 15))


grid.arrange(pic_mh1, pic_mh2, pic_mh3, pic_mt1, pic_mt2, pic_mt3, pic_fh1, pic_fh2, pic_fh3, pic_ft1, pic_ft2, pic_ft3, nrow = 4, ncol = 3)

######


#miRNA aligned to lncRNA:

setwd('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/align_to_lncRNA_important_without_flag/aligned_F4/')
mh1 = readDNAStringSet('male_head_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mh2 = readDNAStringSet('male_head_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mh3 = readDNAStringSet('male_head_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt1 = readDNAStringSet('male_thorax_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt2 = readDNAStringSet('male_thorax_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt3 = readDNAStringSet('male_thorax_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')

fh1 = readDNAStringSet('female_head_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
fh2 = readDNAStringSet('female_head_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
fh3 = readDNAStringSet('female_head_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft1 = readDNAStringSet('female_thorax_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft2 = readDNAStringSet('female_thorax_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft3 = readDNAStringSet('female_thorax_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')


pic_mh1 = ggplot(as.data.frame(mh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male head rep1, miRNA') + theme(text = element_text(size = 15))
pic_mh2 = ggplot(as.data.frame(mh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male head rep2, miRNA') + theme(text = element_text(size = 15))
pic_mh3 = ggplot(as.data.frame(mh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male head rep3, miRNA') + theme(text = element_text(size = 15))
pic_mt1 = ggplot(as.data.frame(mt1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_mt2 = ggplot(as.data.frame(mt2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_mt3 = ggplot(as.data.frame(mt3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep3, miRNA') + theme(text = element_text(size = 15))

pic_fh1 = ggplot(as.data.frame(fh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female head rep1, miRNA') + theme(text = element_text(size = 15))
pic_fh2 = ggplot(as.data.frame(fh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female head rep2, miRNA') + theme(text = element_text(size = 15))
pic_fh3 = ggplot(as.data.frame(fh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female head rep3, miRNA') + theme(text = element_text(size = 15))
pic_ft1 = ggplot(as.data.frame(ft1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_ft2 = ggplot(as.data.frame(ft2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_ft3 = ggplot(as.data.frame(ft3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep3, miRNA') + theme(text = element_text(size = 15))

grid.arrange(pic_mh1, pic_mh2, pic_mh3, pic_mt1, pic_mt2, pic_mt3, pic_fh1, pic_fh2, pic_fh3, pic_ft1, pic_ft2, pic_ft3, nrow = 4, ncol = 3)



#miRNA NOT aligned to lncRNA:

setwd('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/aligned_miRNA_F4/fasta_aligned_to_miRNA/align_to_lncRNA_important_without_flag/not_aligned_to_lncRNA/')
mh1 = readDNAStringSet('male_head_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mh2 = readDNAStringSet('male_head_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mh3 = readDNAStringSet('male_head_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt1 = readDNAStringSet('male_thorax_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt2 = readDNAStringSet('male_thorax_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
mt3 = readDNAStringSet('male_thorax_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')

fh1 = readDNAStringSet('female_head_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
fh2 = readDNAStringSet('female_head_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
fh3 = readDNAStringSet('female_head_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft1 = readDNAStringSet('female_thorax_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft2 = readDNAStringSet('female_thorax_rep2_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')
ft3 = readDNAStringSet('female_thorax_rep3_filtered.fasta.sam.fasta.ALL_ALIGNED_lncRNA_AND_NO.sam.fasta')


pic_mh1 = ggplot(as.data.frame(mh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male head rep1, miRNA') + theme(text = element_text(size = 15))
pic_mh2 = ggplot(as.data.frame(mh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male head rep2, miRNA') + theme(text = element_text(size = 15))
pic_mh3 = ggplot(as.data.frame(mh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male head rep3, miRNA') + theme(text = element_text(size = 15))
pic_mt1 = ggplot(as.data.frame(mt1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_mt2 = ggplot(as.data.frame(mt2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_mt3 = ggplot(as.data.frame(mt3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Male thorax rep3, miRNA') + theme(text = element_text(size = 15))

pic_fh1 = ggplot(as.data.frame(fh1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female head rep1, miRNA') + theme(text = element_text(size = 15))
pic_fh2 = ggplot(as.data.frame(fh2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female head rep2, miRNA') + theme(text = element_text(size = 15))
pic_fh3 = ggplot(as.data.frame(fh3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female head rep3, miRNA') + theme(text = element_text(size = 15))
pic_ft1 = ggplot(as.data.frame(ft1@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep1, miRNA') + theme(text = element_text(size = 15))
pic_ft2 = ggplot(as.data.frame(ft2@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep2, miRNA') + theme(text = element_text(size = 15))
pic_ft3 = ggplot(as.data.frame(ft3@ranges), aes(x = width)) + geom_density() + theme_bw() +
  xlab('miRNA NOT to lncRNA length') + ylab('Density') + ggtitle('Female thorax rep3, miRNA') + theme(text = element_text(size = 15))

grid.arrange(pic_mh1, pic_mh2, pic_mh3, pic_mt1, pic_mt2, pic_mt3, pic_fh1, pic_fh2, pic_fh3, pic_ft1, pic_ft2, pic_ft3, nrow = 4, ncol = 3)




#align pi & si RNA to a genome, overlap with gtf file:

.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
bam_to_df <- function(path){
  bam <- scanBam(path)
  bam_field <- names(bam[[1]])
  list_bam <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  bam_df <- do.call("DataFrame", list_bam)
  names(bam_df) <- bam_field
  bam_df <- data.frame(bam_df)
  #bam_df$type <- ifelse(bam_df$qwidth == 21 | bam_df$qwidth == 22, 'siRNA', 
  #                      ifelse((bam_df$qwidth >= 23 & bam_df$qwidth <= 30), 'piRNA', 'smth'))
  #bam_df <- bam_df[bam_df$flag != 4,]
  return(bam_df)
}


setwd('~/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/notmiRNA/fasta_notmirna/align_to_genome_without_flag/aligned_F4/')
female_thorax_rep1 = bam_to_df('female_thorax_rep1_filtered.fasta.sam.fasta.ALL_ALIGNED_genome_AND_NO.sam')

rna_ranges = GRanges(seqnames = female_thorax_rep1$rname, ranges = IRanges(start = female_thorax_rep1$pos, width = female_thorax_rep1$qwidth + 1))
ann = rtracklayer::import('~/Anopheles/reference/genome/GCF_000005575.2_AgamP3_genomic.gtf')

rna_annotate = subsetByOverlaps(rna_ranges, ann)
length(rna_annotate)/length(rna_ranges)*100

annotated_small = as.data.frame(subsetByOverlaps(ann, rna_ranges))
table(annotated_small$type)

df2 = NULL
df2 <- annotated_small %>% group_by(type) %>% mutate(count_name_occurr = n())
df2 = df2[complete.cases(df2$type),]

library(ggplot2)
ggplot(data = df2, aes(x = reorder(df2$type, -count_name_occurr))) + 
  geom_bar(stat = "count") +
  geom_bar(fill = 'cornflowerblue') + 
  coord_flip() + 
  theme_bw() +
  theme(text = element_text(size = 15)) +
  xlab('Feature type') + ylab('Number of alignments') + ggtitle('Alignments of piRNA and siRNA\nto AgamP3 genome')

table(df2$gene_biotype)  



#check 1U and 10A bias:


female_thorax_rep1$first_letter = unlist(lapply(strsplit(female_thorax_rep1$seq, ''), function(x)x[1]))
table(female_thorax_rep1$first_letter)


female_thorax_rep1$ten_letter = unlist(lapply(strsplit(female_thorax_rep1$seq, ''), function(x)x[10]))
table(female_thorax_rep1$ten_letter)





#check. from the initial trimmed library extract peak 20bp and align to UTR , TE and lncRNA, miRNA, one sample (female_thorax_rep1)
data = readDNAStringSet('~/Anopheles/filtered_from_xRNA_result_copy/Check_peak/female_thorax_rep1_filtered.fasta')
data = data[data@ranges@width == 20,]
writeXStringSet(data, '~/Anopheles/filtered_from_xRNA_result_copy/Check_peak/peak_20.fasta')


w = readDNAStringSet('/mnt/raid/pro_milka/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/notmiRNA/fasta_notmirna/align_to_genome_without_flag/aligned_F4/head_rep1.fasta')
table(w@ranges@width)



peak_18 = readDNAStringSet('/mnt/raid/pro_milka/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/notmiRNA/fasta_notmirna/female_head_rep1_filtered.fasta.sam.fasta')
peak_18 = peak_18[peak_18@ranges@width == 18,]
writeXStringSet(peak_18, '/mnt/raid/pro_milka/Anopheles/filtered_from_xRNA_result_copy/bowtie_to_miRNA_result/notmiRNA/fasta_notmirna/PEAK_18.fasta')




