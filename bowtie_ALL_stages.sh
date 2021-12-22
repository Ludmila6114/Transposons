#!/bin/bash
mkdir Unmapped_fasta
mkdir Unmapped_fasta/final_sam

for i in *.fq; do
	echo $i
	bowtie $i -v 3 --index /mnt/raid/pro_milka/2021_year/GENES_piRNA_automatic/All_small_RNA_except_piRNA --sam -p 64 > $i.sam
	samtools view -Sbh -f 4 $i.sam >> $i.unmapped.bam
	samtools fasta $i.unmapped.bam >> $i.smallRNA_filtered.fasta
done

mv *.fasta Unmapped_fasta
cd /mnt/raid/pro_milka/2021_year/GENES_piRNA_automatic/trimmed/Unmapped_fasta
echo 'OK'

for j in *.fasta; do
	awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) >= 23 && length(seq) <= 30) {print header, seq}}' < $j > $j.piRNA.fasta
	bowtie $j.piRNA.fasta -v 3 -f -m 1 --index /mnt/raid/pro_milka/2021_year/GENES_piRNA_automatic/reference --sam -p 64 > $j.final.sam
done
