library(GenomicRanges)
library(Biostrings)
library(BSgenome)

chromatin_test <- function(start, end, chr){
  if (chr == 'CM017604.2'){
    if (end >= 29391167){
      return('Heterochromatin')
    } else {
      return('Euchromatin')
    } 
  } else if (chr == 'CM017605.2'){
    if (start >= 381633	& end <= 33758602){
      return('Euchromatin')
    } else{
      return('Heterochromatin')
    }
  } else if (chr == 'CM017606.2'){
    if (start >= 118357	& end <= 26407003){
      return('Euchromatin')
    } else{
      return('Heterochromatin')
    }
  } else if (chr == 'CM017607.2'){
    if (end >= 28693856){
      return('Heterochromatin')
    } else{
      return('Euchromatin')
    }
  } else if (chr == 'CM017608.2'){
    if (end >= 25039679){
      return('Heterochromatin')
    } else{
      return('Euchromatin')
    }
  } else {
    return('Heterochromatin')
  }
}
blast_reduce_upd <- function(left, right){
  left <- left[as.numeric(as.character(left$V4)) >= 1000,]
  right <- right[as.numeric(as.character(right$V4)) >= 1000,]
  if (nrow(left) == 0 | nrow(right) == 0){
    return()
  } else{
    left_1 <- left[as.numeric(as.character(left$V9)) < as.numeric(as.character(left$V10)),]
    left_2 <- left[as.numeric(as.character(left$V9)) > as.numeric(as.character(left$V10)),]
    right_1 <- right[as.numeric(as.character(right$V9)) < as.numeric(as.character(right$V10)),]
    right_2 <- right[as.numeric(as.character(right$V9)) > as.numeric(as.character(right$V10)),]
    
    
    if ((nrow(left_1) == 0 & nrow(left_2) != 0 & nrow(right_1) != 0 & nrow(right_2) == 0) | (nrow(left_1) != 0 & nrow(left_2) == 0 & nrow(right_1) == 0 & nrow(right_2) != 0)){
      return('FUCK!')
    }
    else {
      if (nrow(left_1) != 0 & nrow(right_1) != 0){
        k <- 1
        kmin <- nrow(left_1) + 1
        name_for_future_blast_1 <- ''
        
        dist_1 <- vector()
        for (k in 1:nrow(left_1)){
          R_border_1 <- right_1[(as.numeric(as.character(right_1$V9)) >= as.numeric(as.character(left_1$V10[k]))) & (as.character(right_1$V2) == as.character(left_1$V2[k])),]
          if (nrow(R_border_1) == 0){
            return()
          } else {
            dist_1 <- c(dist_1, min((as.numeric(as.character(R_border_1$V9)) - as.numeric(as.character(left_1$V10[k])))))
            if (min(as.numeric(as.character(R_border_1$V9)) - as.numeric(as.character(left_1$V10[k]))) <= min(dist_1)){
              kmin <- k
              start_left_1 <- left_1$V10[kmin]
              name_for_future_blast_1 <- left_1$V2[kmin]
            }
          }
        }
        dist_1 <- min(dist_1)
        #end first if
      } 
      
      
      if (nrow(left_2) != 0 & nrow(right_2) != 0){
        m <- 1
        mmin <- nrow(right_2) + 1
        name_for_future_blast_2 <- ''
        dist_op <- vector()
        for (m in 1:nrow(right_2)){
          
          L_border_2 <- left_2[(as.numeric(as.character(left_2$V10)) >= as.numeric(as.character(right_2$V9[m]))) & (as.character(right_2$V2[m]) == as.character(left_2$V2)),]
          if (nrow(L_border_2) == 0){
            return()
          } else {
            dist_op <- c(dist_op, min(abs(as.numeric(as.character(L_border_2$V10)) - as.numeric(as.character(right_2$V9[m])))))
            if (min(as.numeric(as.character(L_border_2$V10)) - as.numeric(as.character(right_2$V9[m]))) <= min(dist_op)){
              mmin <- m
              start_op_2 <- right_2$V9[mmin]
              name_for_future_blast_2 <- right_2$V2[mmin]
            }
          }
        }
        dist_op <- min(dist_op)
      }
      
      if (exists("dist_1") & exists("dist_op")){
        if (dist_1 < dist_op) {
          return(paste(dist_1, ' ', start_left_1, ' ', name_for_future_blast_1))
        } else{
          return(paste(dist_op, ' ', start_op_2, ' ', name_for_future_blast_2))
        }
        
      } else if (exists("dist_op")){
        #print(paste(dist_op, ' ', start_op_2, ' ', name_for_future_blast_2))
        return (paste(dist_op, ' ', start_op_2, ' ', name_for_future_blast_2))
      } else if (exists("dist_1")){
        #print(paste(dist_1, ' ', start_left_1, ' ', name_for_future_blast_1))
        return (paste(dist_1, ' ', start_left_1, ' ', name_for_future_blast_1))
      } else{
        return()
      }
    }
  }
}
strange_analisys <- function(dist, cord, chr, N){
  print(paste('I will try blast with', dist, ' ', cord, ' ', chr, ' number: ', N)) 
  if (dist < 0.01*GR_NEW[N]$canon){
    print('Dist less then 1% of canon. Blast doesnt make sence')
    assign("strange_nonsence", c(strange_nonsence, N), envir = .GlobalEnv)
  } else if (dist > 100000){
    print('Dist longer then 100k. BLAST DOESNT MAKE SENCE')
    assign("strange_nonsence", c(strange_nonsence, N), envir = .GlobalEnv)
  }
  else{
    particular_chr <- DNAStringSet(genome_9[names(genome_9) == as.character(chr),])
    inside_seq <- DNAStringSet(unlist(particular_chr)[as.integer(cord):(as.integer(cord)+as.integer(dist)-1)])
    writeXStringSet(inside_seq, '~/STRANGE_BETWEEN_FLANCS/inside.fasta')
    inside_blast <- system('blastn -task megablast -query ~/STRANGE_BETWEEN_FLANCS/inside.fasta -db ~/TE_bases/TE_base_all_Justin_REPEATMODELER.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
    inside_blast_REZ <- as.data.frame(do.call('rbind', strsplit(inside_blast, '\t')))
    if (nrow(inside_blast_REZ) == 0){
      print('LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOK')
      print('Area between flancs has no simmilarity with our TE base')
      assign("strange_nondata", c(strange_nondata, N), envir = .GlobalEnv)
    } else {
      print('LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOK_rez')
      for (i in 1:nrow(inside_blast_REZ)){
        print(paste('Blast length = ', inside_blast_REZ$V4[i], 'Name: ', inside_blast_REZ$V2[i], 'parent_TE_name: ', GR_NEW[N]$TE_name))
        if (inside_blast_REZ$V2[i] == GR_NEW[N]$TE_name){
          assign("strange_me", c(strange_me, N), envir = .GlobalEnv)
          break
        } else {
          assign("strange_wit", c(strange_wit, N), envir = .GlobalEnv)
          assign("change_places", c(change_places, as.character(inside_blast_REZ$V2[i]), GR_NEW[N]$TE_name, 'NEXT:  '), envir = .GlobalEnv)
          break
        }
        print(paste('КООРДИНАТЫ БЛАСТА: ', inside_blast_REZ$V9[i], ' ', inside_blast_REZ$V10[i], ' score: ', inside_blast_REZ$V12[i], 'contig: ', inside_blast_REZ$V2[i]))
      }
    }
  }
}


genome_9 <- readDNAStringSet('~/9_POLISHED/9_correct_names.fasta')
New_one_code <- read.csv('~/RepeatMasker_UPD_NEWBASE_FINAL/RepeatMasker_160/elem_sorted/new_elem_sorted.csv', sep = '\t')
Len_table <- read.table('~/RepeatMasker_UPD_NEWBASE_FINAL/RepeatMasker_160/160_chromosome.fasta.out.length')
New_one_code$canon <- Len_table$V2[match(New_one_code$Element, Len_table$V1)]
New_one_code <- New_one_code[complete.cases(New_one_code$canon),]
New_one_code$DIF <- as.numeric(as.character(New_one_code$End.)) - as.numeric(as.character(New_one_code$Beg.))
New_one_code$RATIO <- New_one_code$DIF/New_one_code$canon

for (i in 1:nrow(New_one_code)){
  New_one_code$Correct_merged[i] <- ifelse(length(unique(unlist(strsplit(as.character(New_one_code$ID[i]), '/')))) == 1, TRUE, FALSE)
}
for (i in 1:nrow(New_one_code)){
  if (as.numeric(as.character(New_one_code$RATIO[i])) > 1.5 & length(unique(unlist(strsplit(as.character(New_one_code$ID[i]), '/')))) == 1)
  {
    New_one_code$Correct_merged[i] <- 'PENELOPE-LIKE' 
  }
}

New_one_code

New_one_code <- New_one_code[New_one_code$Correct_merged != FALSE,]
New_one_code <- New_one_code[New_one_code$RATIO > 0.5,]
nrow(New_one_code)

GR_NEW <- GRanges(seqnames = as.character(New_one_code$Query), ranges=IRanges(start=as.numeric(as.character(New_one_code$Beg.)), end=as.numeric(as.character(New_one_code$End.))), TE_name=as.character(New_one_code$Element), TE_family=as.character(New_one_code$Family), TYPE=as.character(New_one_code$Correct_merged), dif=as.numeric(as.character(New_one_code$DIF)), ratio=as.numeric(as.character(New_one_code$RATIO)), canon=as.numeric(as.character(New_one_code$canon)))
GR_NEW <- GR_NEW[seqnames(GR_NEW) %in% c('CM017604.2', 'CM017605.2', 'CM017606.2', 'CM017607.2', 'CM017608.2', 'CM017609.2')]

genome_160 <- readDNAStringSet('~/160_chromosome_assembly/160_chromosome.fasta')
names(genome_160) <- unlist(lapply(strsplit(names(genome_160), ' '), function(x)x[1]))
j <- 1
chromatin <- list()
for (i in 1:length(GR_NEW)){
  print(paste(i, ' ', chromatin_test(start(GR_NEW[i]), end(GR_NEW[i]), as.character(seqnames(GR_NEW[i])))))
  chromatin[j] <- chromatin_test(start(GR_NEW[i]), end(GR_NEW[i]), as.character(seqnames(GR_NEW[i])))
  j <- j+ 1
}
GR_NEW$chromatin <- chromatin
GR_NEW <- GR_NEW[GR_NEW$chromatin == 'Euchromatin',]
length(GR_NEW)
GR_NEW[GR_NEW$TYPE == 'PENELOPE-LIKE',]


TE_160 <- getSeq(genome_160, GR_NEW)
GR_NEW$sequence <- TE_160
GR_NEW$flanc_left <- IRanges( start = start(ranges(GR_NEW))-2000, end=start(ranges(GR_NEW)))
GR_NEW$flanc_rigth <- IRanges( start = end(ranges(GR_NEW)), end=end(ranges(GR_NEW))+2000)
GR_NEW$TE_length <- width(GR_NEW$sequence)
GR_NEW_left <- GRanges(seqnames = seqnames(GR_NEW), ranges=GR_NEW$flanc_left)
GR_NEW$left_fl_seq <- getSeq(genome_160, GR_NEW_left)
GR_NEW_rigth <- GRanges(seqnames = seqnames(GR_NEW), ranges=GR_NEW$flanc_rigth)
GR_NEW$rigth_fl_seq <- getSeq(genome_160, GR_NEW_rigth)

#######################################################
#Check in reference.
annotation <- rtracklayer::import('~/virilis_fly_base_1.07/dvir-all-r1.07.gtf')
#annotation@elementMetadata$type

GR_NEW$condition1_in_ref <- NA
GR_NEW$condition2_in_ref <- NA
GR_NEW$condition1_id<- NA
GR_NEW$condition2_id <- NA

GR_NEW$closest_gene_name <- NA
GR_NEW$closest_gene_id <- NA

for (N in 1:length(GR_NEW)){
  left_flanc <- GR_NEW$left_fl_seq[N]
  names(left_flanc) <- paste(N, ' left_flanc_seq')
  writeXStringSet(left_flanc, '~/Repeat_Masker_results/flanc_analysis/left_flanc.fasta')
  left_blast <- system('blastn -task megablast -query ~/Repeat_Masker_results/flanc_analysis/left_flanc.fasta -db ~/virilis_fly_base_1.07/dvir-all-chromosome-r1.07.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  left_df <- as.data.frame(do.call('rbind', strsplit(left_blast, '\t')))
  right_flanc <- GR_NEW$rigth_fl_seq[N]
  names(right_flanc) <- paste(N, ' rigth_flanc_seq')
  writeXStringSet(right_flanc, '~/Repeat_Masker_results/flanc_analysis/right_flanc.fasta')
  right_blast <- system('blastn -task megablast -query ~/Repeat_Masker_results/flanc_analysis/right_flanc.fasta -db ~/virilis_fly_base_1.07/dvir-all-chromosome-r1.07.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  right_df <- as.data.frame(do.call('rbind', strsplit(right_blast, '\t')))
  #################
  left_df <- left_df[as.numeric(as.character(left_df$V4)) > 1000,]
  right_df <- right_df[as.numeric(as.character(right_df$V4)) > 1000,]
  
  if (nrow(left_df) == 1 & nrow(right_df) == 1) {
    Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10))), end =max(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10)))))
    Fl_ref_2 <- GRanges(seqnames = as.character(right_df$V2), ranges = IRanges(start = min(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10))), end =max(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10)))))
    #Blast results for each flanc
    
    annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
    if (length(subjectHits(annotation_res_1)) != 0){
      GR_NEW$condition1_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_1)])
      GR_NEW$condition1_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_1)])
    }
    annotation_res_2 <- findOverlaps(Fl_ref_2, annotation, ignore.strand = T)
    if (length(subjectHits(annotation_res_2)) != 0){
      GR_NEW$condition2_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_2)])
      GR_NEW$condition2_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_2)])
    }
  }### 1-1 close
  else if (nrow(left_df) == 1){
    right_df <- right_df[as.character(right_df$V2) == as.character(left_df$V2),]
    count <- nrow(right_df)
    if (count == 0){
      next
    }
    if (as.numeric(as.character(left_df$V9)) < as.numeric(as.character(left_df$V10))){
      right_df <- right_df[as.numeric(as.character(right_df$V9)) < as.numeric(as.character(right_df$V10)),]
      if (count != 0 & nrow(right_df) == 0){
        next
      }
      i <- nrow(right_df) + 1
      dist_1_many <- vector()
      for (ii in 1:nrow(right_df)){
        dist_1_many <- c(dist_1_many, (abs(as.numeric(as.character(right_df$V9[ii])) - as.numeric(as.character(left_df$V10)))))
        if (abs(as.numeric(as.character(right_df$V9[ii])) - as.numeric(as.character(left_df$V10))) <= min(dist_1_many)){
          i <- ii
        }
      }
      Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10))), end =max(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10)))))
      Fl_ref_2 <- GRanges(seqnames = as.character(right_df$V2[i]), ranges = IRanges(start = min(as.numeric(as.character(right_df$V9[i])), as.numeric(as.character(right_df$V10[i]))), end =max(as.numeric(as.character(right_df$V9[i])), as.numeric(as.character(right_df$V10[i])))))
      annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_1)) != 0){
        GR_NEW$condition1_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_1)])
        GR_NEW$condition1_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_1)])
      }
      annotation_res_2 <- findOverlaps(Fl_ref_2, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_2)) != 0){
        GR_NEW$condition2_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_2)])
        GR_NEW$condition2_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_2)])
      }
    } 
    else {
      right_df <- right_df[as.numeric(as.character(right_df$V9)) > as.numeric(as.character(right_df$V10)),]
      if (count != 0 & nrow(right_df) == 0){
        next
      }
      dmin <- nrow(right_df) + 1
      dist_1_many_oposite <- vector()
      for (u in 1:nrow(right_df)){
        dist_1_many_oposite <- c(dist_1_many_oposite, abs(as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9[u]))) )
        if (abs(as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9[u]))) <= min(dist_1_many_oposite)){
          dmin <- u
        }
      }
      
      Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10))), end =max(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10)))))
      Fl_ref_2 <- GRanges(seqnames = as.character(right_df$V2[dmin]), ranges = IRanges(start = min(as.numeric(as.character(right_df$V9[dmin])), as.numeric(as.character(right_df$V10[dmin]))), end =max(as.numeric(as.character(right_df$V9[dmin])), as.numeric(as.character(right_df$V10[dmin])))))
      annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_1)) != 0){
        GR_NEW$condition1_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_1)])
        GR_NEW$condition1_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_1)])
      }
      annotation_res_2 <- findOverlaps(Fl_ref_2, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_2)) != 0){
        GR_NEW$condition2_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_2)])  
        GR_NEW$condition2_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_2)])
      }
    }
  }#end_left_1 right_many
  #right_1_left_many
  else if (nrow(right_df) == 1){
    left_df <- left_df[as.character(left_df$V2) == as.character(right_df$V2),]
    count_left_df <- nrow(left_df)
    if (count_left_df == 0){
      next
    }
    if (as.numeric(as.character(right_df$V9)) < as.numeric(as.character(right_df$V10))){
      left_df <- left_df[as.numeric(as.character(left_df$V9)) < as.numeric(as.character(left_df$V10)),]
      if (count_left_df != 0 & nrow(left_df) == 0){
        next
      }
      zmin <- nrow(left_df) + 1
      dist_many_1 <- vector()
      for (p in 1:nrow(left_df)){
        dist_many_1 <- c(dist_many_1, abs(as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10[p]))))
        if (abs(as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10[p]))) <= min(dist_many_1)){
          zmin <- p
        }
      }
      
      Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2[zmin]), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9[zmin])), as.numeric(as.character(left_df$V10[zmin]))), end =max(as.numeric(as.character(left_df$V9[zmin])), as.numeric(as.character(left_df$V10[zmin])))))
      Fl_ref_2 <- GRanges(seqnames = as.character(right_df$V2), ranges = IRanges(start = min(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10))), end =max(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10)))))
      annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_1)) != 0){
        GR_NEW$condition1_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_1)])
        GR_NEW$condition1_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_1)])
      }
      annotation_res_2 <- findOverlaps(Fl_ref_2, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_2)) != 0){
        GR_NEW$condition2_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_2)])  
        GR_NEW$condition2_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_2)])
      }
      
    } 
    else {
      left_df <- left_df[as.numeric(as.character(left_df$V9)) > as.numeric(as.character(left_df$V10)),]
      if (count_left_df != 0 & nrow(left_df) == 0){
        next
      }
      
      lmin <-nrow(left_df) + 1
      dist_many_1_Oposite <- vector()
      for (e in 1:nrow(left_df)){
        dist_many_1_Oposite <- c(dist_many_1_Oposite, abs(as.numeric(as.character(left_df$V10[e])) - as.numeric(as.character(right_df$V9))))
        if (abs(as.numeric(as.character(left_df$V10[e])) - as.numeric(as.character(right_df$V9))) <= min(dist_many_1_Oposite)){
          lmin <- e
        }
      }
      
      Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2[lmin]), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9[lmin])), as.numeric(as.character(left_df$V10[lmin]))), end =max(as.numeric(as.character(left_df$V9[lmin])), as.numeric(as.character(left_df$V10[lmin])))))
      Fl_ref_2 <- GRanges(seqnames = as.character(right_df$V2), ranges = IRanges(start = min(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10))), end =max(as.numeric(as.character(right_df$V9)), as.numeric(as.character(right_df$V10)))))
      annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_1)) != 0){
        GR_NEW$condition1_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_1)])
        GR_NEW$condition1_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_1)])
      }
      annotation_res_2 <- findOverlaps(Fl_ref_2, annotation, ignore.strand = T)
      if (length(subjectHits(annotation_res_2)) != 0){
        GR_NEW$condition2_in_ref[N] <- paste(annotation@elementMetadata$type[subjectHits(annotation_res_2)])      
        GR_NEW$condition2_id[N] <- paste(annotation@elementMetadata$gene_id[subjectHits(annotation_res_2)])
      }
      
    }
  } #end right 1  
  else{
    #######many_to_many
    function_rez <-  blast_reduce_upd(left_df, right_df)
    hard_dist <- function_rez[1]
    if (is.null(hard_dist)){
      next
    } else if (hard_dist == 'FUCK!'){
      next
    } else {
      
      #function_rez <- unlist(strsplit(function_rez, '   '))
      #hard_dist <- as.numeric(as.character(function_rez[1]))
      #coordinates_for_future_blast <-as.numeric(as.character(function_rez[2]) )
      #chromosome <- as.character(function_rez[3])
      GR_NEW$condition1_in_ref[N] <- 'Many_hits'
      GR_NEW$condition2_in_ref[N] <- 'Many_hits'
      GR_NEW$condition1_id[N] <- 'Many_hits'
      GR_NEW$condition2_id[N] <- 'Many_hits'
    }
  }
  ###loop end 
}

table(GR_NEW$condition1_in_ref)
table(GR_NEW$condition2_in_ref)

#Is it possible to TE be inserted in one gene? between genes?
WITHOUT_NA <- GR_NEW[!is.na(GR_NEW$condition1_id) & !is.na(GR_NEW$condition2_id),]
#here we have also Many_hits
length(WITHOUT_NA[WITHOUT_NA$condition2_id == 'Many_hits' | WITHOUT_NA$condition1_id == 'Many_hits',])
length(WITHOUT_NA[WITHOUT_NA$condition1_id != 'Many_hits' & WITHOUT_NA$condition2_id != 'Many_hits',])


#TE_inserted_inside_genes 177.
TE_INSIDE_GENE <- WITHOUT_NA[WITHOUT_NA$condition1_id == WITHOUT_NA$condition2_id,]
TE_INSIDE_GENE <- TE_INSIDE_GENE[TE_INSIDE_GENE$condition1_id != 'Many_hits',]
TE_INSIDE_GENE <- TE_INSIDE_GENE[TE_INSIDE_GENE$ratio > 0.9,]

data_TE_INSIDE_GENE <- as.data.frame(TE_INSIDE_GENE)
#Интересно сгруппировать по генам. Интересно, есть ли какой-то паттерн, какой тип траснпозона любит встраиваться в какие гены.


here <- aggregate(TE_name ~ condition1_id, TE_INSIDE_GENE, length)
harmful_genes <- here[here$TE_name != 1,]$condition1_id

for (i in 1:length(harmful_genes)){
  print(harmful_genes[i])
  print(paste(TE_INSIDE_GENE[TE_INSIDE_GENE$condition1_id == harmful_genes[i],]$TE_name, ' ', seqnames(TE_INSIDE_GENE[TE_INSIDE_GENE$condition1_id == harmful_genes[i],]), ' ', ranges(TE_INSIDE_GENE[TE_INSIDE_GENE$condition1_id == harmful_genes[i],])))
}


length(GR_NEW[!is.na(GR_NEW$condition1_id) & !is.na(GR_NEW$condition2_id),])
length(GR_NEW[is.na(GR_NEW$condition1_id) & is.na(GR_NEW$condition2_id),])

##########################################################################################################

empty <- vector()
same <- vector()
strange <- vector()
no_flanc <- vector()
inversia <- vector()

strange_nonsence <- vector()
strange_nondata <- vector()
strange_me <- vector()
strange_wit <- vector()
change_places <- vector()

INSERTION_REZ <- data.frame(I_Number = 1:length(GR_NEW), Chr=seqnames(GR_NEW), TE_name=GR_NEW$TE_name, TE_family=GR_NEW$TE_family, dif_initial=GR_NEW$dif, ratio_initial=GR_NEW$ratio, canon=GR_NEW$canon, between_flancs=NA, contig=NA, start_for_blast=NA, type=GR_NEW$TYPE, rez=NA)

for (N in 1:length(GR_NEW)){
  left_flanc <- GR_NEW$left_fl_seq[N]
  names(left_flanc) <- paste(N, ' left_flanc_seq')
  writeXStringSet(left_flanc, '~/Repeat_Masker_results/flanc_analysis/left_flanc.fasta')
  left_blast <- system('blastn -task megablast -query ~/Repeat_Masker_results/flanc_analysis/left_flanc.fasta -db ~/9_POLISHED/9_correct_names.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  left_df <- as.data.frame(do.call('rbind', strsplit(left_blast, '\t')))
  right_flanc <- GR_NEW$rigth_fl_seq[N]
  names(right_flanc) <- paste(N, ' rigth_flanc_seq')
  writeXStringSet(right_flanc, '~/Repeat_Masker_results/flanc_analysis/right_flanc.fasta')
  right_blast <- system('blastn -task megablast -query ~/Repeat_Masker_results/flanc_analysis/right_flanc.fasta -db ~/9_POLISHED/9_correct_names.fasta -outfmt 6 -num_threads 128 -evalue 1e-15', intern = T)
  right_df <- as.data.frame(do.call('rbind', strsplit(right_blast, '\t')))
  #################
  left_df <- left_df[as.numeric(as.character(left_df$V4)) > 500,]
  right_df <- right_df[as.numeric(as.character(right_df$V4)) > 500,]
  
  if (nrow(left_df) == 0 | nrow(right_df) == 0){
    print(paste(N, ' there is NO flanc for this TE.'));
    INSERTION_REZ[N,12] <- 'No_flancs_00'
    no_flanc <- c(no_flanc, N)
    next
  } else if (nrow(left_df) == 1 & nrow(right_df) == 1) {
    #####both 1
    if (as.character(left_df$V2) != as.character(right_df$V2)){
      print(paste('1-1 BLAST rez, but contigs are different. Insertion number:', N));
      no_flanc <- c(no_flanc, N)
      INSERTION_REZ[N,12] <- '11_but_names_No_flancs'
      next
    } else if (    (as.numeric(as.character(left_df$V10)) - as.numeric(as.character(left_df$V9)))*(as.numeric(as.character(right_df$V10)) - as.numeric(as.character(right_df$V9)))     < 0   ){
      print(paste('Problems in 1-1 flanc orientation. INVERSIA. Insertion Number: ', N))
      inversia <- c(inversia, N)
      INSERTION_REZ[N,12] <- '11_but_inversia'
      next
    } else if (   as.numeric(as.character(left_df$V9)) < as.numeric(as.character(left_df$V10))     ){
      dist <- as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10))
      INSERTION_REZ[N,8] <- dist
      INSERTION_REZ[N,9] <- as.character(left_df$V2)
      INSERTION_REZ[N,10] <- as.numeric(as.character(left_df$V10))
      INSERTION_REZ[N,12] <- '11_FORWARD'
      
      if (dist >= GR_NEW$TE_length[N] - 50 & dist <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist <= 20){
        print(paste('Empty space 1- 1. Insertion number: ', N, 'Dist = ', dist))
        empty <- c(empty, N)
      } else {
        print(paste('Smth strange. Insertion number: ', N, 'canon:', GR_NEW$canon[N]))
        strange_analisys(dist, as.numeric(as.character(left_df$V10)), as.character(left_df$V2), N)
        strange <- c(strange, N)
      }
    } #end of norm orientation in 1-1 way
    #1-1 but opposite orientation
    else {
      
      dist_oposite <- as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9))
      
      INSERTION_REZ[N,8] <- dist_oposite
      INSERTION_REZ[N,9] <- as.character(left_df$V2)
      INSERTION_REZ[N,10] <- as.numeric(as.character(right_df$V9))
      INSERTION_REZ[N,12] <- '11_REVERSE'
      
      if (dist_oposite >= GR_NEW$TE_length[N] - 50 & dist_oposite <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist_oposite <= 20){
        print(paste('Empty space 1 -1 Op. Insertion number: ', N, 'Dist = ', dist_oposite))
        empty <- c(empty, N)
      } else{
        print(paste('Smth strange. Insertion number: ', N, 'canon:', GR_NEW$canon[N]))
        strange_analisys(dist_oposite, as.numeric(as.character(right_df$V9)), as.character(left_df$V2), N)
        strange <- c(strange, N)
      }
    }#oposite orientation end
  }### 1-1 close
  #left 1, right many
  else if (nrow(left_df) == 1){
    right_df <- right_df[as.character(right_df$V2) == as.character(left_df$V2),]
    count <- nrow(right_df)
    if (count == 0){
      print(paste('PROBLEM WITH NAMES OF CONTIGS. Insertion number: ', N))
      no_flanc <- c(no_flanc, N)
      INSERTION_REZ[N,12] <- 'No_flans_1many_NAMES_filtering'
      next
    }
    if (as.numeric(as.character(left_df$V9)) < as.numeric(as.character(left_df$V10))){
      right_df <- right_df[as.numeric(as.character(right_df$V9)) < as.numeric(as.character(right_df$V10)),]
      if (count != 0 & nrow(right_df) == 0){
        print(paste('INVERSIA. Insertion number: ', N))
        INSERTION_REZ[N,12] <- 'Insersia_1many'
        inversia <- c(inversia, N)
        next
      }
      dist_1_many <- min(abs(as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10))))
      
      INSERTION_REZ[N,8] <- dist_1_many
      INSERTION_REZ[N,9] <- as.character(left_df$V2)
      INSERTION_REZ[N,10] <- as.numeric(as.character(left_df$V10))
      INSERTION_REZ[N,12] <- 'Forward_1many'
      
      if (dist_1_many >= GR_NEW$TE_length[N] - 50 & dist_1_many <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist_1_many <= 20){
        print(paste('Empty space 1 - many. Insertion number: ', N, 'Dist = ', dist_1_many))
        empty <- c(empty, N)
      } else{
        print(paste('Smth strange. Insertion number: ', N, 'canon:', GR_NEW$canon[N]))
        strange_analisys(dist_1_many, as.numeric(as.character(left_df$V10)), as.character(left_df$V2), N)
        strange <- c(strange, N)
      }
    } #end norm orientation
    else {
      right_df <- right_df[as.numeric(as.character(right_df$V9)) > as.numeric(as.character(right_df$V10)),]
      if (count != 0 & nrow(right_df) == 0){
        print(paste('INVERSIA. Insertion number: ', N))
        inversia <- c(inversia, N)
        INSERTION_REZ[N,12] <- 'Insersia_1many'
        next
      }
      
      #################
      
      dmin <- nrow(right_df) + 1
      dist_1_many_oposite <- vector()
      for (u in 1:nrow(right_df)){
        dist_1_many_oposite <- c(dist_1_many_oposite, abs(as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9[u]))) )
        if (abs(as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9[u]))) <= min(dist_1_many_oposite)){
          dmin <- u
        }
      }
      dist_1_many_oposite <- min(dist_1_many_oposite)
      INSERTION_REZ[N,8] <- dist_1_many_oposite
      INSERTION_REZ[N,9] <- as.character(right_df$V2[dmin])
      INSERTION_REZ[N,10] <- as.numeric(as.character(right_df$V9[dmin]))
      INSERTION_REZ[N,12] <- 'REVERSE_1many'
      
      ##################
      if (dist_1_many_oposite >= GR_NEW$TE_length[N] - 50 & dist_1_many_oposite <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist_1_many_oposite <= 20){
        print(paste('Empty space 1 - many OP. Insertion number: ', N, 'Dist = ', dist_1_many_oposite))
        empty <- c(empty, N)
      } else{
        print(paste('Smth strange. Insertion number: ', N, 'canon: ', GR_NEW$canon[N]))
        strange_analisys(dist_1_many_oposite, as.numeric(as.character(right_df$V9[dmin])), as.character(right_df$V2[dmin]), N)
        strange <- c(strange, N)
      }
      
    }
  }#end_left_1 right_many
  #right_1_left_many
  else if (nrow(right_df) == 1){
    left_df <- left_df[as.character(left_df$V2) == as.character(right_df$V2),]
    count_left_df <- nrow(left_df)
    if (count_left_df == 0){
      print(paste('PROBLEM WITH NAMES. Insertion number: ', N))
      INSERTION_REZ[N,12] <- 'No_flancs_names_many_1'
      no_flanc <- c(no_flanc, N)
      next
    }
    if (as.numeric(as.character(right_df$V9)) < as.numeric(as.character(right_df$V10))){
      left_df <- left_df[as.numeric(as.character(left_df$V9)) < as.numeric(as.character(left_df$V10)),]
      if (count_left_df != 0 & nrow(left_df) == 0){
        print(paste('INVERSIA. INSERTION NUMBER: ', N))
        inversia <- c(inversia, N)
        INSERTION_REZ[N,12] <-'Inversia_many1'
        next
      }
      ####################################################################################################
      zmin <- nrow(left_df) + 1
      dist_many_1 <- vector()
      for (p in 1:nrow(left_df)){
        dist_many_1 <- c(dist_many_1, abs(as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10[p]))))
        if (abs(as.numeric(as.character(right_df$V9)) - as.numeric(as.character(left_df$V10[p]))) <= min(dist_many_1)){
          zmin <- p
        }
      }
      dist_many_1 <- min(dist_many_1)
      INSERTION_REZ[N,8] <- dist_many_1
      INSERTION_REZ[N,9] <- as.character(left_df$V2[zmin])
      INSERTION_REZ[N,10] <- as.numeric(as.character(left_df$V10[zmin]))
      INSERTION_REZ[N,12] <- 'Forward_many_1'
      
      ######################################################################################################
      if (dist_many_1 >= GR_NEW$TE_length[N] - 50 & dist_many_1 <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist_many_1 <= 20){
        print(paste('Empty space. many - 1. Insertion number: ', N, 'Dist = ', dist_many_1))
        empty <- c(empty, N)
      } else {
        print(paste('Smth strange. Insertion number: ', N, 'canon: ', GR_NEW$canon[N]))
        strange_analisys(dist_many_1, as.numeric(as.character(left_df$V10[zmin])), as.character(left_df$V2[zmin]), N)
        strange <- c(strange, N)
      }
    } 
    else {
      left_df <- left_df[as.numeric(as.character(left_df$V9)) > as.numeric(as.character(left_df$V10)),]
      if (count_left_df != 0 & nrow(left_df) == 0){
        print(paste('INVERSIA. INSERTION NUMBER: ', N))
        inversia <- c(inversia, N)
        INSERTION_REZ[N,12] <- 'Insersia_many_1'
        next
      }
      ####################################
      dist_many_1_Oposite <- min(abs(as.numeric(as.character(left_df$V10)) - as.numeric(as.character(right_df$V9))))
      INSERTION_REZ[N,8] <- dist_many_1_Oposite
      INSERTION_REZ[N,9] <- as.character(right_df$V2)
      INSERTION_REZ[N,10] <- as.numeric(as.character(right_df$V9))
      INSERTION_REZ[N,12] <- 'REVERSE_many_1'
      
      ###########################################
      
      if (dist_many_1_Oposite >= GR_NEW$TE_length[N] - 50 & dist_many_1_Oposite <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (dist_many_1_Oposite <= 20){
        print(paste('Empty space many - 1 OP. Insertion number: ', N, 'Dist = ', dist_many_1_Oposite))
        empty <- c(empty, N)
      } else {
        print(paste('SMth strange. Insertion number: ', N, 'canon: ', GR_NEW$canon[N]))
        strange_analisys(dist_many_1_Oposite, as.numeric(as.character(right_df$V9)), as.character(right_df$V2), N)
        strange <- c(strange, N)
      }
    }
  } #end right 1  
  else{
    #######many_to_many
    print(paste('I Use blast_reduce_function. Many to many hits: Insertion number: ', N))
    function_rez <-  blast_reduce_upd(left_df, right_df)
    ##################Here is uncorrect output 
    hard_dist <- function_rez[1]
    
    if (is.null(hard_dist)){
      print(paste(N, ' there is NO flanc for this TE after blast filtering.'));
      INSERTION_REZ[N,12] <- 'NO_flancs_after_many_to_many_filtering'
      no_flanc <- c(no_flanc, N)
      next
    } else if (hard_dist == 'FUCK!'){
      print(paste('INVERSIA WITHIN FLANCS!!!!! NUMBER: ', N))
      inversia <- c(inversia, N)
      INSERTION_REZ[N,12] <- 'Inversion_betweeen_flancs'
      next
    } else {
      
      function_rez <- unlist(strsplit(function_rez, '   '))
      hard_dist <- as.numeric(as.character(function_rez[1]))
      coordinates_for_future_blast <-as.numeric(as.character(function_rez[2]) )
      chromosome <- as.character(function_rez[3])
      
      INSERTION_REZ[N,8] <- hard_dist
      INSERTION_REZ[N,9] <- as.character(chromosome)
      INSERTION_REZ[N,10] <- as.numeric(as.character(coordinates_for_future_blast))
      INSERTION_REZ[N,12] <- 'Many_to_many_rez'
      
      
      if (hard_dist >= GR_NEW$TE_length[N] - 50 & hard_dist <= GR_NEW$TE_length[N] + 50){
        print(paste('Insertion without localization changes. Insertion number: ', N))
        same <- c(same, N)
      } else if (hard_dist <= 20){
        print(paste('Empty space many-many. Insertion number: ', N, 'Dist = ', hard_dist))
        empty <- c(empty, N)
      } else{
        print(paste('Smth strange. Insertion number: ', N, 'Dist =', hard_dist, ' canon: ', GR_NEW$canon[N]))
        strange_analisys(hard_dist, coordinates_for_future_blast, chromosome, N)
        strange <- c(strange, N)
      }
    }
  }
  ###loop end 
}

View(INSERTION_REZ)
table(INSERTION_REZ$rez)
length(strange_nondata)
change_places

library(ggplot2)
table(INSERTION_REZ$rez)
sum(table(INSERTION_REZ$rez))
length(empty) + length(same) + length(strange) + length(no_flanc) + length(inversia)
length(inversia)


count_df_160 <- as.data.frame(table(t(INSERTION_REZ$rez)))
count_df_160

ggplot(data = count_df_160, aes(y = Freq, x = factor(Var1))) + geom_boxplot() +   theme(axis.text.x = element_text(angle = 15 , hjust = 1, size = 12))


vec_1 <- c('empty', 'same', 'strange', 'no_flanc', 'inversia')
vec_2 <- c(length(empty), length(same), length(strange), length(no_flanc), length(inversia))

count_2_df_160 <- data.frame(x = vec_1, y = vec_2)
count_2_df_160                           
ggplot(data = count_2_df_160, aes(x = x, y = y)) + geom_boxplot() +   theme(axis.text.x = element_text(angle = 15 , hjust = 1, size = 12)) + xlab('According to distance between flancs after blast') + ylab('count')



Fl_ref_1 <- GRanges(seqnames = as.character(left_df$V2), ranges = IRanges(start = min(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10))), end =max(as.numeric(as.character(left_df$V9)), as.numeric(as.character(left_df$V10)))))
Fl_ref_1
annotation_res_1 <- findOverlaps(Fl_ref_1, annotation, ignore.strand = T)
length(subjectHits(annotation_res_1))
GR_NEW$condition1_in_ref <- paste(annotation@elementMetadata$type[annotation_res_1@subjectHits])


annotation@elementMetadata$gene_id

View(as.data.frame(annotation))


###################################################################################################



ggplot(data = INSERTION_REZ, aes(x = TE_name, y = between_flancs, fill = TE_family)) + 
  geom_boxplot() + 
  theme_bw(base_size = 16) + 
  labs(y = "Distance between flancs", x = "TE name", title = "Distance (TE name)", fill = "TE family") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 15 , hjust = 1, size = 12))


#Все по всем, полторазмерные.

library(ggplot2)
ggplot(data = INSERTION_REZ[INSERTION_REZ$ratio_initial > 0.9,], aes(x = TE_name, y = between_flancs, fill = TE_family)) + 
  geom_boxplot() + 
  theme_bw(base_size = 16) + 
  labs(y = "Distance between flancs", x = "TE family", title = "Distance (TE family)", fill = "TE family") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 15 , hjust = 1, size = 12))


INSERTION_REZ[INSERTION_REZ$TE_family == 'LINE/Penelope' & INSERTION_REZ$ratio_initial >0.9,]$between_flancs

#pen <- c(30, 35, 37, 40, 143, 183, 185, 187, 190, 250,355, 435, 440, 466, 470)
#table(GR_NEW[pen]$TE_name)

ggplot(data = INSERTION_REZ[INSERTION_REZ$ratio_initial > 0.9 & INSERTION_REZ$between_flancs <2000,], aes(x = TE_name, y = between_flancs, fill = TE_family)) + 
  geom_boxplot() + 
  theme_bw(base_size = 16) + 
  labs(y = "Distance between flancs", x = "TE name", title = "Distance (TE name)", fill = "TE family") +
  scale_fill_brewer(palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 15 , hjust = 1, size = 12))



INSERTION_REZ[INSERTION_REZ$ratio_initial > 0.9 & INSERTION_REZ$TE_name == 'dvir-934_full',]$between_flancs

INSERTION_REZ[INSERTION_REZ$ratio_initial > 0.9 & INSERTION_REZ$TE_name == 'Helena-2-2014_full',]$between_flancs


length(empty)
GR_empty <- GR_NEW[empty]
length(GR_empty[GR_empty$ratio > 0.9])/length(empty)
table(GR_NEW[empty]$TE_name)
length(same)
ii <- GR_NEW[no_flanc]
length(ii)
ii <- ii[ii$TE_name == 'Penelope_full',]
length(ii)
ii <- ii[ii$ratio > 0.9,]
GR_same <- GR_NEW[same]
length(GR_same[GR_same$ratio > 0.9])/length(same)
GR_same


INS
length(strange)
GR_strange <- GR_NEW[strange]

unique(GR_NEW[strange]$TE_name)

table(GR_NEW[strange]$TE_name)



length(GR_strange[GR_strange$ratio > 0.9])/length(strange)
table(GR_NEW[strange]$TE_name)

length(no_flanc)
length(inversia)


length(empty) + length(same) + length(no_flanc) + length(inversia) + length(strange)
vector <- c(empty, same, no_flanc, inversia, strange)

length(same)


##############check inversion for contigs. Distance.
###################### optimize длина фланков в остальных случаях.
###################### во фланках могут быть другие транспозонам. Эти могут попадать в ноу фланкс.
GR_NEW[185]
View(INSERTION_REZ)
i <- GR_NEW[GR_NEW$TE_name == 'Penelope_full',]
length(i[i$ratio > 0.9,])



base_to_blast <- readDNAStringSet('~/TE_bases/TE_base_all_Justin_REPEATMODELER.fasta')
base_to_blast['rnd-6_family-6033#RC/Helitron']
