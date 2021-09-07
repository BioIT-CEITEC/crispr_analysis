library(data.table)
library(Biostrings)
library(parallel)


final_blast_counts <- function(insert_files,alignment_counts,all_fasta_reads,counts,output_dir,cores,sample_name){
  # sample_name <- strsplit(alignment_counts, "trimmed_20bpInserts")[[1]][1]
  
  # CONTINUE HERE - FILTERING ALIGNMENT AND SECOND ROUND OPTIMAL ALIGNMENT PROCESSING
  # load unique inserts
  unique_inserts <- fread(insert_files) 
  names(unique_inserts) <- c("gene_id","UID","seq")
  unique_inserts[,merge_id := paste(gene_id,UID, sep="|")] # 
  
  
  ############################################################################
  #table with the list of hits from query (NGS read sequences) to the target 
  #(list of inserts) created by gast 
  a <- Sys.time()
  blast_counts <- fread(alignment_counts, header = F) 
  
  #filtering of only best alignemnt per each hit
  duplicated_counts_tab <- blast_counts[duplicated(V1)]
  unique_counts_tab <- blast_counts[!duplicated(V1)]    
  
  
  # PREPISANIE NASLEDUJUCICH RIADKOV
  # blast_counts[,to_remove := duplicated(V1)]
  # blast_counts[to_remove,to_remove := V3 != max(V3),by = V1]
  # blast_counts[to_remove,to_remove := V4 != max(V4),by = V1]
  # blast_counts[to_remove,to_remove := V5 != min(V5),by = V1]
  # blast_counts <- unique(blast_counts[!to_remove,],by = "V1")
  
  unique_duplicated_counts_tab <- duplicated_counts_tab[,.SD[which(V3 == max(V3) & V4 == max(V4) & V5 == min(V5))],by = V1]
  unique_duplicated_counts_tab <- unique_duplicated_counts_tab[!duplicated(unique_duplicated_counts_tab$V1)] 
  best_alignments_with_duplicates <- rbind(unique_counts_tab,unique_duplicated_counts_tab) 
  setorder(best_alignments_with_duplicates,V1) 
  
  #keep only the hits with alignment longer then 17 nt
  best_alignments_with_duplicates_17aln <- best_alignments_with_duplicates[V4>=17] 
  best_alignments_with_duplicates_17aln <- best_alignments_with_duplicates_17aln[!V1 %in% best_alignments_with_duplicates_17aln[V4==17 & V5==1]$V1] 
  
  
  ############################################################################
  #looking for reads, that are asigned exactly to just one insert sgRNA
  #this table contains same number of hits as best_alignments_with_duplicates_17aln
  unique_tab_blast <- data.table(table(best_alignments_with_duplicates_17aln$V1))[N==1]  
  setnames(unique_tab_blast,"V1","target_name") 
  
  
  ############################################################################
  #sequences that has unique hits are concatenated with their counts from the *.counts files
  #load original fasta file which can be joined with *counts file

  insert_reads <- readLines(all_fasta_reads)  
  insert_reads_tab <- data.table(name = insert_reads[seq(1,length(insert_reads),2)],seq = insert_reads[seq(2,length(insert_reads),2)])
  insert_reads_tab$name <- gsub(">","",insert_reads_tab$name) 
  
  counts_file <- fread(counts)       
  setnames(counts_file,"V1","reads") 
  setnames(counts_file,"V2","N")     
  merged_inserts_reads_seqs_counts <- merge(insert_reads_tab,counts_file,by.x="seq",by.y = "reads")
  merged_inserts_reads_seqs_counts$name <- as.integer(merged_inserts_reads_seqs_counts$name)
  final_tab_for_counting <- best_alignments_with_duplicates_17aln[V1 %in% unique_tab_blast$target_name]
  
  #tab which contains any read seqeunces that occured, with its most probable insert they belong to
  final_tab_for_counting <- merge(merged_inserts_reads_seqs_counts,final_tab_for_counting,by.x="name",by.y="V1")
  #unique_blast_count_tab <- data.table(table(final_tab_for_counting$V2))
  #setnames(unique_blast_count_tab,"V1","merge_id")
  #setorder(unique_blast_count_tab,-N)
  
  
  ############################################################################
  #searching for reads, that had no match with any inserts
  #list of reads that were not mapped in any way to any insert target
  #this number together with unique_vsearch_count and duplication_tab gives the number of all reads
  #vector of reads, that are in final_tab_for_counting and added the reads in duplicates, that should be resolved
  
  
  not_mapped_by_blast_reads <- insert_reads_tab[!name %in% final_tab_for_counting$name]
  not_mapped_by_blast_reads_counts <- data.table(table(not_mapped_by_blast_reads$seq))
  setnames(not_mapped_by_blast_reads_counts,"V1","seq")
  setorder(not_mapped_by_blast_reads_counts,-N)
  
  
  ############################################################################
  #OPTIMAL ALIGNMENT
  # the rest of reads, that were not mapped at least with 17 bp length alignments
  # is forwarded for the analysis  by optimal alignment
  
  not_mapped_by_blast_with_counts <- merge(not_mapped_by_blast_reads_counts[,1],merged_inserts_reads_seqs_counts,by="seq")
  sm <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
  diag(sm) <- 1
  sm[,"N"] <- 1
  sm["N",] <- 1
  
  
  #list of the best alignment , when aligning the sequences of non mapped reads to insert database, means that te first number is number of insert read
  # which fits for the first unmapped sequence in best way
  #this need to be filtered according the scores and laignment lengths
  #only the best hits from each read are taken
  
  optimal_alignment_of_not_mapped <- unlist(mclapply(not_mapped_by_blast_with_counts[N>=10]$seq,function(target) {
    which.max(pairwiseAlignment(pattern = unique_inserts$seq, subject = target,substitutionMatrix = sm,gapOpening = -0.1, gapExtension = -1,type = "local",scoreOnly = T))
  },mc.cores = cores))
  
  res <- pairwiseAlignment(pattern = unique_inserts$seq[optimal_alignment_of_not_mapped],subject = not_mapped_by_blast_with_counts[N>=10]$seq,substitutionMatrix = sm,gapOpening = -0.1, 
  gapExtension = -1,type = "local" )
  
  #res@pattern@range@width
  #res@score
  #get only those alignment whose score and alignment length are at least 18
  #keep the new alignments, these numbers follow the order as in unique inserts table
  #unique_inserts[optimal_alignment_of_not_mapped[res@pattern@range@width >= 17 & res@score >= 14]]
  #keep the same order from not mapped reads
  #not_mapped_by_blast_with_counts[N>=10][res@pattern@range@width >= 17 & res@score >= 14]
  
  optimal_aln_resolved_not_mapped <- cbind(unique_inserts[optimal_alignment_of_not_mapped[res@pattern@range@width >= 17 & res@score >= 14],c(1,2,4)],not_mapped_by_blast_with_counts[N>=10][res@pattern@range@width >= 17 & res@score >= 14])
  summed_optimal_aln_resolved_not_mapped <- optimal_aln_resolved_not_mapped[,list(N = sum(N) ), by = c("merge_id")]
  setnames(summed_optimal_aln_resolved_not_mapped,"merge_id","V2")
  
  final_blast_tab_counts <- rbind(final_tab_for_counting[,c(4,3)],summed_optimal_aln_resolved_not_mapped)
  final_blast_tab_counts <- final_blast_tab_counts[,list(N = sum(N) ), by = "V2"]
  setnames(final_blast_tab_counts,"V2","merge_id")
  
  #list of the inserts, and counts of reads assigned to them, except those reads that were targeted by the same read (the unique targets)
  insert_counts <- merge(unique_inserts,final_blast_tab_counts,by="merge_id")
  insert_counts$merge_id <- NULL
  setorder(insert_counts,-N)
  sum(insert_counts$N)
  fwrite(insert_counts, paste0(output_dir,sample_name,"unique_inserts_tab_counts.tsv"), sep = "\t", row.names = F, quote = F)
  
  
  ############################################################################
  #final rest of unmapped reads
  
  
  final_not_mapped_by_blast_reads <- insert_reads_tab[!seq %in% c(final_tab_for_counting$seq,optimal_aln_resolved_not_mapped$seq)]
  
  final_not_mapped_by_blast_reads_counts <- merge(final_not_mapped_by_blast_reads[,2],merged_inserts_reads_seqs_counts,by="seq")
  
  final_not_mapped_by_blast_reads_counts400 <- final_not_mapped_by_blast_reads_counts[N>=400]
  if (dim(final_not_mapped_by_blast_reads_counts400)[1] >=1){
    setorder(final_not_mapped_by_blast_reads_counts400,-N)
    fa = character(2 * nrow(final_not_mapped_by_blast_reads_counts400))
    fa[c(TRUE, FALSE)] = sprintf(">%s", seq(1,length(final_not_mapped_by_blast_reads_counts400$seq)),1)
    fa[c(FALSE, TRUE)] = final_not_mapped_by_blast_reads_counts400$seq
    writeLines(fa, paste0(output_dir,sample_name,"unique_not_mapped_reads_400.tsv"))
  }
  setorder(final_not_mapped_by_blast_reads_counts,-N)
  fwrite(final_not_mapped_by_blast_reads_counts, paste0(output_dir,sample_name,"unique_not_mapped_reads.tsv"),sep = "\t",row.names = F,quote = F)
  sum(final_not_mapped_by_blast_reads_counts$N)
  
  
  
  ############################################################################
  #inserts that were not mapped
  not_hitted_inserts <- unique_inserts[!merge_id %in% paste(insert_counts$gene_id,insert_counts$UID, sep="|")]
  not_hitted_inserts$merge_id <- NULL
  fwrite(not_hitted_inserts,paste0(output_dir,sample_name,"unique_not_hitted_inserts.tsv"),sep = "\t",row.names = F,quote = F)
  
  
  ############################################################################
  #GRAPHS
  pdf(file = paste0(output_dir,sample_name,"graphs.pdf"))
  hist(final_not_mapped_by_blast_reads_counts$N,breaks = 200,ylim = c(0,1000),main = "non-insert sequences histogram",xlab = "read count of sequence")
  hist(insert_counts[N > 0,]$N,breaks = 300,main = "inserts histogram",xlab = "read count of inserts")
  hist(insert_counts[N > 0 & N < 1000,]$N,breaks = 300,main = "inserts histogram",xlab = "read count of inserts")
  hist(insert_counts[N > 0 & N < 200,]$N,breaks = 300,main = "inserts histogram",xlab = "read count of inserts")
  dev.off()
  
  
  ############################################################################
  #STATISTICS
  stat <- data.table(`number of reads matching any insert` = sum(insert_counts$N), 
                     `number of unmapped reads` = sum(final_not_mapped_by_blast_reads_counts$N), 
                     `original number of unique inserts` = length(unique_inserts$gene_id),
                     `number of targeted inserts` = length(insert_counts$gene_id),
                     `number of not targeted inserts` = length(unique_inserts$gene_id) - length(insert_counts$gene_id))
  fwrite(stat,paste0(output_dir,sample_name,"statistics.tsv"),sep = "\t",row.names = F,quote = F)
  
  Sys.time() - a
}

# debuging
#setwd("/mnt/ssd/ssd_3/temp/vasek/a156_obrdlik_crispr_filter_counts/")
#insert_files <- 'unique_with_revcomp_inserts.csv'
#alignment_counts <- 'Brunello_plasmid01_R1_trimmed_20bpInserts_unique_gapopen0_extend2_reward1_penalty__3.hsblastn.txt'
#all_fasta_reads <- 'Brunello_plasmid01_R1_trimmed_20bpInserts.seqs.fa'
#counts <- 'Brunello_plasmid01_R1_trimmed_20bpInserts.seqs.counts'
#output_dir <- 'results'

# running as Rscript
args = commandArgs(trailingOnly=TRUE)

insert_files <- args[1]
alignment_counts <- args[2]
all_fasta_reads <- args[3]
counts <- args[4]
output_dir <- args[5]
cores <- args[6]
sample <- args[7]

if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
print(insert_files)
print(alignment_counts)
print(all_fasta_reads)
print(counts)
print(output_dir)
print(cores)
print(sample)

final_blast_counts(insert_files,
                   alignment_counts,
                   all_fasta_reads,
                   counts,
                   output_dir,
                   cores,
                   sample)
