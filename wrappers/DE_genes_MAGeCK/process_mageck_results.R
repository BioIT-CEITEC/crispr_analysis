library(data.table)

setwd("/mnt/ssd/ssd_1/sequia/220307__crispr_analysis__2641")
path = "DE_results/MAGeCK/NOTCH1_d21_vs_NOTCH1_d0"
gene_tab = "NOTCH1_d21_vs_NOTCH1_d0.gene_summary.txt"
sg_tab   = "NOTCH1_d21_vs_NOTCH1_d0.sgrna_summary.txt"
norm_tab = "NOTCH1_d21_vs_NOTCH1_d0.normalized.txt"
input_gene_file = paste(path,gene_tab,sep = "/")
input_sg_file = paste(path,sg_tab,sep = "/")
input_norm_file = paste(path,norm_tab,sep = "/")
out_neg_genes = paste(path,"NOTCH1_d21_vs_NOTCH1_d0.gene_summary.neg_rank.tsv",sep = "/")
out_pos_genes = paste(path,"NOTCH1_d21_vs_NOTCH1_d0.gene_summary.pos_rank.tsv",sep = "/")
out_new_genes = paste(path,"NOTCH1_d21_vs_NOTCH1_d0.gene_summary.tsv",sep = "/")
# top_genes

args <- commandArgs(trailingOnly = T)
input_gene_file = args[1]
input_sg_file = args[2]
input_norm_file = args[3]
out_neg_genes = args[4]
out_pos_genes = args[5]
out_new_genes = args[6]

norm_tab = fread(file = input_norm_file, header = T, sep = "\t")
sg_tab = fread(file = input_sg_file, header = T, sep = "\t")
gene_tab = fread(file = input_gene_file, header = T, sep = "\t")
gene_tab[,l2FC:=fifelse(`neg|fdr`<`pos|fdr`, `neg|lfc`, `pos|lfc`)]
gene_tab[,good_sgRNAs:=fifelse(`neg|fdr`<`pos|fdr`, `neg|goodsgrna`, `pos|goodsgrna`)]
gene_tab[,`RRA-score`:=fifelse(`neg|fdr`<`pos|fdr`, `neg|score`, `pos|score`)]
gene_tab[,`Pvalue`:=fifelse(`neg|fdr`<`pos|fdr`, `neg|p-value`, `pos|p-value`)]
gene_tab[,FDR:=fifelse(`neg|fdr`<`pos|fdr`, `neg|fdr`, `pos|fdr`)]
gene_tab[,Direction:=fifelse(l2FC<0, "Down", "Up")]

fwrite(gene_tab[,.(Gene=id,
                   total_sgRNAs=`num`,
                   good_sgRNAs,
                   Direction,
                   l2FC,
                   Pvalue,
                   FDR)][order(FDR)],
       out_new_genes, sep = "\t", quote = F, row.names = F, col.names = T)

fwrite(gene_tab[order(`neg|rank`)][,.(gene=`id`,
                                  total_sgRNAs=`num`,
                                  good_sgRNAs=`neg|goodsgrna`,
                                  l2FC=`neg|lfc`,
                                  RRA_score=`neg|score`,
                                  pvalue=`neg|p-value`,
                                  FDR=`neg|fdr`)], 
       out_neg_genes, sep = '\t', quote = F, row.names = F, col.names = T)

fwrite(gene_tab[order(`pos|rank`)][,.(gene=`id`,
                                 total_sgRNAs=`num`,
                                 good_sgRNAs=`pos|goodsgrna`,
                                 l2FC=`pos|lfc`,
                                 RRA_score=`pos|score`,
                                 pvalue=`pos|p-value`,
                                 FDR=`pos|fdr`)], 
       out_pos_genes, sep = '\t', quote = F, row.names = F, col.names = T)

