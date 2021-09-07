#Script for observe differential expressed sgRNAs and genes from crispr cas9 wide screen experiments
########## 4 inputs needed for the script
##### first input is the library of crispr sgRNAs in format:
# first column is gene which the sgRNA is related to
# second column is sgRNA identifer
# third column is the sequence of sgRNA insert
##### as an second input serves the count table in the frmat following mageck counts table
# first column sgRNA identifier
# second column gene which sgRNA is related to
# any number of samples, each sample is one column and contains counts per each sgRNA
##### the third input is desing matrix for the particular analysis in format:
# first column original name of sample
# second column the name of the sample in counts table
# third column is a condition (if the sample is control or treatment)
#fourth column is patient (from which batch the individual samples came from) - not mandatory
##### the fourth input is the path to the output folder, if the folder does not exist, it is created in the same directory, where the script is executed 

library(edgeR)
library(data.table)
library(ggplot2)
library(gplots)

# setwd("/mnt/ssd/ssd_1/snakemake/")
# args <- c("/mnt/ssd/ssd_3/references/general/CRISPR_Brunello/CRISPR_Brunello_mod.csv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/hsblastn_filter/all_samples_report.tsv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/DE_genes_edgeR/1668_d7_vs_DMSO_d7/1668_d7_vs_DMSO_d7.design_table.tsv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/DE_genes_edgeR/1668_d7_vs_DMSO_d7/",
#           20)

args = commandArgs(trailingOnly=TRUE)

insert_lib_input <-args[1]
counts_tab_input <- args[2]
design_input <- args[3]
output_dir <- args[4]
top_genes <- as.numeric(args[5])

if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

#####################################
insert_lib <- fread(insert_lib_input)
colnames(insert_lib) <- c("Gene","ID","Sequences")
insert_lib <- insert_lib[,c(2,3,1)]

samples <- read.table(counts_tab_input, header = T)
samples$Gene <- NULL
rownames(samples) <- samples$sgRNA
samples$sgRNA <- NULL

design_tab <- read.table(design_input, stringsAsFactors = F, sep="\t", header = T)
design_tab = as.data.frame(as.data.table(design_tab)[name %like% "^[0-9]", name:=paste0("X",name)])
conds <- factor(c(design_tab$condition), levels=c(unique(design_tab$condition)))
coldata <- as.data.frame(conds)
# coldata <- as.data.frame(t(t(conds)))
colnames(coldata) <- "condition"
rownames(coldata) <- design_tab$name # if we read from a sample sheet
coldata <- as.data.frame(coldata)

# Remove samples not in a design (to avoid errors later)
samplesNotInDesign <- colnames(samples)[!(colnames(samples) %in% design_tab$name)]
mrcounts <- samples[, colnames(samples) %in% design_tab$name]

# Drop levels of removed samples
coldata$condition <- droplevels(coldata$condition)

mrcounts <- mrcounts[,match(rownames(coldata), colnames(mrcounts))]

mrcounts_rownames <- copy(mrcounts)
mrcounts_rownames$ID <- row.names(mrcounts_rownames)
input_tab <-  as.data.frame(merge(insert_lib,mrcounts_rownames,by.x = "ID", by.y="ID"))
row.names(input_tab) <- input_tab$ID
#d<-DGEList(counts=mrcounts[order(row.names(mrcounts)),], group=coldata$condition, genes=insert_lib[order(insert_lib$ID),]) # edgeR DGE object
#new
# d <- DGEList(counts=input_tab[c(4,5,6,7)], group=coldata$condition, genes=input_tab[c(1,2,3)]) # edgeR DGE object
d <- DGEList(counts=input_tab[c(4:ncol(input_tab))], group=coldata$condition, genes=input_tab[c(1,2,3)]) # edgeR DGE object
d <- calcNormFactors(d) # Calculate normalization factors
cols <- as.numeric(d$samples$group)+2


pdf(file = file.path(paste0(output_dir,"graphs.pdf")))

par(mfrow=c(2,1))
barplot(colSums(d$counts), las=2, main="Counts per sample",
        col=cols, cex.names=0.5, cex.axis=0.8)
#legend("topright", legend=c("Control", "Treatment"), col=c(3,4), pch=15)
barplot(rowSums(d$counts), las=2, main="Counts per sgRNA",
        axisnames=FALSE, cex.axis=0.8)

logCPM<-cpm(d, log=TRUE, prior.count=1)
logCPMc<-removeBatchEffect(logCPM, coldata$patient)

cols2 <- d$samples$group
#MDS needs at least
if(ncol(samples) > 4 && nrow(design_tab) > 2){
plotMDS(logCPM, col=cols, main="MDS Plot: with batch effect")
#legend("topleft", legend=c("Control", "Drug"), col=c(3,4), pch=15)
plotMDS(logCPMc, col=cols, main="MDS Plot: removed batch effect")
#legend("topleft", legend=c("Inf#1", "Inf#2"), col=c(1,2), pch=15)
}
par(cex.main=0.8)
heatmap.2(cor(logCPM),trace = 'none',density.info = 'none',cexRow = 0.8,cexCol = 0.8,offsetRow = -0.2,offsetCol = -0.2,notecol="black",cellnote = round(cor(logCPM),3),notecex=0.7,
          main = "Correlation of samples with batch effect")

heatmap.2(cor(logCPMc),trace = 'none',density.info = 'none',cexRow = 0.8,cexCol = 0.8,offsetRow = -0.2,offsetCol = -0.2,notecol="black",cellnote = round(cor(logCPMc),3),notecex=0.7,
          main = "Correlation of samples with removed batch effect")
par(cex.main=1.0)
design <- model.matrix(~condition, data=coldata)

#if we have no replicates, we need to somehow fake dispersions --->>> use EXACT test
#Typical values for the common BCV (square-rootdispersion) for datasets arising from well-controlled experiments are 0.4 for human data,
#0.1 for data on genetically identical model organisms or 0.01 for technical replicates.

if(nrow(coldata) <= length(unique(coldata$condition))){
  bcv <- 0.2
  xglm <- estimateDisp(d, design)
  xglm$common.dispersion <- bcv^2
  xglm$trended.dispersion <- bcv^2
  xglm$tagwise.dispersion <- bcv^2
  fit <- glmFit(xglm, design)
  lrt <- glmLRT(fit, coef=ncol(fit$design))
} else{
  xglm <- estimateDisp(d, design,trend.method="loess")
  sqrt(xglm$common.disp)
  #plotBCV(xglm, main="BCV Plot")
  #use the function glmFit to fit the sgRNA-specific models and glmLRT to do the testing
  #between the drug treated and control samples. The top ranked sgRNAs are listed using the
  #topTags function
  fit <- glmFit(xglm, design)
  lrt <- glmLRT(fit, coef=ncol(fit$design))
}
# topTags(lrt) # print to the console, not needed in automatic processing
par(mfrow=c(1,1))
thresh <- 0.05
lfc <- 1
top4 <- topTags(lrt, n=Inf)
top4ids <- rownames(top4$table[abs(top4$table$logFC)>lfc & top4$table$FDR<thresh,])
plotSmear(lrt, de.tags=top4ids, pch=20, cex=0.6,
          main="logFC per single sgRNAs (lfc > 1, fdr < 0.05)")
abline(h=c(-1, 0, 1), col=c("dodgerblue","yellow","dodgerblue"), lty=2)

write.table(top4, paste0(output_dir,"sgRNAs_summary.tsv"), sep = "\t",row.names = F, quote = F)

genesymbols <- d$genes[,3]
#new
#genesymbols <- input_tab$Gene

genesymbollist <- list()
unq <- unique(genesymbols)
unq <- unq[!is.na(unq)]
for(i in unq) {
  sel <- genesymbols==i & !is.na(genesymbols)
  if(sum(sel)>=1)
    genesymbollist[[i]] <- which(sel)
}
camera.res <- camera(xglm, index=genesymbollist, design, contrast=2)

#add sgRNAs info
combine_tab <- as.data.table(d$genes)[,paste(ID,collapse = ";"),by = Gene]
camera.res <- merge(camera.res,combine_tab,by.x = "row.names", by.y="Gene")
colnames(camera.res)[which(names(camera.res) == "V1")] <- "sgRNAs"

#new part
camera.res <- as.data.table(camera.res)[,.(ID=unlist(strsplit(sgRNAs, ";"))),by=.(Row.names,NGenes,Direction,PValue,FDR)]
#add p-values per sgRNA info
#add logFC values per sgRNA info
camera.res <- merge(camera.res,as.data.table(lrt$table,keep.rownames = T), by.x = "ID", by.y = "rn")

#need to solve how to add counts based on ipnut
camera.res <- merge(camera.res,as.data.table(d$counts,keep.rownames = T), by.x = "ID", by.y = "rn")
camera.res <- merge(camera.res,as.data.table(logCPM,keep.rownames = T), by.x = "ID", by.y = "rn")

camera.res <- camera.res[, lapply(.SD, paste0, collapse=";"), by=.(Row.names,NGenes,Direction,PValue.x,FDR), .SDcols=!c("logCPM","LR")]

#end new part


#camera.res$P_values_per_sgRNA <- sapply(genesymbollist, function(x) paste0(lrt$table[x,]$PValue,collapse = ";"))

#add logFC values per sgRNA info
#camera.res$logFC_values_per_sgRNA <- sapply(genesymbollist, function(x) paste0(lrt$table[x,]$logFC,collapse = ";"))

#add raw counts of sgRNAs per sample


#for (sample in colnames(d$counts)){
#  camera.res[, paste0(sample,"_raw_counts")] <- sapply(genesymbollist, function(x) paste0(d$counts[,sample][x],collapse = ";"))
#  camera.res[, paste0(sample,"_logCPM_normalized_counts")] <- sapply(genesymbollist, function(x) paste0(logCPM[,sample][x],collapse = ";"))
#}
colnames(camera.res) <- c("Gene","NsgRNAs","Direction","PValue","FDR","sgRNAs","logFC_per_sgRNA","P_values_per_sgRNA",paste0(colnames(d$counts),"_raw_counts"),paste0(colnames(d$counts),"_logCPM_normalized_counts"))
setorder(camera.res,FDR)
fwrite(camera.res, paste0(output_dir,"gene_summary.tsv"), sep = "\t", row.names = F)

#row.names(d$genes[genesymbollist[["A1BG"]],]) 

for (gene in camera.res[FDR<=0.05, .SD[1:min(top_genes,.N),Gene]]){
  barcodeplot(lrt$table$logFC,index=genesymbollist[[gene]],
              main=paste("Barcodeplot for Gene",gene),
              labels=c("Negative logFC", "Positive logFC"),
              quantile=c(-0.5,0.5))
  # Data 
  b <-logCPMc[genesymbollist[[gene]],]
  b <- reshape2::melt(b)
  #Graph
  print(qplot( x=Var2 , y=value , data=b , geom=c("boxplot","jitter"),  fill=Var2, xlab = "Samples", ylab = "logCPM values for each sgRNA", main = paste0("logCPM values per sgRNA without batch effect for gene ", gene)))
}

dev.off()


