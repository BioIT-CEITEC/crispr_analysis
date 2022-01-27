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

args = commandArgs(trailingOnly=TRUE)

# setwd("/mnt/nfs/shared/CFBioinformatics/all_projects/")
# setwd("/mnt/ssd/ssd_1/snakemake/")
# args <- c("/mnt/ssd/ssd_3/references/general/CRISPR_Brunello/CRISPR_Brunello_mod.csv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/hsblastn_filter/all_samples_report.tsv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/DE_genes_edgeR/1464_d14_vs_1668_d14/1464_d14_vs_1668_d14.design_table.tsv",
#           "stage470_Michal_-_Haspin_CRISPR/CRISPR_general_analysis/DE_genes_edgeR/1464_d14_vs_1668_d14/",
#           "20",
#           "True")

insert_lib_input <-args[1]
counts_tab_input <- args[2]
design_input <- args[3]
output_dir <- args[4]
top_genes <- as.numeric(args[5])
paired = as.logical(toupper(args[6]))

if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

#####################################
insert_lib <- fread(insert_lib_input)[,.(ID=UID,Sequences=seq,Gene=gene_id)]

samples <- read.table(counts_tab_input, header = T)
samples$Gene <- NULL
rownames(samples) <- samples$sgRNA
samples$sgRNA <- NULL

design_tab <- fread(design_input, stringsAsFactors = F, sep="\t", header = T)
design_tab[name %like% "^[0-9]", name:=paste0("X",name)]
design_tab[,N:=.N,by=condition]
if(design_tab[N==1,.N]>0) {
  new_cols = as.data.table(samples)[, design_tab[N == 1, name], with=F]
  colnames(new_cols) = sub('rep1','rep2',colnames(new_cols))
  samples = cbind(samples, as.data.frame(new_cols))
}
design_tab = as.data.frame(rbind(design_tab[,.SD,.SDcols=!c("N")], 
                                 design_tab[N==1, .(sample=sub('rep1','rep2',sample), 
                                                    name=sub('rep1','rep2',name), 
                                                    condition=as.factor(condition), 
                                                    patient=sub('rep1','rep2',patient))]))
rownames(design_tab) <- design_tab$name

# Remove samples not in a design (to avoid errors later)
samplesNotInDesign <- colnames(samples)[!(colnames(samples) %in% design_tab$name)]
mrcounts <- samples[, colnames(samples) %in% design_tab$name]
mrcounts <- mrcounts[,match(rownames(design_tab), colnames(mrcounts))]
mrcounts$ID <- row.names(mrcounts)
input_tab <-  as.data.frame(merge(insert_lib,mrcounts,by.x = "ID", by.y="ID"))
row.names(input_tab) <- input_tab$ID

d <- DGEList(counts=input_tab[c(4:ncol(input_tab))], group=design_tab$condition, genes=input_tab[c(1,2,3)]) # edgeR DGE object
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
logCPMc<-removeBatchEffect(logCPM, design_tab$patient)

plotMDS(logCPM, col=cols, main="MDS Plot: with batch effect")
#legend("topleft", legend=c("Control", "Drug"), col=c(3,4), pch=15)
plotMDS(logCPMc, col=cols, main="MDS Plot: removed batch effect")
#legend("topleft", legend=c("Inf#1", "Inf#2"), col=c(1,2), pch=15)

par(cex.main=0.8)
heatmap.2(cor(logCPM),trace = 'none',density.info = 'none',cexRow = 0.8,cexCol = 0.8,offsetRow = -0.2,offsetCol = -0.2,notecol="black",cellnote = round(cor(logCPM),3),notecex=0.7,
          main = "Correlation of samples with batch effect")

heatmap.2(cor(logCPMc),trace = 'none',density.info = 'none',cexRow = 0.8,cexCol = 0.8,offsetRow = -0.2,offsetCol = -0.2,notecol="black",cellnote = round(cor(logCPMc),3),notecex=0.7,
          main = "Correlation of samples with removed batch effect")
par(cex.main=1.0)

if(paired){
  print("Using paired design in edgeR")
  design = model.matrix(~patient+condition, data=design_tab)
}else{
  print("Using simple design in edgeR")
  design = model.matrix(~condition, data=design_tab)
}

#if we have no replicates, we need to somehow fake dispersions --->>> use EXACT test
#Typical values for the common BCV (square-rootdispersion) for datasets arising from well-controlled experiments are 0.4 for human data,
#0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
if(any(as.data.table(design_tab)[,.N,by=condition]$N == 1)) {
  print("This condition should never be true (only one replicate)")
  bcv <- 0.2
  xglm <- estimateDisp(d, design)
  xglm$common.dispersion <- bcv^2
  xglm$trended.dispersion <- bcv^2
  xglm$tagwise.dispersion <- bcv^2
  
  fit <- glmFit(xglm, design)
  lrt <- glmLRT(fit, coef=ncol(fit$design))
} else{
  xglm <- estimateDisp(d, design)
  # sqrt(xglm$common.disp)
  # plotBCV(xglm, main="BCV Plot")
  
  #use the function glmFit to fit the sgRNA-specific models and glmLRT to do the testing
  #between the drug treated and control samples. The top ranked sgRNAs are listed using the
  #topTags function
  fit <- glmFit(xglm, design)
  lrt <- glmLRT(fit)
}
par(mfrow=c(1,1))
thresh <- 0.05
lfc <- 1
top4 <- topTags(lrt, n=Inf)
top4ids <- rownames(top4$table[abs(top4$table$logFC)>lfc & top4$table$FDR<thresh,])
plotSmear(lrt, de.tags=top4ids, pch=20, cex=0.6,
          main=paste0("logFC per single sgRNAs (lfc > ",lfc,", fdr < ",thresh,")"))
abline(h=c(-1, 0, 1), col=c("dodgerblue","yellow","dodgerblue"), lty=2)

write.table(top4, paste0(output_dir,"sgRNAs_summary.tsv"), sep = "\t",row.names = F, quote = F)

genesymbols <- d$genes[,3]
genesymbollist <- list()
unq <- unique(genesymbols)
unq <- unq[!is.na(unq)]
for(i in unq) {
  sel <- genesymbols==i & !is.na(genesymbols)
  if(sum(sel)>=1)
    genesymbollist[[i]] <- which(sel)
}

camera.res <- camera(xglm, index=genesymbollist, design)

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

colnames(camera.res) <- c("Gene","NsgRNAs","Direction","PValue","FDR","sgRNAs","logFC_per_sgRNA","P_values_per_sgRNA",paste0(colnames(d$counts),"_raw_counts"),paste0(colnames(d$counts),"_logCPM_normalized_counts"))
setorder(camera.res,FDR)
fwrite(camera.res, paste0(output_dir,"gene_summary.tsv"), sep = "\t", row.names = F)

if(camera.res[FDR<=0.05, .N] > 0) {
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
}
dev.off()


