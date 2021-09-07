library(data.table)
library(gplots)

# setwd("/mnt/ssd/ssd_1/snakemake/")
# args <- c("stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/all_samples_report",
#           "/mnt/ssd/ssd_3/references/general/CRISPR_Brunello/CRISPR_Brunello_mod.csv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D22_WT_rep2/D22_WT_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D7_dA5_rep2/D7_dA5_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D7_WT_rep2/D7_WT_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D0_WT_rep1/D0_WT_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D0_dA5_rep1/D0_dA5_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D22_dA5_rep1/D22_dA5_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D0_WT_rep2/D0_WT_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D22_WT_rep1/D22_WT_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D22_dA5_rep2/D22_dA5_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D7_dA5_rep1/D7_dA5_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D7_WT_rep1/D7_WT_rep1_unique_inserts_tab_counts_resolved_duplicates.tsv",
#           "stage372_CRISPR_test/CRISPR_general_analysis/hsblastn_filter/D0_dA5_rep2/D0_dA5_rep2_unique_inserts_tab_counts_resolved_duplicates.tsv")

args = commandArgs(trailingOnly=TRUE)

output_file <-args[1]
unique_inserts <- fread(args[2], header=T)
unique_inserts[,seq := NULL]
sample <- fread(args[3], header=T)
sample[,seq := NULL]
sample_name <- basename(dirname(args[3])) # we use sample name which is the name of the parent directory
created_names <- sample_name
setnames(sample,"N",sample_name)

all_samples <- sample[unique_inserts,on=.(gene_id,UID)]

if(length(args) > 3){ 
  for (i in 4:length(args)){
    sample <- fread(args[i], header=T)
    sample[,seq := NULL]
    sample_name <- basename(dirname(args[i])) # we use sample name which is the name of the parent directory
    created_names <- c(created_names, sample_name)
    setnames(sample,"N",sample_name)
    all_samples <- sample[all_samples,on=.(gene_id,UID)]
  } 
}
all_samples[is.na(all_samples)] <- 0
created_names <- created_names[order(created_names)]
setnames(all_samples, c("gene_id", "UID"), c("Gene", "sgRNA"))
all_samples <- all_samples[order(sgRNA)][,.SD,.SDcols=c("sgRNA","Gene",created_names)]
write.table(all_samples, paste0(output_file,".tsv"), sep = "\t", row.names = F, quote = F, col.names=T)

pdf(paste0(output_file,".pdf"))

a = "PDF report"
b = "This report was done for samples in following table."
plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(0.5,0.7,a, pos=1,cex = 2)
text(0.5,0.3,b, pos=1,cex = 0.9)

textplot(created_names, cex = 0.7)

###################################################
##Normalized read count distribution of all samples

a = "Normalized read count distribution of all samples"
b = "The following figure shows the distribution of median-normalized read counts in all samples."
plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(0.5,0.7,a, pos=1,cex = 1.1)
text(0.5,0.3,b, pos=1,cex = 0.7)


boxplot(log2(all_samples[,.SD,.SDcols=created_names]+1), pch='.', las=2, ylab='log2(read counts)', cex.axis=0.8)


###############################
###DISTRIBUTIOIN OF READ COUNTS
a = "The following figure shows the histogram of median-normalized read counts in all samples."
plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(0.5,0.5,a, pos=1,cex = 0.7)

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F");
tabsmat=as.matrix(log2(all_samples[,.SD,.SDcols=created_names]+1))
colnames(tabsmat) = created_names
samplecol = colors[((1:ncol(tabsmat)) %% length(colors)) ]
if(ncol(tabsmat)>=1){
  histlist=lapply(1:ncol(tabsmat),function(X){ return (hist(tabsmat[,X],plot=F,breaks=40)) })
  xrange=range(unlist(lapply(histlist,function(X){X$mids})))
  yrange=range(unlist(lapply(histlist,function(X){X$counts})))
  hst1=histlist[[1]]
  plot(hst1$mids,hst1$counts,type='b',pch=20,xlim=c(0,xrange[2]*1.2),ylim=c(0,yrange[2]*1.2),xlab='log2(counts)',ylab='Frequency',main='Distribution of read counts',col = samplecol[1] )
}

if(ncol(tabsmat)>=2){ 
  for(i in 2:ncol(tabsmat)){
    hstn=histlist[[i]]
    lines(hstn$mids,hstn$counts,type='b',pch=20,col=samplecol[i])
  }
}
legend('topright',colnames(tabsmat),pch=20,lwd=1,col=samplecol)



#########################
###PCA COMPONENT ANALYSIS
if(length(args) > 3){
  a = "Principle Component Analysis"
  b = "The following figure shows the first 2 principle components (PCs) from the Principle Component Analysis
  (PCA)"
  plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(0.5,0.7,a, pos=1,cex = 1.1)
  text(0.5,0.3,b, pos=1, cex = 0.7)
  
  panel.plot<-function(x,y,textnames=names(x),...){
    par(new=TRUE)
    m<-cbind(x,y)
    plot(m,pch=20,xlim = range(x)*1.1,ylim=range(y)*1.1,...)
    text(x,y,textnames,...)
  }
  
  slmat = as.matrix(log2(all_samples[,.SD,.SDcols=created_names]+1))
  slmat_log=log2(slmat+1)
  ctfit_tx<<-prcomp(t(slmat_log),center=TRUE)
  
  samplecol=colors[((1:ncol(slmat)) %% length(colors)) ]
  
  if(length(samplecol)>2){
    pairs(ctfit_tx$x[,1:3],panel=panel.plot,textnames=rownames(ctfit_tx$x),main='First 3 principle components',col=samplecol)
  }else{
    if(length(samplecol)>1){
      pairs(ctfit_tx$x[,1:2],panel=panel.plot,textnames=rownames(ctfit_tx$x),main='First 2 principle components',col=samplecol)
    }
  }

  #############################
  ##PCA VARIANCE
  a = "The percentage of variances explained by the top PCs."
  plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(0.5,0.5,a, pos=1,cex = 0.7)
  
  
  varpca=ctfit_tx$sdev^2
  varpca=varpca/sum(varpca)*100
  if(length(varpca)>10){
    varpca=varpca[1:10]
  }
  plot(varpca,type='b',lwd=2,pch=20,xlab='PCs',ylab='% Variance explained')

  #######################################
  ######Heatmap analysis
  a = "Sample clustering"
  b = "The following figure shows the sample clustering result."
  plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(0.5,0.7,a, pos=1,cex = 1.1)
  text(0.5,0.3,b, pos=1,cex = 0.7)
  
  slmat=as.matrix(all_samples[,.SD,.SDcols=created_names])
  slmat_log=log2(slmat+1)
  
  result=tryCatch({
    heatmap.2(cor(slmat_log),trace = 'none',density.info = 'none',cexRow = 0.8,cexCol = 0.8,offsetRow = -0.2,offsetCol = -0.2,notecol="black",cellnote = round(cor(slmat_log),3),notecex=0.7)
  }, error=function(e){
    heatmap(cor(slmat_log),scale='none',cexRow = 0.8,cexCol = 0.8,cex.axis=0.8)
  })
}
dev.off()



