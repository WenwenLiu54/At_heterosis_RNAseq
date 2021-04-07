#------------------------------------------------------------------------------------

# RNAseq_part3_tximport.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-11-17

#------------------------------------------------------------------------------------

setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r")

#------------------------------------------------------------------
# library packages
library(edgeR)
library(limma)
library(tximport)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))

#------------------------------------------------------------------
# Defined functions
### filter cpm
edgeR_filter<-function(mat,cutoff=40,...){
  rowSums(edgeR::cpm(mat,...)>1)>=cutoff
}

#------------------------------------------------------------------
#load data file names
samples<-list.files('/data2/usr/LiuWW/project_1/try/RNAseq_part2_salmon_output_2/')
samples
files<-paste0('/data2/usr/LiuWW/project_1/try/RNAseq_part2_salmon_output_2/', samples, '/', samples, '_transcripts_quant/quant.sf')
quant <- read.delim("/data2/usr/LiuWW/project_1/try/RNAseq_part2_salmon_output_2/10C1/10C1_transcripts_quant/quant.sf")

##trans2genes
trans<-as.vector(quant$Name)
genes<-substr(trans,1,9)
trans2genes<-data.frame(TXNAME=trans,GENEID=genes)

###gene level
txi_genes_lengthScaledTPM <- tximport(files, type = "salmon", tx2gene = trans2genes, countsFromAbundance = "lengthScaledTPM")
load("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/samples1.RData")
colnames(txi_genes_lengthScaledTPM$counts)<-samples1
colnames(txi_genes_lengthScaledTPM$abundance)<-samples1
save(txi_genes_lengthScaledTPM, file='./txi_genes_lengthScaledTPM.RData')

###transcript level
txi_trans <- tximport(files, type = "salmon", tx2gene = NULL, countsFromAbundance = "lengthScaledTPM",txOut = TRUE)
colnames(txi_trans$counts)<-samples1
colnames(txi_trans$abundance)<-samples1
save(txi_trans,file='./txi_trans.RData')

#------------------------------------------------------------------
#get counts and TPM
est_genes_counts_lengthScaledTPM <- txi_genes_lengthScaledTPM$counts
est_trans_counts <- txi_trans$counts
est_genes_TPM_lengthScaledTPM <- txi_genes_lengthScaledTPM$abundance
est_trans_TPM <- txi_trans$abundance

save(est_genes_counts_lengthScaledTPM,file='./est_genes_counts_lengthScaledTPM.RData')
save(est_trans_counts,file='./est_trans_counts.RData')
save(est_genes_TPM_lengthScaledTPM,file='./est_genes_TPM_lengthScaledTPM.RData')
save(est_trans_TPM,file='./est_trans_TPM.RData')
write.table(est_genes_counts_lengthScaledTPM, file="est_genes_counts_lengthScaledTPM.txt", sep="\t")
write.table(est_genes_TPM_lengthScaledTPM, file="est_genes_TPM_lengthScaledTPM.txt", sep="\t")
write.table(est_trans_counts, file="est_trans_counts.txt", sep="\t")
write.table(est_trans_TPM, file="est_trans_TPM.txt", sep="\t")

#------------------------------------------------------------------
# Filter low expressed transcripts
# Low expressed transcripts were filtered based on counts per million (CPM). An expressed transcript must have 3 or more samples out of 189 ≥ 1 CPM. Low expressed genes were the genes with only low expressed transcripts.
trans<-rownames(est_trans_counts)
genes<-unique(substr(trans,1,9))
filter.idx<-edgeR_filter(est_trans_counts,cutoff = 3)
table(filter.idx)
trans_high<-names(filter.idx)[filter.idx==TRUE]
genes_high<-unique(substr(trans_high,1,9))
save(trans_high,file='./trans_high.RData')
save(genes_high,file='./genes_high.RData')
