#------------------------------------------------------------------------------------

# get_expr_data_for_wgcna.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-12-4

#------------------------------------------------------------------------------------

setwd("G:/project1/my_result/RNAseq_part3_r/final/WGCNA/expr_data")

#load data
load("G:/project1/my_result/RNAseq_part3_r/final/RData/sample.63.RData")
load("G:/project1/my_result/RNAseq_part3_r/final/RData/sample.189.RData")
load("G:/project1/my_result/RNAseq_part3_r/final/RData/tpm_average.RData")
TPM <- read.table(file="G:/project1/my_result/RNAseq_part3_r/final/est_genes_TPM_lengthScaledTPM.genes_high.txt", header = TRUE, sep = "\t", row.names = NULL)

#TPM of each organ and genotype at each stage with 3 biological replicates
colnames(TPM) <- c("GeneID", sample.189)
TPM <- TPM[order(TPM$GeneID),]

mat_189 <- matrix(sample.189, ncol = 9, byrow = TRUE)

sample.18.C <- as.vector(t(mat_189[c(1:6),c(1:3)]))
sample.18.F <- as.vector(t(mat_189[c(1:6),c(4:6)]))
sample.18.P <- as.vector(t(mat_189[c(1:6),c(7:9)]))

sample.45.C <- as.vector(t(mat_189[c(7:21),c(1:3)]))
sample.45.F <- as.vector(t(mat_189[c(7:21),c(4:6)]))
sample.45.P <- as.vector(t(mat_189[c(7:21),c(7:9)]))

TPM_shoot_C <- TPM[, c("GeneID", sample.18.C)]
TPM_shoot_F <- TPM[, c("GeneID", sample.18.F)]
TPM_shoot_P <- TPM[, c("GeneID", sample.18.P)]

TPM_leaf_C <- TPM[, c("GeneID", sample.45.C)]
TPM_leaf_F <- TPM[, c("GeneID", sample.45.F)]
TPM_leaf_P <- TPM[, c("GeneID", sample.45.P)]

#get the average TPM of biological replicates
TPMave_63 <- tpm_average[, c("GeneID", sample.63)]
TPMave_63 <- TPMave_63[order(TPMave_63$GeneID),]
mat_63 <- matrix(sample.63, ncol = 3, byrow = TRUE)

sample.6.C <- as.vector(t(mat_63[c(1:6),1]))
sample.6.F <- as.vector(t(mat_63[c(1:6),2]))
sample.6.P <- as.vector(t(mat_63[c(1:6),3]))

sample.15.C <- as.vector(t(mat_63[c(7:21),1]))
sample.15.F <- as.vector(t(mat_63[c(7:21),2]))
sample.15.P <- as.vector(t(mat_63[c(7:21),3]))

TPMave_shoot_C <- TPMave_63[, c("GeneID", sample.6.C)]
TPMave_shoot_F <- TPMave_63[, c("GeneID", sample.6.F)]
TPMave_shoot_P <- TPMave_63[, c("GeneID", sample.6.P)]

TPMave_leaf_C <- TPMave_63[, c("GeneID", sample.15.C)]
TPMave_leaf_F <- TPMave_63[, c("GeneID", sample.15.F)]
TPMave_leaf_P <- TPMave_63[, c("GeneID", sample.15.P)]

#save data
save.image("raw_expr.RData")

#############filter genes############
#define function
f<-function(x) sum(x >= 1)
#For shoot：among the average TPM values at 6 stages,at least one value >= 1
#For leaf：among the average TPM values at 15 stages,at least one value >= 1
filein1 <- list(TPMave_shoot_C, TPMave_shoot_F, TPMave_shoot_P, TPMave_leaf_C, TPMave_leaf_F, TPMave_leaf_P)
filein2 <- list(TPM_shoot_C, TPM_shoot_F, TPM_shoot_P, TPM_leaf_C, TPM_leaf_F, TPM_leaf_P)
fileout1 <- c("TPM_shoot_C.txt", "TPM_shoot_F.txt", "TPM_shoot_P.txt", "TPM_leaf_C.txt", "TPM_leaf_F.txt", "TPM_leaf_P.txt")
fileout2 <- c("TPM_log_shoot_C.txt", "TPM_log_shoot_F.txt", "TPM_log_shoot_P.txt", "TPM_log_leaf_C.txt", "TPM_log_leaf_F.txt", "TPM_log_leaf_P.txt")

for (i in c(1:6)){
  num <- apply(filein1[[i]][,-1], 1, f)
  id <- which(num >= 1)
  ############get expression matrix (1) TPM value#############
  filter <- filein2[[i]][id,]
  write.table(filter, file = fileout1[i], row.names = FALSE, sep = "\t")
  ############get expression matrix (2) log2(TPM+1) value#############
  filter <- filein2[[i]][id,]
  filter_log <- cbind(filter[,1],log2(filter[,-1] + 1))
  write.table(filter_log, file = fileout2[i], row.names = FALSE, sep = "\t")
}