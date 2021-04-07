#------------------------------------------------------------------------------------

# limma_voom_PvsC.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-11-27

#------------------------------------------------------------------------------------

setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r")

library(limma)
library(edgeR)
library(stringr)

load("./genes_high.RData")
load("./samples2.RData")

#expression matrix
count <- read.table(file="est_genes_counts_lengthScaledTPM.txt",header=T,sep="\t", row.names = 1, check.names = FALSE)
count <- count[,samples2]
count_filter <- count[genes_high,]
write.table(count_filter, file="est_genes_counts_lengthScaledTPM.genes_high.txt", sep="\t")

#make expriment design matrix
data <- read.table("samp.inf.txt",head=TRUE,check.names = FALSE)

x1 <- data$Sample[1:81]
x1_new <- paste (0, x1, sep = "")
x2 <- data$tissue[1:81]
x3 <- data$group[1:81]
x3_new <- paste (0, x3, sep = "")
x4 <- data$strain[1:81]
x5 <- data[82:189,]

samp_info = data.frame(Sample = x1_new, tissue = x2, group = x3_new, strain = x4)
samp_info <- rbind(samp_info, x5)

#make the orders of sample name and group name consistent
temp <- samp_info[c(1,2,3,7,8,9,4,5,6),]
for (i in c(2:21)){
  part1 <- c((9*i-8):(9*i-6))
  part2 <- c((9*i-5):(9*i-3))
  part3 <- c((9*i-2):(9*i))
  new <- samp_info[c(part1, part3, part2),]
  temp <- rbind(temp, new)
}

write.table(temp, file="samp.inf.1.txt", row.names = FALSE, sep="\t")

data1 <- read.table("samp.inf.1.txt", row.names="Sample", head=TRUE, check.names = FALSE)
group <- factor(data1$group,levels=levels(data1$group))
design <- model.matrix(~0+group)
rownames(design) <- rownames(data1)


#make contrast matrix
contr.matrix <- makeContrasts(
    
  contrast_03PvsC=group03P-group03C,

  contrast_04PvsC=group04P-group04C,
    
  contrast_05PvsC=group05P-group05C,
    
  contrast_06PvsC=group06P-group06C,
     
  contrast_07SPvsC=group07SP-group07SC,
     
  contrast_08SPvsC=group08SP-group08SC,
     
  contrast_07LPvsC=group07LP-group07LC,
     
  contrast_08LPvsC=group08LP-group08LC,
     
  contrast_09PvsC=group09P-group09C,
     
  contrast_10PvsC=group10P-group10C,
    
  contrast_11PvsC=group11P-group11C,
    
  contrast_12PvsC=group12P-group12C,
     
  contrast_13PvsC=group13P-group13C,
     
  contrast_14PvsC=group14P-group14C,
    
  contrast_15PvsC=group15P-group15C,
    
  contrast_16PvsC=group16P-group16C,
    
  contrast_17PvsC=group17P-group17C,
    
  contrast_18PvsC=group18P-group18C,
    
  contrast_19PvsC=group19P-group19C,
     
  contrast_20PvsC=group20P-group20C,
    
  contrast_21PvsC=group21P-group21C,
  
  levels = design)
contr.matrix

v <- voom(count_filter, design, plot=TRUE)
# step 1
fit <- lmFit(v, design)
# step 2
fit1 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit1)
# step 3
results <- decideTests(fit2)
colnames(results)


DE<-topTable(fit2,coef=1,adjust='BH',number = 'all')
write.table(DE, "03PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=2,adjust='BH',number = 'all')
write.table(DE, "04PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=3,adjust='BH',number = 'all')
write.table(DE, "05PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=4,adjust='BH',number = 'all')
write.table(DE, "06PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=5,adjust='BH',number = 'all')
write.table(DE, "07SPvsC.txt",sep="\t")

DE<-topTable(fit2,coef=6,adjust='BH',number = 'all')
write.table(DE, "08SPvsC.txt",sep="\t")

DE<-topTable(fit2,coef=7,adjust='BH',number = 'all')
write.table(DE, "07LPvsC.txt",sep="\t")

DE<-topTable(fit2,coef=8,adjust='BH',number = 'all')
write.table(DE, "08LPvsC.txt",sep="\t")

DE<-topTable(fit2,coef=9,adjust='BH',number = 'all')
write.table(DE, "09PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=10,adjust='BH',number = 'all')
write.table(DE, "10PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=11,adjust='BH',number = 'all')
write.table(DE, "11PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=12,adjust='BH',number = 'all')
write.table(DE, "12PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=13,adjust='BH',number = 'all')
write.table(DE, "13PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=14,adjust='BH',number = 'all')
write.table(DE, "14PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=15,adjust='BH',number = 'all')
write.table(DE, "15PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=16,adjust='BH',number = 'all')
write.table(DE, "16PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=17,adjust='BH',number = 'all')
write.table(DE, "17PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=18,adjust='BH',number = 'all')
write.table(DE, "18PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=19,adjust='BH',number = 'all')
write.table(DE, "19PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=20,adjust='BH',number = 'all')
write.table(DE, "20PvsC.txt",sep="\t")

DE<-topTable(fit2,coef=21,adjust='BH',number = 'all')
write.table(DE, "21PvsC.txt",sep="\t")

dev.off()