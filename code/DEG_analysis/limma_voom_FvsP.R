#------------------------------------------------------------------------------------

# limma_voom_FvsP.R
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
    
  contrast_03FvsP=group03F-group03P,

  contrast_04FvsP=group04F-group04P,
    
  contrast_05FvsP=group05F-group05P,
    
  contrast_06FvsP=group06F-group06P,
     
  contrast_07SFvsP=group07SF-group07SP,
     
  contrast_08SFvsP=group08SF-group08SP,
     
  contrast_07LFvsP=group07LF-group07LP,
     
  contrast_08LFvsP=group08LF-group08LP,
     
  contrast_09FvsP=group09F-group09P,
     
  contrast_10FvsP=group10F-group10P,
    
  contrast_11FvsP=group11F-group11P,
    
  contrast_12FvsP=group12F-group12P,
     
  contrast_13FvsP=group13F-group13P,
     
  contrast_14FvsP=group14F-group14P,
    
  contrast_15FvsP=group15F-group15P,
    
  contrast_16FvsP=group16F-group16P,
    
  contrast_17FvsP=group17F-group17P,
    
  contrast_18FvsP=group18F-group18P,
    
  contrast_19FvsP=group19F-group19P,
     
  contrast_20FvsP=group20F-group20P,
    
  contrast_21FvsP=group21F-group21P,
  
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
write.table(DE, "03FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=2,adjust='BH',number = 'all')
write.table(DE, "04FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=3,adjust='BH',number = 'all')
write.table(DE, "05FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=4,adjust='BH',number = 'all')
write.table(DE, "06FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=5,adjust='BH',number = 'all')
write.table(DE, "07SFvsP.txt",sep="\t")

DE<-topTable(fit2,coef=6,adjust='BH',number = 'all')
write.table(DE, "08SFvsP.txt",sep="\t")

DE<-topTable(fit2,coef=7,adjust='BH',number = 'all')
write.table(DE, "07LFvsP.txt",sep="\t")

DE<-topTable(fit2,coef=8,adjust='BH',number = 'all')
write.table(DE, "08LFvsP.txt",sep="\t")

DE<-topTable(fit2,coef=9,adjust='BH',number = 'all')
write.table(DE, "09FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=10,adjust='BH',number = 'all')
write.table(DE, "10FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=11,adjust='BH',number = 'all')
write.table(DE, "11FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=12,adjust='BH',number = 'all')
write.table(DE, "12FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=13,adjust='BH',number = 'all')
write.table(DE, "13FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=14,adjust='BH',number = 'all')
write.table(DE, "14FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=15,adjust='BH',number = 'all')
write.table(DE, "15FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=16,adjust='BH',number = 'all')
write.table(DE, "16FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=17,adjust='BH',number = 'all')
write.table(DE, "17FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=18,adjust='BH',number = 'all')
write.table(DE, "18FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=19,adjust='BH',number = 'all')
write.table(DE, "19FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=20,adjust='BH',number = 'all')
write.table(DE, "20FvsP.txt",sep="\t")

DE<-topTable(fit2,coef=21,adjust='BH',number = 'all')
write.table(DE, "21FvsP.txt",sep="\t")

dev.off()