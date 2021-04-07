#------------------------------------------------------------------------------------

# wgcna_individual_shoot_C.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-12-14

#------------------------------------------------------------------------------------

setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA2")

library(WGCNA)
library(reshape2)
library(stringr)

#############set parameters############# 
options(stringsAsFactors = FALSE)

enableWGCNAThreads()

type = "signed"

corType = "bicor"

corFnc = ifelse(corType=="pearson", "cor", "bicor")

maxPOutliers = ifelse(corType=="pearson",1,0.05)

robustY = ifelse(corType=="pearson",T,F)

exprMat <- "expr_data/TPM_log_shoot_C.txt"

#output dir
if (file.exists(str_sub(basename(paste0(exprMat)), 1, -5)) == FALSE){
  dir.create(str_sub(basename(paste0(exprMat)), 1, -5))  
}
outdir <- str_sub(basename(paste0(exprMat)), 1, -5)


#############input data###############
dataExpr <- read.table(exprMat, sep='\t', row.names = 1, header=TRUE, quote="", comment="", check.names=FALSE)

dim(dataExpr)

head(dataExpr)[,1:8]


#############filter genes#############
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

write.table(dataExprVar, file = paste(outdir, "TPM_log_shoot_C_filter.txt", sep = "/"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep = "\t")
#############transposition#############
dataExpr <- as.data.frame(t(dataExprVar))
dim(dataExpr)

#############detect missing value#############
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

head(dataExpr)[,1:8]


#############detect whether outlier samples existed#############
pdf(file = paste(outdir, "sampleTree.pdf", sep = "/"), height=6, width=32)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

#############pickSoftThreshold#############
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, 
                        powerVector=powers, 
                        networkType=type,
                        corFnc = corFnc,
                        removeFirst = FALSE,
                        verbose=5)

#plot
pdf(file = paste(outdir, "pickSoftThreshold.pdf", sep = "/"), height=6, width=8)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")

abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off() 

power = sft$powerEstimate
power

if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
power


#############network construction#############
# One-step network construction and module detection
net = blockwiseModules(dataExpr,  
                       
                       maxBlockSize = nGenes, 
                       
                       reassignThreshold = 0, 
                       
                       corType = corType, 
                       maxPOutliers = maxPOutliers,
                       
                       power = power, 
                       networkType = type,
                       
                       TOMType = type,
                       
                       saveTOMs = TRUE,  
                       loadTOM = FALSE,
                       saveTOMFileBase = paste0(outdir, ".tom"),
                       
                       deepSplit = 2,
                       detectCutHeight = 0.995, 
                       minModuleSize = 100,  
                       
                       pamStage = TRUE, 
                       pamRespectsDendro = TRUE,   
                       
                       mergeCutHeight = 0.25, 
                      
                       numericLabels = TRUE, 
                       
                       verbose = 3)

table(net$colors)

#############gene number in each module#############
module_size <- as.data.frame(table(net$colors))
colnames(module_size) <- c("module_ID", "gene_number")
write.table(module_size, file = paste(outdir, "module_size.txt", sep = "/"), row.names=FALSE, col.names = TRUE, quote = FALSE, sep="\t")

#############Plot the dendrogram and the module colors#############
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
pdf(file = paste(outdir, "geneTree.pdf", sep = "/"), height=6, width=8) 
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main =paste("Cluster Dendrogram", "of", outdir, sep = " "))
dev.off()

#############get module-gene list################
write.table(paste(colnames(dataExpr), moduleLabels, moduleColors, sep = "\t"), file=paste(outdir, "gene_module_color.txt", sep = "/"), row.names=FALSE, col.names=FALSE, quote=FALSE)

#############TOMplot#######
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = 18, networkType = type, TOMType = type, corType = corType); 
# Taking the dissimilarity to a power
plotTOM = dissTOM^7; 
# setting the diagonal to NA also improves the clarity of the plot
diag(plotTOM) = NA; 
library(gplots)
png(file = paste(outdir, "TOM_heatmap.png", sep = "/"),width = 800,height = 800)
TOMplot(plotTOM, geneTree, moduleColors, col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()