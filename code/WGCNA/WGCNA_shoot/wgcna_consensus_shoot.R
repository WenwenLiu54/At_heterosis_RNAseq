#------------------------------------------------------------------------------------

# wgcna_consensus_shoot.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-12-15

#------------------------------------------------------------------------------------

#################
# This script provides the code to generate Consensus CoExpression analysis of shoot RNA-seq

# reference:
# Matt Zinkgraf, US Forest Service and UC Davis Computer Science, new phy paper
#################

setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA2")

library(edgeR)
library(WGCNA)
require(plyr)
library(tidyr)
library(RColorBrewer)
library(preprocessCore)
library(fields)
source("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/references/ConsensusCoExpression-master/R/networkFunctions-extras-05.R")
source("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/references/ConsensusCoExpression-master/R/functions.R")

#################
# set parameters
#################

allowWGCNAThreads(n=63)

options(stringsAsFactors = FALSE)

type = "signed"

corType = "bicor"

corFnc = ifelse(corType=="pearson", "cor", "bicor")

maxPOutliers = ifelse(corType=="pearson",1,0.05)

robustY = ifelse(corType=="pearson",T,F)

# output dir

outdir1 <- "shoot_consensus_result"
if (file.exists(outdir1) == FALSE){
  dir.create(outdir1)  
}

outdir2 <- "shoot_consensus_result/process"
if (file.exists(outdir2) == FALSE){
  dir.create(outdir2)  
}


#################
# load gene expression data
#################

# Col
Data1 <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_C/TPM_log_shoot_C_filter.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE, row.names=1, check.names = FALSE)

# Per
Data2 <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_P/TPM_log_shoot_P_filter.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE, row.names=1, check.names = FALSE)

# F1
Data3 <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_F/TPM_log_shoot_F_filter.txt", sep="\t", stringsAsFactors=FALSE, head=TRUE, row.names=1, check.names = FALSE)


# combine datasets
data <- merge(Data1, Data2, by.x="row.names", by.y="row.names")
data <- merge(data, Data3, by.x="Row.names", by.y="row.names")
row.names(data) <- data[,1]
data <- data[,-1]


#################
# load module information from individual analyses
#################

# Col
C_output <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_C/gene_module_color.txt", sep="\t")
row.names(C_output) <- C_output[,1]
C_mods <- C_output[row.names(data),3]

# Per
P_output <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_P/gene_module_color.txt", sep="\t")
row.names(P_output) <- P_output[,1]
P_mods <- P_output[row.names(data),3]

# F1
F_output <- read.table("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/TPM_log_shoot_F/gene_module_color.txt", sep="\t")
row.names(F_output) <- F_output[,1]
F_mods <- F_output[row.names(data),3]

save(data, C_mods, P_mods, F_mods, file = paste(outdir2, "Grand_analysis_expressions_mods.RData", sep = "/"))


#################
# load expression and modules
#################

load(file = paste(outdir2, "Grand_analysis_expressions_mods.RData", sep = "/"))

lib_names <- names(data)

# build multiset
setLabels = c("Col","Per","F1")

multiExpr = list(Col = list(data=t(data[,1:18])), Per = list(data=t(data[,19:36])), F1 = list(data=t(data[,37:54])))

multiColor = list(Col = C_mods, Per = P_mods, F1 = F_mods)

# Basic sizes of data
size = checkSets(multiExpr);
nSamples = size$nSamples;
nGenes = size$nGenes;
nSets = size$nSets;


#################
# Get soft threshold for individual experiments
#################

powers = c(c(1:10), seq(from = 12, to=30, by=2))
powerTables = vector(mode= "list", length = nSets)
for(set in 1:nSets){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, 
                                                     powerVector = powers, 
                                                     networkType = type,
                                                     corFnc = corFnc,
                                                     removeFirst = FALSE,
                                                     verbose = 2)[[2]])
}

# Save the results
save(powerTables, file = paste(outdir2, "scaleFreeAnalysis-powerTables.RData", sep = "/"))

load(file = paste(outdir2, "scaleFreeAnalysis-powerTables.RData", sep = "/"))
collectGarbage()

# Re-format results for plotting
meanK = modelFit = matrix(0, length(powers), nSets)
for (set in 1:nSets)
{
  modelFit[, set] = -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2];
  meanK[,set] = powerTables[[set]]$data[,5];
}
# Plot scatterplots of topology indices vs. soft-thresholding power
colors = c("black", "red", "blue")


#################
# Plot Figure S1: Soft threshold analysis of individual data sets used in the consensus network
#################

pdf(file = paste(outdir1, "FigS1_MultiExp_scaleFreeTopologyAnalysis.pdf", sep = "/"), wi = 8, h=6)
par(mfrow = c(1,2))

# Plot of Scale independence
plot(powers, modelFit[, 1],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     pch = 21, col = 1, bg = 1,
     main = "Scale independence",
     ylim = range(modelFit))
addGrid()
abline(h=0.85, col="red") 
for(set in 2:nSets){
  points(powers, modelFit[, set], pch = 21, col = colors[set], bg = colors[set]);
  legendClean("bottomright", legend = setLabels, pch = 21, col = colors);
}
  
# Plot of mean connectivity
plot(powers, meanK[, 1],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     main = "Mean connectivity", pch = 21, col = 1, bg = 1)
addGrid()
for(set in 2:nSets){
  points(powers, meanK[, set], pch = 21, col = colors[set], bg = colors[set]);
  legendClean("topright", legend = setLabels, pch = 21, col = colors);
}

dev.off()

#################
# calculate consensus network
#################

STPowers = c(18,18,18) 

TOMinfo_vas = blockwiseIndividualTOMs(multiExpr,
                                      
                                      maxBlockSize = 40000,
                                      
                                      corType = corType,
                                      maxPOutliers = maxPOutliers,
                                      
                                      power = STPowers,
                                      networkType = type,
                                      
                                      TOMType = type,
                  
                  					  saveTOMs = TRUE,
                                      individualTOMFileNames = paste(outdir2, "individualTOM-Set%s-Block%b.RData", sep = "/"),
                                      
                                      nThreads = 63)

print(system.time( {
    mods_vas = blockwiseConsensusModules(
    multiExpr,
    
    maxBlockSize = 40000, 
    
    individualTOMInfo = TOMinfo_vas,
    
    corType = corType,
    maxPOutliers = maxPOutliers,
    
    power = STPowers,
    networkType = type, 
    checkPower = TRUE,
    
    TOMType = type,
    
    consensusQuantile = 0,    
    
    saveConsensusTOMs = FALSE,
    
    deepSplit = 2,  
    detectCutHeight = 0.99, 
    minModuleSize = 50,  
    
    pamRespectsDendro = FALSE,
    
    reassignThresholdPS = 0,
    
    mergeCutHeight = 0.25, 
    
    numericLabels = TRUE,
    
    nThreads = 63,
    verbose = 3, 
    indent = 2,
    
    saveTOMs = FALSE,
    getTOMScalingSamples = TRUE,
)
} ) )

collectGarbage()

# Save the results for future use
save(mods_vas, TOMinfo_vas, file = paste(outdir2, "Grand_mods.RData", sep = "/"))

load(file = paste(outdir2, "Grand_mods.RData", sep = "/"))

# merge similar modules
mergeCut = 0.25 

merge = mergeCloseModules(multiExpr, 
                         mods_vas$unmergedColors, 
                         cutHeight = mergeCut, 
                         consensusQuantile = 0.25, 
                         getNewUnassdME = TRUE, 
                         relabel = TRUE)
labels = merge$colors
moduleColors = labels2colors(labels)

# Save the resulting labels for future use
save(merge, labels, moduleColors, file = paste(outdir2, "merge-labels.RData", sep = "/"))
load(file = paste(outdir2, "merge-labels.RData", sep = "/"))

################
# Plot Figure 2A: Dendrogram and Consensus module identification
################

pdf(file = paste(outdir1, "Fig2A_consensusDendro.pdf", sep = "/"), wi=10, 5)
plotDendroAndColors(mods_vas$dendrograms[[1]], labels2colors(labels[mods_vas$goodGenes]),
                    "Consensus",
                    main ="Consensus modules",
                    dendroLabels = FALSE,
                    addGuide = TRUE, hang = 0.01, 
                    abHeight = 0.99,  # same as above
                    guideHang = 0.03, marAll =c(1,5,1,1))
dev.off()


################
# Create Eigengene network for consensus modules
################

MEs0 = multiSetMEs(multiExpr, universalColors = labels2colors(labels))

MEs <- orderMEs(rbind(MEs0[[1]]$data,MEs0[[2]]$data,MEs0[[3]]$data))
row.names(MEs) <- lib_names

# Output module eigengenes
write.table(MEs, file = paste(outdir1, "CM_Eigengene.txt", sep = "/"), sep="\t", row.names = TRUE, quote = FALSE)

# Cluster module eigengenes
if (corType == "pearson") {
  MEDiss = 1-cor(MEs)
}else{
  MEDiss = 1-bicor(MEs)
}     
METree = hclust(as.dist(MEDiss), method = "average")


###############
# Plot Figure 2B: Consensus and Individual Eigenegene dendrograms
################

pdf(file = paste(outdir1, "Fig2B_Experiments_ME_Dendro.pdf", sep = "/"), w=10, h=2.8)
par(mfrow = c(1, 4), oma = c(1, 1, 1, 1))

plot(as.dendrogram(METree), horiz=FALSE, main = "Consensus", xlab = "", sub = "", cex = 0.6)

dataset_labels <- c("Col-0", "Per-1", "F1") 

for(r in 1:3)
{
  if (corType == "pearson") {
    tmpMEDiss = 1-cor(orderMEs(MEs0[[r]]$data))
  }else{
    tmpMEDiss = 1-bicor(orderMEs(MEs0[[r]]$data))
  }  
  # Cluster module eigengenes
  tmpMETree = hclust(as.dist(tmpMEDiss), method = "average");
  plot(as.dendrogram(tmpMETree), horiz=FALSE, main = dataset_labels[r], xlab = "", sub = "", cex = 0.6)
}
dev.off()


################
# Plot Figure 2C: ME correlations for consensus and individual experiments
################
pdf(file = paste(outdir1, "Fig2C_Experiment_MEcor.pdf", sep = "/"), w=16, h=4.3)
par(mfrow = c(1, 4), oma = c(2, 2, 2, 2), mar = c(3, 3, 3, 3), cex = 0.5)
# order color
c <- unlist(as.dendrogram(METree))
oMEs0 <- orderMEs(MEs0, order=c)
PlotExpPCsCor(oMEs0, dataset_labels, ColorLabels = TRUE, IncludeSign = TRUE, setMargins = FALSE, IncludeGrey = FALSE, plotConsensus = TRUE)
dev.off()


################
# Plot Figure 2D: Identify global and experiment specific co-expression relationships
################
pdf(file=paste(outdir1, "Fig2D_Consensus_module_ExperMod_heatmap.pdf", sep = "/"), w=7, h=2)
par(mfrow = c(1, nSets), oma = c(1, 1, 1, 1), mar = c(1, 1, 1, 1), cex = 0.5)
for(l in 1:nSets)
{
  # Convert the numeric module labels to color labels
  TWColors = multiColor[[l]][mods_vas$goodGenes]
  TWModules <- unique(TWColors)
  consColors = labels2colors(labels[mods_vas$goodGenes])
  consModules <- unique(consColors)
  
  # Numbers of individual and consensus modules
  nTWMods = length(TWModules)
  nConsMods = length(consModules)
  
  # Initialize tables of p-values and of the corresponding counts
  pTable = matrix(0, nrow = nTWMods, ncol = nConsMods);
  CountTbl = matrix(0, nrow = nTWMods, ncol = nConsMods);
  # Execute all pairwaise comparisons
  for (fmod in 1:nTWMods)
  {
    for (cmod in 1:nConsMods)
    {
      TWMembers = (TWColors == TWModules[fmod]);
      consMembers = (consColors == consModules[cmod]);
      pTable[fmod, cmod] = -log10(fisher.test(TWMembers, consMembers, alternative = "greater")$p.value)
      CountTbl[fmod, cmod] = sum(TWColors == TWModules[fmod] & consColors == consModules[cmod])
    }
  }
  
  # Truncate small p values
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
  pTable[pTable>100 ] = 100 
  
  # Marginal counts (really module sizes)
  TWModTotals = apply(CountTbl, 1, sum)
  consModTotals = apply(CountTbl, 2, sum)
  
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  gry <- which(consModules=="grey")
  labeledHeatmap(Matrix = t(pTable[,-gry]),
                 yLabels = paste(" ", consModules[-gry]),
                 xLabels = paste(" ", TWModules),
                 colorLabels = TRUE,
                 ySymbols = consModules[-gry],
                 xSymbols = TWModules,
                 colors = blueWhiteRed(100)[50:100],
                 main = dataset_labels[l],
                 cex.text = 0.8, cex.lab = 0.8, setStdMargins = FALSE)
}
dev.off()


################
# CM eigengene heatmap
################
CM_Eigen <- read.table(file = paste(outdir1, "CM_Eigengene.txt", sep = "/"), sep="\t", row.names = NULL, head = TRUE)

colnames(CM_Eigen) <- gsub("ME", "", colnames(CM_Eigen))
CM_order <- gsub("ME", "", colnames(oMEs0[["Col"]][["data"]]))

CM_Eigen <- CM_Eigen[, c("row.names", CM_order[which(CM_order != "grey")])]

CM_Eigen_1 <- gather(data = CM_Eigen, key = module, value = Eigen_expression, -row.names)

CM_Eigen_1$row.names <- gsub("C", "_C", CM_Eigen_1$row.names)
CM_Eigen_1$row.names <- gsub("P", "_P", CM_Eigen_1$row.names)
CM_Eigen_1$row.names <- gsub("F", "_F", CM_Eigen_1$row.names)

CM_Eigen_2 <- separate(CM_Eigen_1, col = row.names, into = c("timepoint", "sample"))

CM_Eigen_2$module_sample <- mapply(function(x, y){return(paste(x, y, sep = "_"))}, CM_Eigen_2$module, CM_Eigen_2$sample)

CM_Eigen_3 <- CM_Eigen_2[, c(1,4,5)]

CM_Eigen_4 <- spread(data = CM_Eigen_3, key = module_sample, value = Eigen_expression)

m <- matrix(CM_Eigen_3$module_sample, ncol = nSamples[1], byrow = TRUE)
m1 <- m[,c(1:3)]

heat_order <- as.vector(t(m1))

CM_Eigen_5 <- CM_Eigen_4[,c("timepoint", heat_order)]

# plot CM Eigengene heatmap
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(require(circlize)) 
suppressPackageStartupMessages(library(dendextend))
library(tidyverse)

mat <- as.matrix(CM_Eigen_5[, -1])
rownames(mat) <-  CM_Eigen_5[,"timepoint"]

sample_info <- data.frame("heat_order" = heat_order, module = str_sub(heat_order, 1, -4), genotype = str_sub(heat_order, -2, -2))
table(sample_info$module)
table(sample_info$genotype)

top_an.1 <- HeatmapAnnotation(df = sample_info[, c(2,3)], 
                              show_annotation_name = TRUE, 
                              col = list("module" = c("brown" = "brown", "greenyellow" = "greenyellow", "blue" = "blue", "green" = "green", "yellow" = "yellow", "black" = "black", "magenta" = "magenta", "turquoise" = "turquoise", "red" = "red", "pink" = "pink", "purple" = "purple"), 
                                        "genotype"=c("C"="gold","P"="lightblue","F"="limegreen")
                                        ), 
                            show_legend = TRUE,  
                            which = "column", 
                            annotation_height=unit(c(1, 4), "cm"), 
                            gap=unit(1, "mm"))

pdf(file = paste(outdir1, "Fig2E_CM_Eigengene_heatmap.pdf", sep = "/"), width = 20, height = 12)
h <- Heatmap(mat, col = colorRamp2(c(min(mat), (min(mat)+max(mat))/2, max(mat)), c("blue", "white", "red")), name = "CM_heatmap", cluster_rows = FALSE, cluster_columns = FALSE, row_title = "timepoint",  column_names_side = "top",row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize =12), column_names_gp = gpar(fontsize = 12), column_title_side= "bottom", heatmap_legend_param = list(title= "eigengene expression", title_position ="leftcenter-rot", legend_height=unit(2,"cm"), legend_direction="vertical"), row_names_max_width = unit(100, "mm"), heatmap_width = unit(40, "cm"), heatmap_height = unit(16, "cm"), column_split = factor(rep(CM_order[which(CM_order != "grey")], each=9), levels = CM_order), column_gap = unit(.5, "mm"), top_annotation = top_an.1)
draw(h, heatmap_legend_side = "left")
dev.off()



#################
# kME
# reference: 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-RelateModsToTraits.R
#################

probes = colnames(multiExpr[[1]]$data)

moduleColors = labels2colors(labels)

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleColors, excludeGrey = TRUE)

kME = list()

for (set in 1:nSets) {
  if (corType == "pearson") {
    kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
  }else{
    kME[[set]] = bicorAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
  }
}

# Code for kME
if (corType == "pearson") {
  kMEmat = rbind(kME[[1]]$cor, kME[[1]]$p, kME[[2]]$cor, kME[[2]]$p, kME[[3]]$cor, kME[[3]]$p)
}else{
  kMEmat = rbind(kME[[1]]$bicor, kME[[1]]$p, kME[[2]]$bicor, kME[[2]]$p, kME[[3]]$bicor, kME[[3]]$p)
}    

MEnames = colnames(consMEs.unord[[1]]$data)
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes
colnames(kMEmat) = spaste(c("kME.Col.", "p.kME.Col.", "kME.Per.", "p.kME.Per.", "kME.F1.", "p.kME.F1."), rep(gsub("ME", "", MEnames), rep(6, nMEs)))

info = data.frame(Probe = probes, 
                  ModuleLabel = labels,
                  ModuleColor = labels2colors(labels),
                  kMEmat)

write.csv(info, file = paste(outdir1, "consensusAnalysis-CombinedNetworkResults_kME.csv", sep = "/"), row.names = FALSE, quote = FALSE)


################
# CM expression heatmap
################

pdf(file = paste(outdir1, "Fig2F_CM_expression_heatmap.pdf", sep = "/"), height=20, width=120)
par(mfrow=c(3, length(CM_order[which(CM_order != "grey")])), mar=c(1, 3, 8, 1), cex=2)
for (set in c(1:nSets)) {
  for (which.module in CM_order[which(CM_order != "grey")])
  {
    which.module
    plotMat(t(scale(multiExpr[[set]]$data[, moduleColors==which.module])), nrgcols=1000, clabels=rownames(multiExpr[[set]]$data), rcols = which.module, main= paste (paste("CM", which.module, sep = "-"), setLabels[set], sep = " of "), cex.main=1.5, cex.lab=2, cex.axis=2)
  }
}
dev.off()

#################
# CM GO
#################
# get GeneID list
module_df <- info[,c(1:3)]
colnames(module_df) <- c("GeneID", "module_ID", "module_color")
rownames(module_df) <- NULL
ID_list <- list()
module_df$GeneID <- as.character(module_df$GeneID)
for (which.module in c(CM_order[which(CM_order != "grey")],"grey")){
  ID_list[[which.module]] <- module_df$GeneID[which(module_df$module_color == which.module)]
}

# use clusterProfiler for GO analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))

# enricher
GO_list <- list()
for (which.module in c(CM_order[which(CM_order != "grey")],"grey")){
  GENE_ID <- ID_list[[which.module]]
  ego <- enrichGO(gene = GENE_ID, 
                  OrgDb = "org.At.tair.db",
                  keyType="TAIR",
                  ont="BP",
                  pAdjustMethod = "fdr" 
  )
  ego <- simplify(ego)
  GO_list[[which.module]] <- ego@result
}

# Save the result for future use
save(ID_list, GO_list, file = paste(outdir2, "ID_GO_list.RData", sep = "/"))
load(file = paste(outdir2, "ID_GO_list.RData", sep = "/"))

# filter by p.adjust
for (which.module in c(CM_order[which(CM_order != "grey")],"grey")) {
  GO_list[[which.module]] <- GO_list[[which.module]][which(GO_list[[which.module]][,6] < 0.05),]
}
# extract Term name and p.adjust
for (which.module in c(CM_order[which(CM_order != "grey")],"grey")) {
  GO_list[[which.module]] <- GO_list[[which.module]][,c(2,6)]
  colnames(GO_list[[which.module]]) <- c("Term", which.module)
}
# merge
temp <- GO_list[[1]]
for (i in names(GO_list)[-1]) {
  df2 <- GO_list[[i]]
  temp <- merge(x = temp, y = df2, by = "Term", all.x = TRUE, all.y = TRUE)  
}

save(temp, file = paste(outdir2, "GO_merge.RData", sep = "/"))
write.table(temp, file = paste(outdir1, "shoot_CM_GO_merge.txt", sep = "/"),
            sep = "\t", row.names = FALSE, col.names = TRUE)
# plot heat map
temp[,-1][!is.na(temp[,-1])] <- -log10(temp[,-1][!is.na(temp[,-1])])
##sort
library(tidyverse)
temp <- dplyr::arrange(temp, desc(brown), desc(greenyellow), desc(blue), desc(green), desc(yellow), desc(black), desc(magenta), desc(turquoise), desc(red), desc(pink), desc(purple), desc(grey))
##plot
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(require(circlize)) 
suppressPackageStartupMessages(library(dendextend))
mat <- as.matrix(temp[, -1])
rownames(mat) <-  temp[,"Term"]

##save plot
top_an.2 <- HeatmapAnnotation(module = c(CM_order[which(CM_order != "grey")],"grey"), 
                              show_annotation_name = TRUE, 
                              col = list(module = c("brown" = "brown", "greenyellow" = "greenyellow", "blue" = "blue", "green" = "green", "yellow" = "yellow", "black" = "black", "magenta" = "magenta", "turquoise" = "turquoise", "red" = "red", "pink" = "pink", "purple" = "purple", "grey" = "grey")), 
                              show_legend = TRUE,  
                              which = "column" 
)

pdf(file = paste(outdir1, "Fig2G_shoot_CM_GO.pdf", sep = "/"), width = 6, height = 16)
h <- Heatmap(mat, col = colorRamp2(c(1.3, 5), c("yellow", "red")), name = "-log10(p.adjust)", na_col = "grey", cluster_rows = FALSE, cluster_columns = FALSE, row_title = "GO Terms",  column_names_side = "top",row_title_gp = gpar(fontsize = 10), column_title_gp = gpar(fontsize = 3), row_names_gp = gpar(fontsize =2), column_names_gp = gpar(fontsize = 3), column_title_side= "bottom", heatmap_legend_param = list(title= "-log10(p.adjust)", title_position ="leftcenter-rot", legend_height=unit(2,"cm"), legend_direction="vertical"), row_names_max_width = unit(100, "mm"), heatmap_width = unit(6, "cm"), heatmap_height = unit(40, "cm"), top_annotation = top_an.2)
draw(h, heatmap_legend_side = "left")
dev.off()

#################
# hub gene GO
#################
# hub gene ID
hub_ID_list <- list()
for (set in c("Col", "Per", "F1")) {
  hub_ID_list[[set]] <- list()
  for (which.module in CM_order[which(CM_order != "grey")]) {
    kME <- paste(paste0("kME.", set), which.module, sep = ".")
    p.kME <- paste(paste0("p.kME.", set), which.module, sep = ".")
    
    a <- which(info$ModuleColor == which.module)
    b <- which(info[, print(kME)] > 0.9)
    c <- which(info[, print(p.kME)] < 1e-06)
    
    id <- intersect(a, b)
    id <- intersect(id, c)
    
    hub_ID_list[[set]][[which.module]] <- info$Probe[id]
  }
}

# union hub gene ID of 3 genotypes for each CM
hub_ID_list_union <- list()
for (which.module in CM_order[which(CM_order != "grey")]) {
  hub_ID_list_union[[which.module]] <- union(union(hub_ID_list[["Col"]][[which.module]], hub_ID_list[["Per"]][[which.module]]), hub_ID_list[["F1"]][[which.module]])  
}

# union hub gene ID of 3 genotypes (all CMs)
GENE_ID <- c()
for (which.module in CM_order[which(CM_order != "grey")]){
  GENE_ID <- c(GENE_ID, hub_ID_list_union[[which.module]])
}
ego <- enrichGO(gene = GENE_ID,
                OrgDb = "org.At.tair.db",
                keyType="TAIR",
                ont="BP",
                pAdjustMethod = "fdr",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)
ego <- simplify(ego)
ego_result <- ego@result
write.table(ego_result, file = paste(outdir1, "shoot_CM_hub_union_GO.txt", sep = "/"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)
emapplot(ego, 
         showCategory = 164,
         color = "p.adjust",
         layout = "nicely",
         by = count,
         includeAll = TRUE,
         line_scale = 0.8,
         pie_scale = 0.5)

#################
# TOMplot
#################
setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA2/shoot_consensus_result")
for (set in 1 : nSets) {
  TOM = TOMsimilarityFromExpr(multiExpr[[set]]$data, power = STPowers[set], networkType = type, TOMType = type, corType = corType)
  dissTOM = 1-TOM 
  plotTOM = dissTOM^7; 
  diag(plotTOM) = NA; 
  
  geneTree = mods_vas$dendrograms[[1]] 
  moduleColors = labels2colors(labels[mods_vas$goodGenes])
  
  library(gplots)
  png(file = paste0("TOM_heatmap_", set, ".png"),width = 800,height = 800)
  TOMplot(plotTOM, geneTree, moduleColors, col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
  dev.off()
}