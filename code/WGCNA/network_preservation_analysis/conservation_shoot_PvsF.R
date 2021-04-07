#------------------------------------------------------------------------------------

# conservation_shoot_PvsF.R
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2019-12-29

#------------------------------------------------------------------------------------

setwd("/data2/usr/LiuWW/project_1/try/RNAseq_part3_r/WGCNA/conservation_z_score")

library(WGCNA) 
library(flashClust) 


###########prepare data###########

datExprC = read.table("../TPM_log_shoot_C/TPM_log_shoot_C_filter.txt", sep="\t", head=TRUE, row.names=1, check.names=FALSE)
datExprP = read.table("../TPM_log_shoot_P/TPM_log_shoot_P_filter.txt", sep="\t", head=TRUE, row.names=1, check.names=FALSE)
datExprF = read.table("../TPM_log_shoot_F/TPM_log_shoot_F_filter.txt", sep="\t", head=TRUE, row.names=1, check.names=FALSE)

datColor = read.table("../TPM_log_shoot_F/gene_module_color.txt", sep="\t", check.names=FALSE)

modulecolor_ref = datColor[,3]

###########PvsF conservation analysis###########

multiExprA = list(A1 = list(data = t(datExprF)), A2 = list(data = t(datExprP)))

multiColorA = list(A1 = modulecolor_ref)

mp = modulePreservation(multiExprA, multiColorA, referenceNetworks=1, verbose=3, networkType="signed",
                      nPermutations=50, maxGoldModuleSize=100, maxModuleSize=10000)

ref = 1
test = 2
Obs.PreservationStats = mp$preservation$observed[[ref]][[test]]
Z.PreservationStats = mp$preservation$Z[[ref]][[test]]

###########output###########

write.table(Obs.PreservationStats, file = "shoot_PvsF_Obs.PreservationStats.txt", sep="\t")
write.table(Z.PreservationStats, file = "shoot_PvsF_Z.PreservationStats.txt", sep="\t")


pdf("shoot_PvsF_conservation_z_score.pdf", height = 6, width = 10)

moduleSize = Obs.PreservationStats$moduleSize

modColors = rownames(Obs.PreservationStats)

selectModules = !(modColors %in% c("grey", "gold"))


point_label = modColors[selectModules]


cor.kIM=Obs.PreservationStats$cor.kIM  
Zsummary=Z.PreservationStats$Zsummary.pres   

### plot
par(mfrow=c(1,2), mar = c(4.5,4.5,2.5,1))
#### plot cor.kIM vs module size
plot(moduleSize[selectModules], cor.kIM[selectModules], col = 1,
     bg = modColors[selectModules], pch = 21, main = "cor.kIM Preservation",
     cex = 2, ylab = "cor.kIM", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules], cor.kIM[selectModules], point_label, cex=1, offs=0.03)

#### plot Zsummary vs module size
plot(moduleSize[selectModules], Zsummary[selectModules], col = 1,
     bg = modColors[selectModules], pch = 21, main = "Zsummary Preservation",
     cex = 2, ylab = "Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules], Zsummary[selectModules], point_label, cex=1, offs=0.03)
# Add threshold lines for Zsummary
abline(h = 0)
abline(h = 2, col = "blue", lty = 2)
abline(h = 10, col = "red", lty = 2)
dev.off()
