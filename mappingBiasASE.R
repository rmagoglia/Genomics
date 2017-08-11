#!/usr/bin/env Rscript

# Assesses mapping bias between two mappings of the same sample
# Hybrid specific (plot ASE rather than counts)

args = commandArgs(TRUE)
fileA = args[1]
fileB = args[2]
pdfFile = args[3]

dataA = read.table(fileA, header=TRUE)
dataB = read.table(fileB, header=TRUE)

countsA_ref = c()
countsA_alt = c()
countsB_ref = c()
countsB_alt = c()

# Only look at genes with 10+ counts in each mapping, change 0's to 1's to avoid undefined values
for(i in seq(1, length(dataA$ref_counts))){
  if(((dataA$ref_counts[i] + dataA$alt_counts[i]) > 9) & (dataB$ref_counts[i] + dataB$alt_counts[i]) > 9){
    if(dataA$ref_counts[i] == 0){
      countsA_ref = c(countsA_ref, 1)
    } else {
      countsA_ref = c(countsA_ref, dataA$ref_counts[i])
    }
    if(dataA$alt_counts[i] == 0){
      countsA_alt = c(countsA_alt, 1)
    } else {
      countsA_alt = c(countsA_alt, dataA$alt_counts[i])
    }
    if(dataB$ref_counts[i] == 0){
      countsB_ref = c(countsB_ref, 1)
    } else {
      countsB_ref = c(countsB_ref, dataB$ref_counts[i])
    }
    if(dataB$alt_counts[i] == 0){
      countsB_alt = c(countsB_alt, 1)
    } else {
      countsB_alt = c(countsB_alt, dataB$alt_counts[i])
    }
    
  }
}


aseA = log2(countsA_ref / countsA_alt)
aseB = log2(countsB_ref / countsB_alt)

maxASE = max(aseB)
minASE = min(aseA)

spear = round(cor.test(aseA, aseB, method='spearman')$estimate, digits=4)
pear = round(cor.test(aseA, aseB, method='pearson')$estimate, digits=4)

pdf(pdfFile)
plot(aseA, aseB, main="Mapping Bias: Log2(ASE Ratio)", ylab=fileB, xlab=fileA)
abline(h=0,v=0,a=0,b=1)
legend(minASE*.9, maxASE*.9, legend = paste("r = ", spear, '\n', "rho = ", pear, sep=""), bty="n")

hist(aseA, main="ASE Distribution", xlab=fileA)
hist(aseB, main = "ASE Distribution", xlab = fileB)
dev.off()

cat(paste("Correlations:\n\tSpearman = ", spear, "\n\tPearson = ", pear, "\n", sep=""))

totalRef = sum(countsA_ref, na.rm=TRUE) + sum(countsB_ref, na.rm=TRUE)
totalAlt = sum(countsA_alt, na.rm=TRUE) + sum(countsB_alt, na.rm=TRUE)
total = totalRef + totalAlt
percentRef = round(((totalRef / total)*100), digits=4)
percentAlt = round(((totalAlt / total)*100), digits=4)

cat(paste("Reference counts = ", as.character(totalRef), " (", as.character(percentRef), "%)\n", sep=""))
cat(paste("Alternative counts = ", as.character(totalAlt), " (", as.character(percentAlt), "%)\n", sep=""))

refBiasedA = length(which(aseA > 0))
refBiasedPercentA = round(refBiasedA / length(aseA), digits=4)*100
altBiasedA = length(which(aseA < 0))
altBiasedPercentA = round(altBiasedA / length(aseA), digits=4)*100

cat(paste("Reference biased (reference mapped)= ", as.character(refBiasedA), " (", as.character(refBiasedPercentA), "%)\n", sep=""))
cat(paste("Alternative biased (reference mapped)= ", as.character(altBiasedA), " (", as.character(altBiasedPercentA), "%)\n", sep=""))


refBiasedB = length(which(aseB > 0))
refBiasedPercentB = round(refBiasedB / length(aseB), digits=4)*100
altBiasedB = length(which(aseB < 0))
altBiasedPercentB = round(altBiasedB / length(aseB), digits=4)*100

cat(paste("Reference biased (alternative mapped)= ", as.character(refBiasedB), " (", as.character(refBiasedPercentB), "%)\n", sep=""))
cat(paste("Alternative biased (alternative mapped)= ", as.character(altBiasedB), " (", as.character(altBiasedPercentB), "%)\n", sep=""))
