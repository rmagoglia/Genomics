#!/usr/bin/env Rscript

# Assesses mapping bias between two mappings of the same sample

args = commandArgs(TRUE)
fileA = args[1]
fileB = args[2]
pdfFile = args[3]

dataA = read.table(fileA, header=TRUE)
dataB = read.table(fileB, header=TRUE)

countsA_ref = dataA$ref_counts
countsA_alt = dataA$alt_counts
countsB_ref = dataB$ref_counts
countsB_alt = dataB$alt_counts

maxRef = max(c(countsA_ref, countsB_ref))
maxAlt = max(c(countsA_alt, countsB_alt))

spearRef = round(cor.test(countsA_ref, countsB_ref, method='spearman')$estimate, digits=4)
pearRef = round(cor.test(countsA_ref, countsB_ref, method='pearson')$estimate, digits=4)
  
spearAlt = round(cor.test(countsA_alt, countsB_alt, method='spearman')$estimate, digits=4)
pearAlt = round(cor.test(countsA_alt, countsB_alt, method='pearson')$estimate, digits=4)

pdf(pdfFile)
plot(countsA_ref, countsB_ref, main="Mapping Bias: Reference Counts", ylab=fileB, xlab=fileA)
abline(h=0,v=0,a=0,b=1)
legend(100, maxRef*.9, legend = paste("r = ", spearRef, '\n', "rho = ", pearRef, sep=""), bty="n")

plot(countsA_alt, countsB_alt, main="Mapping Bias: Alternative Counts", ylab=fileB, xlab=fileA)
abline(h=0,v=0,a=0,b=1)
legend(100, maxAlt*.9, legend = paste("r = ", spearAlt, '\n', "rho = ", pearAlt, sep=""), bty="n")
dev.off()

cat(paste("Correlations:\n\tReference counts (spearman) = ", spearRef, "\n\tReference counts (pearson) = ", pearRef, "\n\tAlternative counts (spearman) = ", spearAlt, "\n\tAlternative counts (pearson) = ", pearAlt, "\n", sep=""))

totalRef = sum(countsA_ref, na.rm=TRUE) + sum(countsB_ref, na.rm=TRUE)
totalAlt = sum(countsA_alt, na.rm=TRUE) + sum(countsB_alt, na.rm=TRUE)
total = totalRef + totalAlt
percentRef = round(((totalRef / total)*100), digits=4)
percentAlt = round(((totalAlt / total)*100), digits=4)

cat(paste("Reference counts = ", as.character(totalRef), " (", as.character(percentRef), "%)\n", sep=""))
cat(paste("Alternative counts = ", as.character(totalAlt), " (", as.character(percentAlt), "%)\n", sep=""))