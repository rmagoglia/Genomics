#!/usr/bin/env Rscript

# Takes an parental and hybrid ASE files and plots cis versus trans effects on expression

args = commandArgs(TRUE)
parentA = args[1]
parentB = args[2]
hybrid = args[3]
pdfFile = args[4]

# Read in files and initialize counts
dataA = read.table(parentA, header=TRUE)
dataB = read.table(parentB, header=TRUE)
dataHy = read.table(hybrid, header=TRUE)

countsA_ref = c()
countsB_alt = c()
countsHy_ref = c()
countsHy_alt = c()

# Only look at genes with 5+ counts in each parental sample, 10+ reads in hybrid, change 0's to 1's to avoid undefined values
for(i in seq(1, length(dataA$ref_counts))){
  if((dataA$ref_counts[i] > 4) & (dataB$alt_counts[i] > 4) & (dataHy$alt_counts[i] + dataHy$ref_counts[i]> 9)){
    if(dataA$ref_counts[i] == 0){
      countsA_ref = c(countsA_ref, 1)
    } else {
      countsA_ref = c(countsA_ref, dataA$ref_counts[i])
    }
    
    if(dataB$alt_counts[i] == 0){
      countsB_alt = c(countsB_alt, 1)
    } else {
      countsB_alt = c(countsB_alt, dataB$alt_counts[i])
    }
    
    if(dataHy$ref_counts[i] == 0){
      countsHy_ref = c(countsHy_ref, 1)
    } else {
      countsHy_ref = c(countsHy_ref, dataHy$ref_counts[i])
    }
    
    if(dataHy$alt_counts[i] == 0){
      countsHy_alt = c(countsHy_alt, 1)
    } else {
      countsHy_alt = c(countsHy_alt, dataHy$alt_counts[i])
    }
  }
}


# Need to normalize for library size (use total counts as a proxy)
totalA = sum(countsA_ref, na.rm=TRUE)
totalB = sum(countsB_alt, na.rm=TRUE)
totalHy = sum(countsHy_alt, na.rm=TRUE) + sum(countsHy_ref, na.rm=TRUE)

# Calculate correlations
spearRef = round(cor.test(countsA_ref, countsHy_ref, method='spearman')$estimate, digits=4)
pearRef = round(cor.test(countsA_ref, countsHy_ref, method='pearson')$estimate, digits=4)

spearAlt = round(cor.test(countsHy_alt, countsB_alt, method='spearman')$estimate, digits=4)
pearAlt = round(cor.test(countsHy_alt, countsB_alt, method='pearson')$estimate, digits=4)

# Get max values to scale plots
maxRef = max(countsA_ref/totalA)
maxAlt = max(countsB_alt/totalB)
minRef = min(countsHy_ref/totalHy)
minAlt = min(countsHy_alt/totalHy)


pdf(pdfFile)
plot(countsHy_alt*2/totalHy, countsB_alt/totalB, main="Tetraploid Expression: Species B", ylab="Expression in Parent B", xlab="Expression in Hybrid")
abline(h=0,v=0,a=0,b=1)
legend(minAlt*0.9, maxAlt*.9, legend = paste("r = ", spearAlt, '\n', "rho = ", pearAlt, sep=""), bty="n")

plot(countsHy_ref*2/totalHy, countsA_ref/totalA, main="Tetraploid Expression: Species A", ylab="Expression in Parent A", xlab="Expression in Hybrid")
abline(h=0,v=0,a=0,b=1)
legend(minRef*0.9, maxRef*.9, legend = paste("r = ", spearRef, '\n', "rho = ", pearRef, sep=""), bty="n")
dev.off()




