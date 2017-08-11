#!/usr/bin/env Rscript

# Takes parental expression files and a list of genes and determines whether there is enrichment in either species for that subset of genes

args = commandArgs(TRUE)
fileA = args[1]
fileB = args[2]
geneFile = args[3]

dataA = read.table(fileA, header=TRUE)
dataB = read.table(fileB, header=TRUE)
genes = read.table(geneFile, header=FALSE, col.names=c('gene'))

countsA_ref = c()
countsB_alt = c()
genesPassed = c()

# Only look at genes with 5+ counts in each sample, change 0's to 1's to avoid undefined values
for(i in seq(1, length(dataA$ref_counts))){
  if((dataA$ref_counts[i] > 4) & (dataB$alt_counts[i] > 4)){
    genesPassed = c(genesPassed, as.character(dataA$gene[i]))
    
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
  }
}

# Need to normalize for library size (use total counts as a proxy)
totalA = sum(countsA_ref, na.rm=TRUE)
totalB = sum(countsB_alt, na.rm=TRUE)

diffExp = log2((countsA_ref/totalA) / (countsB_alt/totalB))

background = c()
signal = c()

for(i in seq(1, length(genesPassed))){
  if(genesPassed[i] %in% genes$gene){
    signal = c(signal, diffExp[i])
  } else {
    background = c(background, diffExp[i])
  }
}

if(is.null(signal) | is.null(background)){
  cat("Not enough genes in one of background or test set.\n")
} else {
  # Use a Wilcoxon rank test to determine if there is a bias
  sig = wilcox.test(background, signal)$p.value
  
  backgroundMed = median(background)
  backgroundMean = mean(background)
  signalMed = median(signal)
  signalMean = mean(signal)
  
  cat(paste("Genes being tested: ", geneFile, "\n", sep=""))
  cat(paste("Samples being tested:", fileA, fileB, "\n\n", sep=" "))
  cat(paste("Background genes' median expression ratio: ", backgroundMed, "\n", sep=""))
  cat(paste("Signal genes' median expression ratio: ", signalMed, "\n\n", sep=""))
  cat(paste("Background genes' mean expression ratio: ", backgroundMean, "\n", sep=""))
  cat(paste("Signal genes' mean expression ratio: ", signalMean, "\n\n", sep=""))
  cat(paste("Wilcoxon rank test p-value: ", sig, "\n", sep=""))
  
}

