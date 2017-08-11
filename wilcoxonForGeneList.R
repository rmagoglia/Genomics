#!/usr/bin/env Rscript

# Takes an ASE file and a list of genes and determines whether there is enrichment in either species for that subset of genes

args = commandArgs(TRUE)
aseFile = args[1]
geneFile = args[2]

data = read.table(aseFile, header=TRUE)
genes = read.table(geneFile, header=FALSE, col.names=c('gene'))

counts_ref = c()
counts_alt = c()
genesPassed = c()

# Only look at genes with 10+ counts in each mapping, change 0's to 1's to avoid undefined values
for(i in seq(1, length(data$ref_counts))){
  if((data$ref_counts[i] + data$alt_counts[i]) > 9){
    
    genesPassed = c(genesPassed, as.character(data$gene[i]))
    
    if(data$ref_counts[i] == 0){
      counts_ref = c(counts_ref, 1)
    } else {
      counts_ref = c(counts_ref, data$ref_counts[i])
    }
    if(data$alt_counts[i] == 0){
      counts_alt = c(counts_alt, 1)
    } else {
      counts_alt = c(counts_alt, data$alt_counts[i])
    }
  }
}

ase = log2(counts_ref / counts_alt)

backgroundAse = c()
signalAse = c()

for(i in seq(1, length(genesPassed))){
  if(genesPassed[i] %in% genes$gene){
    signalAse = c(signalAse, ase[i])
  } else {
    backgroundAse = c(backgroundAse, ase[i])
  }
}

if(is.null(signalAse) | is.null(backgroundAse)){
  cat("Not enough genes in either background or test set.\n")
} else {
  # Use a Wilcoxon rank test to determine if there is a bias
  sig = wilcox.test(backgroundAse, signalAse)$p.value
  
  backgroundMed = median(backgroundAse)
  backgroundMean = mean(backgroundAse)
  signalMed = median(signalAse)
  signalMean = mean(signalAse)
  
  cat(paste("Genes being tested: ", geneFile, "\n", sep=""))
  cat(paste("Sample being tested: ", aseFile, "\n\n", sep=""))
  cat(paste("Background genes' median ASE value: ", backgroundMed, "\n", sep=""))
  cat(paste("Signal genes' median ASE value: ", signalMed, "\n\n", sep=""))
  cat(paste("Background genes' mean ASE value: ", backgroundMean, "\n", sep=""))
  cat(paste("Signal genes' mean ASE value: ", signalMean, "\n\n", sep=""))
  cat(paste("Wilcoxon rank test p-value: ", sig, "\n", sep=""))
  
}

