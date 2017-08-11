#!/usr/bin/env Rscript

# Plots ASE per gene in order of genomic position, each chromosome a different color; also tests sliding window bin of genes for significant deviation from the background (using a Wilcoxon rank test)

args = commandArgs(TRUE)
aseFile = args[1]
posFile = args[2]
pdfFile = args[3]
window = as.numeric(args[4])

libSize = as.numeric(scan(aseFile, what="character", n=4)[4])

data = read.table(aseFile, header=TRUE)
pos = read.table(posFile, header=FALSE, col.names=c('gene', 'chrom', 'position'))

dataPos = data[order(data$chrom,pos$position),]

ase = c()

for(i in seq(1, length(data$ref_counts))){
  if((dataPos$ref_counts[i] + dataPos$alt_counts[i]) > 9){
    
    if(dataPos$ref_counts[i] == 0){
      counts_ref = 1
    } else {
      counts_ref = dataPos$ref_counts[i]
    }
    if(dataPos$alt_counts[i] == 0){
      counts_alt = 1
    } else {
      counts_alt = dataPos$alt_counts[i]
    }
    
    ase = c(ase, log2(counts_ref / counts_alt))
  } else {
    ase = c(ase, NA)
  }
}

medians = c()
averages = c()
stats = c()

for(i in seq(1, (length(ase)-(window-1)))){
  values = ase[i:(i+(window-1))]
  background = c(ase[0:(i-1)], ase[i+window:length(ase)])
  stats = c(stats, tryCatch(wilcox.test(background, values)$p.value, error=function(e) NA))
  averages = c(averages, mean(values, na.rm=TRUE))
  medians = c(medians, median(values, na.rm=TRUE))
}

bonf = -log10(0.05/length(ase))

# Use the default ase values from Peter's script since we only need an overview of the trends
pdf(pdfFile)
par(mfrow=c(4,1))
plot(ase, col=dataPos$chrom, main=paste("ASE by Genomic Position for", aseFile), ylab="Log2(ASE Ratio)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
abline(h=0)
plot(averages,col=dataPos$chrom, main=paste("Average ASE in sliding window of", window, "genes"), ylab="Average Log2(ASE Ratio)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
abline(h=0)
plot(medians,col=dataPos$chrom, main=paste("Median ASE in sliding window of", window, "genes"), ylab="Median Log2(ASE Ratio)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
abline(h=0)
plot(-log10(stats), col=dataPos$chrom, main=paste("Significant deviation from background in sliding window of", window, "genes\n(Wilcoxon Rank Sum Test)"), ylab="-Log10(p-value)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
abline(h=bonf, lty=2)


for(i in levels(dataPos$chrom)){
  par(mfrow=c(4,1))
  plot(ase[which(dataPos$chrom==i)], col=dataPos$chrom, main=paste("Chromosome", i, ": ASE Ratios"), ylab="Log2(ASE Ratio)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=0)
  plot(log2(dataPos$ref_counts[which(dataPos$chrom==i)]/libSize), col=dataPos$chrom, main=paste("Chromosome", i, ": Reference reads"), ylab="Log2(Normalized Read Counts)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=mean(log2(dataPos$ref_counts[which(dataPos$chrom==i)]/libSize)),na.rm=TRUE)
  plot(log2(dataPos$alt_counts[which(dataPos$chrom==i)]/libSize), col=dataPos$chrom, main=paste("Chromosome", i, ": Alternate reads"), ylab="Log2(Normalized Read Counts)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=mean(log2(dataPos$alt_counts[which(dataPos$chrom==i)]/libSize)),na.rm=TRUE)
  
  top = bonf
  if(max(-log10(stats[which(dataPos$chrom==i)]), na.rm=TRUE)> bonf){
    top=max(-log10(stats[which(dataPos$chrom==i)]))
  }
  
  plot(-log10(stats[which(dataPos$chrom==i)]), col=dataPos$chrom, main=paste("Significant deviation from background in sliding window of", window, "genes\n(Wilcoxon Rank Sum Test)"), ylab="-Log10(p-value)", xaxt='n', xlab="Genomic Position", ylim=c(0, top+1), pch=20,cex=.5)
  abline(h=bonf, lty=2)
}

dev.off()
