#!/usr/bin/env Rscript

# Plots differential expression per gene in order of genomic position, each chromosome a different color; also tests sliding window bin of genes for significant deviation from the background (using a Wilcoxon rank test)

args = commandArgs(TRUE)
fileA = args[1]
fileB = args[2]
posFile = args[3]
pdfFile = args[4]
window = as.numeric(args[5])

dataA = read.table(fileA, header=TRUE)
dataB = read.table(fileB, header=TRUE)
pos = read.table(posFile, header=FALSE, col.names=c('gene', 'chrom', 'position'))

dataAPos = dataA[order(dataA$chrom,pos$position),]
dataBPos = dataB[order(dataA$chrom,pos$position),]

libSizeA = as.numeric(scan(fileA, what="character", n=4)[4])
libSizeB = as.numeric(scan(fileB, what="character", n=4)[4])
  
dex = c()

# Only look at genes with 5+ counts in each sample, change 0's to 1's to avoid undefined values
for(i in seq(1, length(dataA$ref_counts))){
  if((dataAPos$ref_counts[i] > 4) & (dataBPos$alt_counts[i] > 4)){
    if(dataAPos$ref_counts[i] == 0){
      countsA_ref = 1
    } else {
      countsA_ref = dataAPos$ref_counts[i]
    }
    
    if(dataBPos$alt_counts[i] == 0){
      countsB_alt = 1
    } else {
      countsB_alt = dataBPos$alt_counts[i]
    }
    dex = c(dex, (log2((countsA_ref/libSizeA) / (countsB_alt/libSizeB))))
  } else {
    dex = c(dex, NA)
  }
}

print(length(dex))
print(length(dataBPos$alt_counts))

mediansd = c()
averagesd = c()
statsd = c()

for(i in seq(1, (length(dex)-(window-1)))){
  values = dex[i:(i+(window-1))]
  background = c(dex[0:(i-1)], dex[i+window:length(dex)])
  statsd = c(statsd, tryCatch(wilcox.test(background, values)$p.value, error=function(e) NA))
  averagesd = c(averagesd, mean(values, na.rm=TRUE))
  mediansd = c(mediansd, median(values, na.rm=TRUE))
}

bonf = -log10(0.05/length(dex))

# Use the default dex values from Peter's script since we only need an overview of the trends
pdf(pdfFile)

par(mfrow=c(4,1))
plot(dex, col=dataAPos$chrom, main=paste("Differential expression by Genomic Position for\n", fileA, "&", fileB), ylab="Log2(Expression Ratio)", xaxt='n', xlab="Genomic Position")
abline(h=0)
plot(averagesd,col=dataAPos$chrom, main=paste("Average differential expression in sliding window of", window, "genes"), ylab="Average Log2(Expression Ratio)", xaxt='n', xlab="Genomic Position")
abline(h=0)
plot(mediansd,col=dataAPos$chrom, main=paste("Median differential expression in sliding window of", window, "genes"), ylab="Median Log2(Expression Ratio)", xaxt='n', xlab="Genomic Position")
abline(h=0)
plot(-log10(statsd), col=dataAPos$chrom, main=paste("Significant deviation from background in sliding window of", window, "genes\n(Wilcoxon Rank Sum Test)"), ylab="-Log10(p-value)", xaxt='n', xlab="Genomic Position")
abline(h=bonf, lty=2)

for(i in levels(dataAPos$chrom)){
  par(mfrow=c(4,1))
  plot(dex[which(dataAPos$chrom==i)], col=dataAPos$chrom, main=paste("Chromosome", i, ": Expression Ratios"), ylab="Log2(Expression Ratio)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=0)
  plot(log2(dataAPos$ref_counts[which(dataAPos$chrom==i)]/libSizeA), col=dataAPos$chrom, main=paste("Chromosome", i, ": Parent A Expression"), ylab="Log2(Normalized Read Counts)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=mean(log2(dataAPos$ref_counts[which(dataAPos$chrom==i)]/libSizeA)))
  plot(log2(dataBPos$alt_counts[which(dataAPos$chrom==i)]/libSizeB), col=dataAPos$chrom, main=paste("Chromosome", i, ": Parent B Expression"), ylab="Log2(Normalized Read Counts)", xaxt='n', xlab="Genomic Position", pch=20,cex=.5)
  abline(h=mean(log2(dataBPos$alt_counts[which(dataAPos$chrom==i)]/libSizeB)))
  
  top = bonf
  if(max(-log10(statsd[which(dataAPos$chrom==i)]), na.rm=TRUE)> bonf){
    top=max(-log10(statsd[which(dataAPos$chrom==i)]), na.rm=TRUE)
  }
  
  plot(-log10(statsd[which(dataAPos$chrom==i)]), col=dataAPos$chrom, main=paste("Significant deviation from background in sliding window of", window, "genes\n(Wilcoxon Rank Sum Test)"), ylab="-Log10(p-value)", xaxt='n', xlab="Genomic Position", ylim=c(0, top+1), pch=20,cex=.5)
  abline(h=bonf, lty=2)
}


dev.off()
