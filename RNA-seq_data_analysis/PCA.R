# Load up the relevant packages
library(DESeq2)
library(ggplot2)

# Loading up the count table and metadata files
countData <- read.csv("Counts.csv", row.names = 1, header = TRUE, sep = ",")
head(countData)
groups <-read.csv("Metadata.csv",  header = T, sep = ",")
head(groups)

#Filtering out genes w/ below 5 counts avg:
keep <- countData[rowMeans(countData) >= 5,]

# Formatting the data into a dds
dds <- DESeqDataSetFromMatrix(keep, colData = Metadata, design = ~ Batch + Phenotype)

# This is the transformation that the PCA plot will be using
vsd <- vst(dds, blind=TRUE)

# Making the PCA plot
pcaData <- plotPCA(vsd, intgroup=c("Batch", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Batch)) +
  geom_point(size=4) +
  theme(legend.text=element_text(size=16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
