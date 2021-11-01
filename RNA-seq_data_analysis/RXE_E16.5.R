# Load up the relevant packages
# This script uses genes that have TPM > 1 as input
library(dplyr)
library(ggplot2)

# Load the TPMs file
TPMs <- read.csv("~/TPMs_Greater_than_1.csv", row.names=1)

# Need to log transform. Using pseudocount (+ 1) to avoid invalid log transformation on 0s
transformed <- log2(TPMs + 1)
colnames(transformed) <- c("Gene_ID","TPM_DEF.1","TPM_DEF.2","TPM_DEF.3","TPM_DEF.4","TPM_DEF.5","TPM_HET.1","TPM_HET.2","TPM_HET.3","TPM_HET.4","TPM_HET.5")

# Read in the reference table full of gene information and then combine with the log-transformed, filtered genes
ref <- read.csv("~/Mouse_Gene_Reference_GRCm38.csv", row.names=NULL)
ref$Gene_ID <- sub("\\..*", "", ref$Gene_ID)
head(ref)

transformed <- cbind(rownames(transformed), data.frame(transformed, row.names=NULL))
merged <- as.data.frame(merge(transformed, ref, by= "Gene_ID"))
head(merged)

# Separate out X chromosome genes and average the TPMs
Xchr <- filter(merged, merged$Chromosome == "chrX")
unique(Xchr$Chromosome)
Xvalues <- Xchr[,2:11]
head(Xvalues)
Xmeans <- as.data.frame(colMeans(Xvalues))
head(Xmeans)


# Separate out autosomal genes and average the TPMs
target <- c("chrX", "chrY", "chrM")
autosomes <- filter(merged, !merged$Chromosome %in% target)
autosomes <- autosomes %>% filter(!is.na(Chromosome))
unique(autosomes$Chromosome)
autosomes$Chromosome <- "autosome"
Avalues <- autosomes[,2:11]
Ameans <- as.data.frame(colMeans(Avalues))
head(Ameans)

# Calculate the RXE
RXE <- Xmeans - Ameans


# Time to make the plot!
# Fix the labels
RXELabels <- c("Def","Def","Def","Def","Def","Het","Het","Het","Het","Het")

newRXE <- cbind(RXELabels, rownames(RXE), data.frame(RXE, row.names=NULL))
colnames(newRXE) <- c("Genotype", "Sample", "RXE")

# Reorder so that the Defs aren't displayed first
level_order <- c('Het', 'Def')

p <- ggplot(newRXE, aes(x= factor(Genotype, level = level_order), y=RXE)) + 
  geom_boxplot() +
  xlab("Genotype")
p
p + geom_jitter(shape=16, position=position_jitter(0.2))
