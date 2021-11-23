# Load up the relevant packages
# This script uses genes that have TPM > 1 as input
library(dplyr)
library(ggplot2)

# Load the TPMs file
TPMs <- read.csv("~/TPMs_Greaterthan1.csv", row.names=1)

# Need to log transform. Using pseudocount (+ 1) to avoid invalid log transformation on 0s
transformed <- log2(TPMs + 1)

# Read in the reference table full of gene information and then combine with the log-transformed, filtered genes
ref <- read.csv("~/Mouse_Gene_Reference_GRCm38.csv", row.names=NULL)
ref$Gene_ID <- sub("\\..*", "", ref$Gene_ID)
head(ref)

transformed <- cbind(rownames(transformed), data.frame(transformed, row.names=NULL))
colnames(transformed) <- c("Gene_ID","TPM_F.DEF.1","TPM_F.DEF.2","TPM_F.DEF.3","TPM_F.DEF.4","TPM_F.HET.1_downsampled","TPM_F.HET.2_downsampled","TPM_F.HET.3_downsampled","TPM_F.HET.4_downsampled","TPM_F.WT.1","TPM_F.WT.2")
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
RXELabels <- c("Def","Def","Def","Def","Wt/Het","Wt/Het","Wt/Het","Wt/Het","Wt/Het","Wt/Het")

newRXE <- cbind(RXELabels, rownames(RXE), data.frame(RXE, row.names=NULL))
colnames(newRXE) <- c("Genotype", "Sample", "RXE")

# Reorder so that the Defs aren't displayed first
level_order <- c('Wt/Het', 'Def')

p <- ggplot(newRXE, aes(x= factor(Genotype, level = level_order), y=RXE)) + 
  geom_boxplot() +
  xlab("Genotype")
p
p + geom_jitter(shape=16, position=position_jitter(0.2))
