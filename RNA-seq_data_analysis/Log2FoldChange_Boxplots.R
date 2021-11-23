# Load up the relevant packages
library(dplyr)
library(ggplot2)

# Load in the DESeq2 output
results <- read.csv("results_Taf4b_Fertile_vs_Infertile.csv")

# Consolidate into autosomes vs X chromosome

Xchr <- filter(results, results$Chromosome == "chrX")
unique(Xchr$Chromosome)

target <- c("chrX", "chrY", "chrM")
autosomes <- filter(results, !results$Chromosome %in% target)
autosomes <- autosomes %>% filter(!is.na(Chromosome))
unique(autosomes$Chromosome)
autosomes$Chromosome <- "autosomes"

reformatted <- rbind(autosomes, Xchr)

# Setting the level order so that R doesn't do it wrong
level_order <- c('autosomes','chrX')

# Create the plot
p <- ggplot(reformatted, aes(x= factor(Chromosome, level = level_order), y=log2FoldChange)) + 
 geom_boxplot(outlier.shape = NA) +
 coord_cartesian(ylim=c(-1.5, 1.5)) +
 xlab("Chromosome")
p
