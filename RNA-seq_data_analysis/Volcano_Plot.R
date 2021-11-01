# Load up the relevant packages
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load in the DESeq2 results
results <- read.csv("results_Taf4b_Fertile_vs_Infertile.csv")

# Second, filter for genes with TPM > 1 average
TPMs <- read.csv("~/TPMs_Greater_than_1.csv")
results2 <- merge(results, TPMs, by="Gene_ID")

# Third, make sure we're specifying the DEGs 
cleaned <- filter(results2, results2$padj != "NA")
cleaned2 <- filter(cleaned, cleaned$Gene_Type == "protein_coding")
cleaned2$Significant <- ifelse(cleaned2$padj < 0.05, "p-adj. < 0.05", "Not Sig")

# Take a peek at the file to see what are the top 5 DEGs (based on p-adj)
view(cleaned2)

# Third, make the plot
plot <- ggplot(cleaned2, aes(x=log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("gray50", "red")) +
  theme(legend.position = "bottom") 

plot

# Fourth, add the labels (these are just the labels for the E16.5 RNA-seq DEGs)
target <- c("Trim52","Magea6","Pramel3","Sfpq","Nav3","Taf4b")

DEGs <- filter(cleaned2, Gene_Name %in% target)

plot + geom_text_repel(
  data = DEGs,
  aes(label = Gene_Name),
  size = 5,
  box.padding = unit(0.45, "lines"),
  point.padding = unit(0.4, "lines"),
  fontface = 2
)
