# The goal of this script is to determine the X:A ratio in E16.5 oocytes between Taf4b-heterozygous versus Taf4b-deficient samples. 
# This procedure is used after determining which genes have an average TPM > 1

# Pull up the TPMs file
TPMs <- read.csv("~/TPMs_Greater_than_1.csv")

# Average the genotypes
TPMs$HetMeans <- rowMeans(subset(TPMs, select = c("TPM_HET.1","TPM_HET.2","TPM_HET.3","TPM_HET.4","TPM_HET.5")))
TPMs$DefMeans <- rowMeans(subset(TPMs, select = c("TPM_DEF.1","TPM_DEF.2","TPM_DEF.3","TPM_DEF.4","TPM_DEF.5")))


# Set aside the averaged values into a new file, read in a reference table, and merge the two together
Filtered_Avg_TPMs <- cbind(TPMs[,1], TPMs[,12:13])
colnames(Filtered_Avg_TPMs) <- c("Gene_ID","HetMeans","DefMeans")

ref <- read.csv("~/Mouse_Gene_Reference_GRCm38.csv", row.names=NULL)
ref$Gene_ID <- sub("\\..*", "", ref$Gene_ID)
head(ref)

merged <- as.data.frame(merge(Filtered_Avg_TPMs, ref, by= "Gene_ID", all.x = T))
head(merged)

# Set aside the important information
All <- merged[,1:4]
head(All)


# Load up the relevant packages
library(dplyr)
library(pairwiseCI)

# Get the gene expression information about the X chromosome
Xchr <- filter(All, All$Chromosome == "chrX")
unique(Xchr$Chromosome)

# Get the gene expression information about the autosomes by removing the X, Y, and mitochondrial genes
target <- c("chrX", "chrY", "chrM")
autosomes <- filter(All, !All$Chromosome %in% target)
autosomes <- autosomes %>% filter(!is.na(Chromosome))
unique(autosomes$Chromosome)
autosomes$Chromosome <- "autosome"


# Put the X chromosome and autosomal data back together
reformatted <- rbind(autosomes, Xchr)
head(reformatted)

# Use pairwiseCI on the Het data first
CI_Het <- pairwiseCI(HetMeans ~ Chromosome, data = reformatted ,alternative = "two.sided", conf.level = 0.95, method = "Median.ratio")
CI_Het

# Use pairwiseCI on the Def data 
CI_Def <- pairwiseCI(DefMeans ~ Chromosome, data = reformatted ,alternative = "two.sided", conf.level = 0.95, method = "Median.ratio")
CI_Def


# Load up the packages needed to make a plot
library(plyr)
library(lattice)
library(reshape)
library(ggplot2)

# Now turn the output [not shown] into image
d <- data.frame(Time=c("Het", "Def"), 
                median=  c(0.9673, 0.8471),
                lower=   c(0.7798, 0.7091), 
                upper=   c(1.186, 1.018))
                
d$order <- factor(d$Time, levels = c("Het", "Def"), ordered = TRUE)

# Plot with connecting line
p <- ggplot(d, aes(x=order,y=median,group = 1)) + 
  scale_y_discrete(limit=c(0.5,1,1.5,2.0))+
  geom_line() + geom_hline(yintercept=1, color="light blue")+  
  geom_hline(yintercept=0.5, color="red")+
  geom_point() +
  geom_pointrange(data=d, mapping=aes(y=median, ymin=upper, ymax=lower),color='blue', size=1) + 
  labs(x = "Genotype", y = "X:A ratio")+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size = 14,face = "bold"))
p
