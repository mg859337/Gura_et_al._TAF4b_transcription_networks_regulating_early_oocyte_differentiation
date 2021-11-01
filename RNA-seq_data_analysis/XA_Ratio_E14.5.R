# The goal of this script is to determine the X:A ratio in E14.5 oocytes between Taf4b-wt and Taf4b-heterozygous versus Taf4b-deficient samples. 
# This procedure is used after determining which genes have an average TPM > 1

# Pull up the TPMs file
TPMs <- read.csv("~/TPMs_Greaterthan1.csv")


# Average the genotypes
TPMs$WtHetMeans <- rowMeans(subset(TPMs, select = c("TPM_F.HET.1_downsampled","TPM_F.HET.2_downsampled","TPM_F.HET.3_downsampled","TPM_F.HET.4_downsampled","TPM_F.WT.1","TPM_F.WT.2")))
TPMs$DefMeans <- rowMeans(subset(TPMs, select = c("TPM_F.DEF.1","TPM_F.DEF.2","TPM_F.DEF.3","TPM_F.DEF.4")))


# Set aside the averaged values into a new file, read in a reference table, and merge the two together
Filtered_Avg_TPMs <- cbind(TPMs[,1], TPMs[,12:13])
colnames(Filtered_Avg_TPMs) <- c("Gene_ID","WtHetMeans","DefMeans")

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

# Use pairwiseCI on the WT/HET data first
CI_WtHet <- pairwiseCI(WtHetMeans ~ Chromosome, data = reformatted ,alternative = "two.sided", conf.level = 0.95, method = "Median.ratio")
CI_WtHet

# Use pairwiseCI on the DEF data 
CI_Def <- pairwiseCI(DefMeans ~ Chromosome, data = reformatted ,alternative = "two.sided", conf.level = 0.95, method = "Median.ratio")
CI_Def


# Load up the packages needed to make a plot
library(plyr)
library(lattice)
library(reshape)
library(ggplot2)

# Now turn the output [not shown] into image
d <- data.frame(Time=c("Wt/Het", "Def"), 
                median=  c(0.7974, 0.8257),
                lower=   c(0.6752, 0.6734), 
                upper=   c(0.9619, 0.9769))

d$order <- factor(d$Time, levels = c("Wt/Het", "Def"), ordered = TRUE)

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
