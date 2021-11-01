# Load up DESeq2
library(DESeq2)

# Read in the counts and metadata files
countData <- read.delim("Counts.csv", row.names = 1, header = TRUE, sep = ",")
head(countData)
groups <-read.csv("Metadata.csv", row.names = NULL, header = T, sep = ",")
head(groups)

#Filtering out genes w/ below 5 counts (total):
keep <- countData[rowMeans(countData) >= 5,]

# Create the dds object and run DESeq2, "Batch" being the mouse collection date and "Phenotype" distinguishing the fertile (WT & HET) and the infertile (DEF) samples
dds <- DESeqDataSetFromMatrix(keep, colData = groups, design = ~ Batch + Phenotype)
dds <- DESeq(dds)
res <- results(dds) 

#Make the res into a data frame:
resdata<- as.data.frame(res)

#Convert row names into first column:
resdata <- cbind(rownames(resdata), data.frame(resdata, row.names=NULL))

#Assigning gene names reference file (in OSCAR) to human or mouse ref:
ref <- read.csv("~/Mouse_Gene_Reference_GRCm38.csv", row.names=NULL)
ref$Gene_ID <- sub("\\..*", "", ref$Gene_ID)
head(ref)

#Merging resdata and ref:
merged <- as.data.frame(merge(resdata, ref, by.x = "rownames(resdata)",by.y = "Gene_ID", all.x = T))

# Rearranging the data to my liking:
merged2 <- merged[, c(8:17, 1:7)]
colnames(merged2) <- c("Chromosome", "Database", "ID_Type", "Start", "End","Strand","Evidence_Level", "Gene_Type", "Gene_Name", "MGI_ID", "Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# Save the results:
write.csv(merged2, row.names= FALSE, file = "results_Taf4b_Fertile_vs_Infertile.csv")
