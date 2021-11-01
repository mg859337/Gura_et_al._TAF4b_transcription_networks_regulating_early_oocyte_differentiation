# Get the information about the sample gtfs (files that StringTie provides as output)
FileNames <- list.files(path="~/gtfs/")
name <- gsub(".gtf","",FileNames)

# Create a folder called "TPM_Calc"
# A gtf organizing loop
j <- 1
m <- 1+length(FileNames)
while(j<m){
  temp <- read.delim(list.files()[j], header = F, row.names = NULL, comment.char="#", stringsAsFactors=FALSE)
  temp2 <- subset(temp, V3 != "exon")
  temp3 <- read.table(text = as.character(temp2[,9]), sep = ";")
  
  geneid <- as.data.frame(gsub(".* ","",temp3[,1]))
  transid <- as.data.frame(gsub(".* ","",temp3[,2]))
  refgen <- as.data.frame(gsub(".* ","",temp3[,3]))
  cov <- as.data.frame(gsub(".* ","",temp3[,4]))
  fpkm <- as.data.frame(gsub(".* ","",temp3[,5]))
  tpm <- as.data.frame(gsub(".* ","",temp3[,6]))
  
  temp4 <- cbind(geneid, transid, refgen, cov, fpkm, tpm)
  
  colnames(temp4) <- c("Gene_ID", "Transcript_ID", "Ref_Gene_Name", "Coverage", "FPKM","TPM")
  
  write.csv(temp4, file = paste("~/TPM_Calc/", name[j], "_organized.csv", sep = ""), row.names = F)
  j <- j + 1
}

# Merging the transcripts into just gene names based on Ensembl ID

Files <- list.files(path="~/TPM_Calc/")
First <- read.csv(file=paste("~/TPM_Calc/", Files[1], sep=""), header=T)
FirstFix <- aggregate(cbind(Coverage, FPKM, TPM) ~ Ref_Gene_Name + Gene_ID, data = First, FUN = sum)

Second <- read.csv(file=paste("~/TPM_Calc/", Files[2], sep=""), header=T)
SecondFix <- aggregate(cbind(Coverage, FPKM, TPM) ~ Ref_Gene_Name + Gene_ID, data = Second, FUN = sum)

dataMerge <- merge(FirstFix, SecondFix, by=c("Gene_ID","Ref_Gene_Name"), all=T)

# Fixing the column names issue
x <- colnames(dataMerge)
fixednames <- gsub(".x", paste("_", name[1]),x)
fixednames <- gsub(".y", paste("_", name[2]),fixednames)
fixednames <- gsub(" ", "", fixednames)

colnames(dataMerge) <- fixednames

# Same steps as before but in a loop
for(i in 3:length(Files)){ 
  ReadInMerge <- read.csv(file=paste("~/TPM_Calc/", Files[i], sep=""),
                          header=T)
  FixInMerge <- aggregate(cbind(Coverage, FPKM, TPM) ~ Ref_Gene_Name + Gene_ID, data = ReadInMerge, FUN = sum)
  dataMerge <- merge(dataMerge, FixInMerge, by=c("Gene_ID","Ref_Gene_Name"), all = T)
  
  fixednames <- colnames(dataMerge)
  fixednames <- gsub(".x", paste("_", name[i-1]),fixednames)
  fixednames <- gsub(".y", paste("_", name[i]),fixednames)
  fixednames <- gsub(" ", "", fixednames)
  colnames(dataMerge) <- fixednames
  
}

write.csv(dataMerge, "All_TPM_Maths.csv", row.names = F)

solo <- dataMerge[,names(dataMerge) %in% colnames(dataMerge)[grepl("TPM",colnames(dataMerge))]]
final <- cbind(dataMerge[,2:1], solo)

write.csv(final, "All_TPMs.csv", row.names = F)
