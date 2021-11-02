library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(monocle3)

## Load up data

zhao.data <- Read10X(data.dir = "~/scRNA-seq/Zhao/10x_data/Aggr_outs/count/filtered_feature_bc_matrix/")

zhao.data[c("Dazl", "Stra8", "Taf4b"), 1:30]

zhao <- CreateSeuratObject(counts = zhao.data, project = "germ_cells", min.cells = 3, min.features =  200) 
zhao

## Add cell metadata

Metadata_barcodes <- read.delim("~/scRNA-seq/Zhao/10x_data/Aggr_outs/Metadata_barcodes.tsv", row.names=1)

zhao <- AddMetaData(zhao, Metadata_barcodes)
head(zhao@meta.data)

# Continue pre-processing

zhao[["percent.mt"]] <- PercentageFeatureSet(zhao, pattern = "^mt-")
head(zhao@meta.data, 5)

VlnPlot(zhao, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(zhao, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zhao, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Subset the data, scale, and normalize

zhao <- subset(zhao, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 2500 & nCount_RNA < 30000 & percent.mt < 5)

VlnPlot(zhao, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## Do some clustering and filter for Dazl

Dazl_expression <- GetAssayData(object = zhao, assay = "RNA", slot = "data")["Dazl",]
pos_ids <- names(which(Dazl_expression>0))
pos_cells <- subset(zhao,cells=pos_ids)


#### Create a Monocle CDS Object

# Create an expression matrix
expression_matrix <- pos_cells@assays$RNA@counts

# Get cell metadata
cell_metadata <- pos_cells@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(pos_cells@assays$RNA), row.names = rownames(pos_cells@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
my.cds <- preprocess_cds(my.cds, num_dim = 100)
plot_pc_variance_explained(my.cds)


## Step 2: Reduce the dimensions using UMAP
my.cds <- reduce_dimension(my.cds, reduction_method = "UMAP")
plot_cells(my.cds)
plot_cells(my.cds, genes="Taf4b")
plot_cells(my.cds, genes=c("Dazl", "Stra8", "Taf4b", "Figla"))

## Step 3: Cluster the cells
my.cds <- cluster_cells(my.cds)
plot_cells(my.cds, color_cells_by="Time", group_cells_by="cluster",label_cell_groups = F,label_groups_by_cluster = F)

my.cds <- learn_graph(my.cds)

my.cds <- order_cells(my.cds)
plot_cells(my.cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

pr_graph_test_res <- graph_test(my.cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

my_genes <- c("Figla", "Stra8","Taf4","Taf4b")
my_germ_cds <- my.cds[rowData(my.cds)$gene_short_name %in% my_genes]

plot_genes_in_pseudotime(my_germ_cds,
                         color_cells_by="Pseudotime",
                         min_expr=0.5)

