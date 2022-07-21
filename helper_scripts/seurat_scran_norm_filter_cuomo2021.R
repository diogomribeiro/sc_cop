#16-Mar-2021 Diogo Ribeio @ UNIL
# Script to normalise and/or filter single cell RNA-seq entries using Seurat

library(data.table)
library(Seurat) 
library(ggplot2)
library(scater)
library(scran)
library(SingleCellExperiment)

args = commandArgs(trailingOnly=TRUE)

############
# Read input file
############

inFile = args[1] # "cuomo_matrix.bed.gz"
outFile = args[2]

raw_counts = fread(inFile, header = T, sep = "\t")
head(raw_counts[,1:10])

infoColumns = raw_counts[,1:6]

# # Get gene name from info (e.g. L=NA;T=protein_coding;R=chrchr1:134901-139379;N=AL627309.1)
# raw_counts$geneEntry = data.table(unlist(lapply(raw_counts$info, function(x) unlist(strsplit(x,";"))[4])))$V1
# genes = data.table(unlist(lapply(raw_counts$geneEntry, function(x) unlist(strsplit(x,"="))[2])))$V1
# raw_counts$geneEntry = NULL
genes = raw_counts$gene

# Use only gene and values
raw_counts = raw_counts[,7:ncol(raw_counts)]
raw_counts = cbind(genes,raw_counts)

###########
# Create Seurat object
###########
# Converting data.table into matrix (so that it is read by Seurat)
mat = as.matrix(raw_counts)
colnames(mat) = colnames(raw_counts)
rownames(mat) = raw_counts$gene
mat = mat[,2:ncol(mat)] # remove "gene" column
head(mat[,1:10])
rm(raw_counts) # this is to save computer RAM
gc() # garbage collection to save computer RAM

# Create seurat object
# Note: this command can also be used to prefilter things like the # cells per gene
cuomo2021 <- CreateSeuratObject(mat, min.cells = 1, min.genes = 1, project = "data", names.delim = 4)
gc()

# Check if all is fine
cuomo2021 # this should say: X features (i.e. genes) across X samples (i.e. cells) within 1 assays. If not, something went wrong 
cuomo2021@assays$RNA@counts # dots mean 0
head(cuomo2021@meta.data)
table(cuomo2021@meta.data$orig.ident) # be sure that there is only one item here (share-seq)

rm(mat) # to save RAM
gc()

FeatureScatter(cuomo2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(text = element_text(size=20))

###########
# Actually normalise and scale data
###########
cuomo2021.scran <- SingleCellExperiment( assays = list(counts = cuomo2021@assays$RNA@counts), rowData = rownames(cuomo2021))
clusters <- quickCluster(cuomo2021.scran)
cuomo2021.scran <- computeSumFactors(cuomo2021.scran, clusters=clusters)

cuomo2021.scran <- logNormCounts(cuomo2021.scran, pseudo.count = 1) # by default this is log2

gc()

###########
# Write matrix
###########

infoColumns$name = genes
dt = data.table(as.matrix(cuomo2021.scran@assays@data$logcounts))
gc()
dt$name = rownames(cuomo2021.scran@assays@data$logcounts)

mergedData = merge(infoColumns, dt, by = "name")
mergedData$name = NULL

rm(dt)
rm(infoColumns)
rm(cuomo2021)
gc()
mergedData = mergedData[do.call(order, mergedData[, c("#chr", "start", "end")]), ]

write.table(mergedData, file = outFile, quote = F, sep = "\t", row.names = F)
