#17-Mar-2021 Diogo Ribeiro @ UNIL
# Script to split matrix into certain donor-experiment combinations

library(data.table)
library(Seurat) 

############
# Read input
############

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("All arguments need to be provided", call.=FALSE)
}

inFile = args[1] #cuomo2021_scrna_scranN_regionF_codingF.bed.gz
medataFile = args[2] #CellSampleMetadata_reformat.txt
outFolder = args[3]

metaData = fread(medataFile, header = T, sep = "\t")

raw_counts = fread(inFile, header = T, sep = "\t")

head(raw_counts[,1:10])
infoColumns = raw_counts[,1:6]
genes = raw_counts$gene

# Use only gene and values
raw_counts = raw_counts[,7:ncol(raw_counts)]
raw_counts = cbind(genes,raw_counts)

# For per cell metadata
metaData = metaData[cell_name %in% colnames(raw_counts)]
rownames(metaData) = metaData$cell_name

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
cuomo2021 <- CreateSeuratObject(mat, min.cells = 1, min.genes = 1, project = "data", names.field = 2, meta.data = metaData)
cuomo2021
head(cuomo2021@meta.data)

rm(mat)
gc()

infoColumns$name = genes

### Split by donor-run
dt = data.table(table(cuomo2021@meta.data$donor_id.run_id))[order(N)]
cuomo2021.donor_run <- SplitObject(cuomo2021, split.by = "donor_id.run_id")

count = 0
## Write files
for (entry in cuomo2021.donor_run){
  count = count + 1
  m = data.table(as.matrix(entry@assays$RNA@counts))
  m$name = rownames(entry@assays$RNA@counts) # all datasets have same nrows
  mergedData = merge(infoColumns, m, by = "name")
  mergedData$name = NULL

  mergedData = mergedData[do.call(order, mergedData[, c("#chr", "start", "end")]), ]

  write.table(mergedData, file = paste(outFolder,"/",unique(entry@meta.data$donor_id.run_id),".bed",sep=""), quote = F, sep = "\t", row.names = F)
  
}
