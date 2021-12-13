#12-Oct-2021 Diogo Ribeiro @ UNIL

library(data.table)

args = commandArgs(trailingOnly=TRUE)

inFile = args[1] # /cuomo2021_scrna_scranN_regionF_codingF.bed.gz
rdsFile = args[2]

data = fread(inFile, header = T, sep = "\t")
gbm = readRDS(rdsFile)

gbm@meta.data
table(gbm@meta.data$Phase)

dt = data.table(gbm@meta.data)

g2mcells = dt[Phase == "G2M"]$cell_name
scells = dt[Phase == "S"]$cell_name
g1cells = dt[Phase == "G1"]$cell_name

g2mData = data[, ..g2mcells]
ncol(g2mData)
nrow(g2mData)

g2mData = data[, ..g2mcells]
sData = data[, ..scells]
g1Data = data[, ..g1cells]

geneData = data[,1:6]
g2mData = cbind(geneData, g2mData)
g1Data = cbind(geneData, g1Data)
sData = cbind(geneData, sData)

write.table(g2mData,"cuomo2021_scrna_scranN_regionF_codingF_g2mCells.bed", quote=F,row.names=F,sep="\t")
write.table(sData,"cuomo2021_scrna_scranN_regionF_codingF_sCells.bed", quote=F,row.names=F,sep="\t")
write.table(g1Data,"cuomo2021_scrna_scranN_regionF_codingF_g1Cells.bed", quote=F,row.names=F,sep="\t")
