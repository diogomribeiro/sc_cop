
library(data.table)

inFile = "~/Downloads/GSE118723_scqtl-counts.txt"
data = fread(inFile, sep = "\t")

genesList = fread("~/Downloads/GSE118723_genes-pass-filter.txt.gz")
table(genesList$V2)

cellsList = fread("~/Downloads/GSE118723_quality-single-cells.txt.gz")
table(cellsList$V2)

nrow(data)
data = data[gene %in% genesList[V2 == TRUE]$V1]

dt = data.table(t(data))

cols = data.table(colnames(data))

cols$sample = data.table(unlist(lapply(cols$V1, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
cols$cellp1 = data.table(unlist(lapply(cols$V1, function(x) unlist(strsplit(x,"[.]"))[2])))$V1
cols$cellp2 = data.table(unlist(lapply(cols$V1, function(x) unlist(strsplit(x,"[.]"))[3])))$V1
cols$cell = paste(cols$cellp1,cols$cellp2,sep="-")

cols = rbind(data.table(V1 = "gene",sample = "NA",cellp1 = "NA", cellp2 = "NA",cell = "NA"),cols[cell %in% cellsList[V2 == TRUE]$V1])

filteredData = data[,colnames(data) %in% cols$V1, with = FALSE]

ncol(filteredData)
nrow(filteredData)

filteredData[1:10,1:9]

# create a matrix for each individual, on which to run CODer, see how this was done for Cuomo dataset!

for (s in unique(cols$sample)){
  if (s != "NA"){
    cells = c("gene",cols[sample == s]$V1)
    d = filteredData[,colnames(filteredData) %in% cells, with = FALSE]
    print(s)
    print(dim(d))
    write.table(d,paste0("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/sarkar2019/matrices/",s,".bed"),row.names=F,quote=F,sep="\t")
    
  }
}

write.table(cols,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/sarkar2019/sarkar2019_metadata.txt",quote=F,row.names=F,sep="\t")

