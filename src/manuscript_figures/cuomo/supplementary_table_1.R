
library(data.table)

inFile = "../data/cuomo2021_sc_cops_final_dataset.bed.gz"

data <- fread(inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

data$chr = data$`#chr`

data$gene1_name = data.table(unlist(lapply(data$centralInfo, function(x) unlist(strsplit(x,"[=]"))[5])))$V1
data$gene2_name = data.table(unlist(lapply(data$cisInfo, function(x) unlist(strsplit(x,"[=]"))[5])))$V1

data$gene1_tss = data$centralStart
data[centralStrand == "-"]$gene1_tss = data[centralStrand == "-"]$centralEnd

data$gene2_tss = data$cisStart
data[cisStrand == "-"]$gene2_tss = data[cisStrand == "-"]$cisEnd

tssData = unique(data[,.(chr, cisPheno, gene2_tss)])
d1 = unique(data[,.(chr, centralPhenotype,gene1_tss)])
colnames(d1) = c("chr","cisPheno","gene2_tss")
tssData = unique(rbind(tssData, d1 ))

geneNames = data[,.(centralPhenotype,gene1_name)]
colnames(geneNames) = c("cisPheno","gene2_name")
geneNames = unique(rbind(geneNames, data[,.(cisPheno,gene2_name)]))

tssData = merge(tssData,geneNames,by = "cisPheno")

pruned = data[,.(pairID, corr, dataset)]

resultDT = data.table()
for (id in unique(pruned$pairID)){
  dt = pruned[pairID == id]
  nExp = nrow(dt)
  meanCorr = mean(dt$corr)
  dt$ind = data.table(unlist(lapply(dt$dataset, function(x) unlist(strsplit(x,"[-]"))[1])))$V1

  resultDT = rbind(resultDT, data.table(pair_id = id, n_ind = length(unique(dt$ind)), n_exp = nrow(dt),
                                        mean_corr = round(mean(dt$corr),2), min_corr = round(min(dt$corr),2),
                                        max_corr = round(max(dt$corr),2), ind_exp = toString(dt$dataset)))
}

resultDT

resultDT$gene1 = data.table(unlist(lapply(resultDT$pair_id, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
resultDT$gene2 = data.table(unlist(lapply(resultDT$pair_id, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

mergedData = merge(resultDT,tssData, by.x = "gene2", by.y = "cisPheno")
colnames(tssData) = c("gene1","gene1_tss","gene1_name")
mergedData = merge(mergedData,tssData, by = "gene1",all.x = T)
