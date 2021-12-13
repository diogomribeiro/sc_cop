#13-May-2021 Diogo Ribeiro @ UNIL
# Exploring overlap of proteomics data with COPs and non-COPs

library(data.table)
library(ggplot2)

options(scipen=1)
args = commandArgs(trailingOnly=TRUE)

##########
# Input parameters
##########

proteinFile = args[1] #"mirauta_2020_elife/elife-57390-fig1-data3-v3.xls"
proteinData = fread( proteinFile, stringsAsFactors = FALSE, header = T, sep="\t")

wantedGenes = proteinData$ensembl_gene_id

listDonorRun = fread( "../data/cuomo2021_list_individual_experiment.tsv", stringsAsFactors = FALSE, header = F, sep="\t")

## SINGLE CELL
copFile = "../data/cuomo2021_sc_cops_control.bed.gz"

## BULK
# copFile = "../data/cuomo2021_bulk_cops_control.bed.gz"

##########
# List of donors in common
##########
listDonorRun$donor = data.table(unlist(lapply(listDonorRun$V1, function(x) unlist(strsplit(x,"[-]"))[1])))$V1

proteomicsDonors = data.table(colnames(proteinData))

proteomicsDonors$donor = data.table(unlist(lapply(proteomicsDonors$V1, function(x) unlist(strsplit(x,"[-]"))[2])))$V1
proteomicsDonors$donor = data.table(unlist(lapply(proteomicsDonors$donor, function(x) unlist(strsplit(x,"[@]"))[1])))$V1

donorMerge = merge(listDonorRun,proteomicsDonors, by = "donor", all.x = T)

donorMerge = donorMerge[!is.na(V1.y)]

##########
# Process protein data (sum intensities for several isoforms of same gene)
##########
ncol(proteinData)

# Remove HPSI0314i-bubh_3 (control)
# Remove Strontium
dropCols = grep("Strontium", colnames(proteinData))
proteinData[, (dropCols) := NULL]
dropCols = grep("bubh", colnames(proteinData))
proteinData[, (dropCols) := NULL]
ncol(proteinData)

geneCol = grep("ensembl_gene_id", colnames(proteinData))
wantedCols = grep("HPSI", colnames(proteinData))

d = data.table(ensembl_gene_id = proteinData$ensembl_gene_id)
geneLabels = unique(d[order(ensembl_gene_id)])

finalDT = data.table()
for (c in wantedCols){
    co = c(c,geneCol)
    dt = proteinData[,..co]
    if (colnames(dt)[1] %in% donorMerge$V1.y){
      print(colnames(dt)[1])
    res = data.table(aggregate(dt[,1], list(intensity = dt$ensembl_gene_id), sum ))
    finalDT = cbind(finalDT, data.table(res[,2]))
  }
}
finalDT = cbind(finalDT,geneLabels)

geneNames = finalDT$ensembl_gene_id
finalDT$ensembl_gene_id = NULL

##########################
# Gene pair part
##########################

rowmeans = data.table(rowMeans(finalDT, na.rm = T))
rowmeans$gene = geneNames

copData = fread( copFile, stringsAsFactors = FALSE, header = T, sep="\t")
copData = copData[,.(pairID, significant, nullId)]

copData$gene1 = data.table(unlist(lapply(copData$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
copData$gene2 = data.table(unlist(lapply(copData$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

colnames(rowmeans) = c("intensity1","gene1")
d = merge(copData, rowmeans, by = "gene1", all.x =T)
colnames(rowmeans) = c("intensity2","gene2")
mergedData = merge(d, rowmeans, by = "gene2", all.x =T)

mergedData = mergedData[gene1 %in% geneNames][gene2 %in% geneNames]

mergedData$nullId = NULL
mergedData = unique(mergedData)

table(mergedData$significant)

cor1 = cor.test(mergedData[significant == 1]$intensity1, mergedData[significant == 1]$intensity2, method = "spearman")
cor2 = cor.test(mergedData[significant == 0]$intensity1, mergedData[significant == 0]$intensity2, method = "spearman")
text = paste("COP spearman R: ",round(cor1$estimate,2), "P-value:", format.pval(cor1$p.value,1), "N:",length(unique(mergedData[significant ==1]$pairID)),
             "\nnon-COP spearman R: ",round(cor2$estimate,2), "P-value:", format.pval(cor2$p.value,1), "N:",length(unique(mergedData[significant ==0]$pairID)))

ggplot( mergedData, aes(x = log10(intensity1), y = log10(intensity2), color = as.factor(significant), fill = as.factor(significant))) +
  geom_point(alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = text, vjust = 1.1, hjust = 1, size = 7, fontface = "bold"  ) +
  geom_smooth(method = "lm") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, 
        legend.title=element_blank(), legend.position = "None", legend.text = element_text(size = 22))

