
options(scipen=1)

library(ggplot2)
library(data.table)

peakAssociationFDRCutoff = 0.05
peakCorrCutoff = 0.05

### Read gene-enhancer interactions
## Rep3
inFile = "/scratch/dribeir1/single_cell/multi_omics/coex_peak_F0.5.tsv"
# inFile = "/scratch/dribeir1/single_cell/multi_omics/coex_peak_F0.5_coding_nofilt.tsv"
# inFile = "/scratch/dribeir1/single_cell/multi_omics/coex_peak_100kb.tsv"
## Rep2
# inFile = "/scratch/dribeir1/single_cell/multi_omics/rep2/coex_peak_F0.5.tsv"
# inFile = "/scratch/dribeir1/single_cell/multi_omics/rep2/coex_peak_100kb.tsv"

## Epimap
# inFile = "zcat /scratch/dribeir1/single_cell/multi_omics/vs_epimap/epimap_lcl.bed.gz"
# peakGeneData = fread( inFile, stringsAsFactors = FALSE, header = F, sep="\t")
# colnames(peakGeneData) = c("chr","start","end","gene","corr")

geneModels = fread( "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/gencode/gencode_v19/genes_tss.bed", stringsAsFactors = FALSE, header = F, sep="\t")
geneModels$gene = data.table(unlist(lapply(geneModels$V4, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
geneModels = geneModels[,.(V1,V2,gene)]
colnames(geneModels) = c("gene_chr","gene_tss","gene")

peakGeneData = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")
summary(peakGeneData$corr)

peakGeneData$gene = data.table(unlist(lapply(peakGeneData$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
peakGeneData$peak = data.table(unlist(lapply(peakGeneData$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

peakGeneData$fdr = p.adjust(peakGeneData$corrPval, method = "BH")
peakGeneData = peakGeneData[,.(gene,peak,corr,corrSign,corrPval,fdr)]

peakGeneData$chr = data.table(unlist(lapply(peakGeneData$peak, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
peakGeneData$start = as.numeric(data.table(unlist(lapply(peakGeneData$peak, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
peakGeneData$end = as.numeric(data.table(unlist(lapply(peakGeneData$peak, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)

peakGeneData$midpoint = (peakGeneData$start + peakGeneData$end) / 2

mergedData = merge(peakGeneData, geneModels, by = "gene", all.x = T)

mergedData$distance = abs(mergedData$midpoint - mergedData$gene_tss)

mergedData$group = "not significant"
# mergedData[corr > 0.1]$group = "significant"
mergedData[corr > peakCorrCutoff][fdr < peakAssociationFDRCutoff]$group = "significant"

ggplot( mergedData, aes(x = distance, fill = group ))  +
  geom_density(alpha = 0.5) +
  theme_linedraw() + 
  # xlim(c(0,1000000)) + # for epimap
  scale_fill_brewer(palette = "Set2") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)


nrow(mergedData[group == "significant"][distance > 100000]) / nrow(mergedData[group == "significant"])
