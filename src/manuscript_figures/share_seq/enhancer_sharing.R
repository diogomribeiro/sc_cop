#07-May-2021
# Analysis of results of linking genes to enhancers in SHARE-seq

options(scipen=1)

library(ggplot2)
library(data.table)

peakAssociationFDRCutoff = 0.05
peakCorrCutoff = 0.05

### Read gene-enhancer interactions
inFile = "../data/share_seq_gene_enhancer_associations.tsv.gz"

peakGeneData = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")

summary(peakGeneData$corr)

peakGeneData$gene = data.table(unlist(lapply(peakGeneData$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
peakGeneData$peak = data.table(unlist(lapply(peakGeneData$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

peakGeneData$fdr = p.adjust(peakGeneData$corrPval, method = "BH")
peakGeneData = peakGeneData[,.(gene,peak,corr,corrSign,corrPval,fdr)]

peakGeneData[corr > peakCorrCutoff][fdr < peakAssociationFDRCutoff]

### Read gene-gene co-expression
inFile = "../data/shareseq_sc_cops_control.bed.gz"

geneGeneData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

cops = geneGeneData[significant == 1]
cops$gene1 = data.table(unlist(lapply(cops$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
cops$gene2 = data.table(unlist(lapply(cops$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
cops = cops[,.(pairID, gene1, gene2, corr, distance)]

noncops = geneGeneData[significant == 0]
noncops$gene1 = data.table(unlist(lapply(noncops$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
noncops$gene2 = data.table(unlist(lapply(noncops$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
noncops = noncops[,.(pairID, gene1, gene2, corr, distance)]

## Merge peak-gene and gene-gene data
colnames(peakGeneData) = c("gene1","peak1","peak_corr1","peak_corrSign1","peak_corrPval1","peak_corrFDR1")
copsMerge1 = merge(cops, peakGeneData, by = "gene1", allow.cartesian=TRUE)
colnames(peakGeneData) = c("gene2","peak2","peak_corr2","peak_corrSign2","peak_corrPval2","peak_corrFDR2")
copsMerge = merge(copsMerge1, peakGeneData, by = "gene2", allow.cartesian = TRUE)
copsMerge = unique(copsMerge[peak1 == peak2])

colnames(peakGeneData) = c("gene1","peak1","peak_corr1","peak_corrSign1","peak_corrPval1","peak_corrFDR1")
noncopsMerge1 = merge(noncops, peakGeneData, by = "gene1", allow.cartesian=TRUE)
colnames(peakGeneData) = c("gene2","peak2","peak_corr2","peak_corrSign2","peak_corrPval2","peak_corrFDR2")
noncopsMerge = merge(noncopsMerge1, peakGeneData, by = "gene2", allow.cartesian = TRUE)
noncopsMerge = unique(noncopsMerge[peak1 == peak2])

copsMerge[peak_corrSign1 == "-"]$peak_corr1 = -copsMerge[peak_corrSign1 == "-"]$peak_corr1
copsMerge[peak_corrSign2 == "-"]$peak_corr2 = -copsMerge[peak_corrSign2 == "-"]$peak_corr2
noncopsMerge[peak_corrSign1 == "-"]$peak_corr1 = -noncopsMerge[peak_corrSign1 == "-"]$peak_corr1
noncopsMerge[peak_corrSign2 == "-"]$peak_corr2 = -noncopsMerge[peak_corrSign2 == "-"]$peak_corr2

## Test whether the same enhancer is correlated with both genes on the pair

cor.test(copsMerge$peak_corr1, copsMerge$peak_corr2)
cor.test(noncopsMerge$peak_corr1, noncopsMerge$peak_corr2)

length(unique(copsMerge$pairID))
length(unique(noncopsMerge$pairID))

## Total number of peaks
copN = data.table(table(copsMerge$pairID))
noncopN = data.table(table(noncopsMerge$pairID))

copN$dataset = "COP"
noncopN$dataset = "Non-COP"
m = rbind(copN, noncopN)
test = wilcox.test(copN$N, noncopN$N)
text1 = paste("COP mean:", round(mean(m[dataset == "COP"]$N),1) )
text2 = paste("Non-COP mean:", round(mean(m[dataset == "Non-COP"]$N),1) )
text3 = paste("Wilcoxon p-value", format.pval(test$p.value,1 ) )
ggplot( rbind(copN,noncopN), aes(x = N, y = dataset, fill = dataset ))  +
  geom_boxplot( width = 0.3, size = 1.5,outlier.shape = NA) +
  annotate("text", x = Inf, y = Inf, label = text1, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = text2, hjust = 1.05, vjust = 3.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = text3, hjust = 1.05, vjust = 5.5, size = 7, fontface = "bold"  ) +
  ggtitle("# tested enhancers per gene pair") +
  coord_flip() +
  xlab("# enhancers") +
  scale_fill_manual( values = c("#66C2A5","#969696") ) +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)


totCOP = length(unique(copsMerge$pairID))
totNonCOP = length(unique(noncopsMerge$pairID))

## Total number of associated peaks
copN = data.table(table(copsMerge[peak_corr1 > peakCorrCutoff][peak_corrFDR1 < peakAssociationFDRCutoff][peak_corr2 > peakCorrCutoff][peak_corrFDR2 < peakAssociationFDRCutoff]$pairID))
noncopN = data.table(table(noncopsMerge[peak_corr1 > peakCorrCutoff][peak_corrFDR1 < peakAssociationFDRCutoff][peak_corr2 > peakCorrCutoff][peak_corrFDR2 < peakAssociationFDRCutoff]$pairID))
test = wilcox.test(copN$N, noncopN$N)

copN = rbind(copN, data.table(V1 = rep("other",totCOP-nrow(copN)),N = 0))
noncopN = rbind(noncopN, data.table(V1 = rep("other",totNonCOP-nrow(noncopN)),N = 0))

copN$dataset = "COP"
noncopN$dataset = "Non-COP"
mergedData = rbind(copN,noncopN)

summary(mergedData[dataset == "COP"]$N)

text1 = paste("COP mean:", round(mean(mergedData[dataset == "COP"]$N),1) )
text2 = paste("Non-COP mean:", round(mean(mergedData[dataset == "Non-COP"]$N),1) )
text3 = paste("Wilcoxon p-value", format.pval(test$p.value,1, eps=-Inf ) )
ggplot( mergedData, aes(x = N, y = dataset, fill = dataset ))  +
  # geom_violin() +
  geom_boxplot( width = 0.5, size = 1.5) +
  # geom_jitter( size = 0.5, alpha = 0.5, height = 0.1, width = 0.5) +
  annotate("text", x = round(mean(mergedData[dataset == "COP"]$N),1), y = "COP", label = round(mean(mergedData[dataset == "COP"]$N),1), hjust = -2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  annotate("text", x = round(mean(mergedData[dataset == "Non-COP"]$N),1), y = "Non-COP", label = round(mean(mergedData[dataset == "Non-COP"]$N),1), hjust = -2, vjust = 1.5, size = 8, fontface = "bold"  ) +
  # annotate("text", x = Inf, y = Inf, label = text3, hjust = 1.05, vjust = 5.5, size = 7, fontface = "bold"  ) +
  # ggtitle("# associated enhancers per gene pair") +
  coord_flip() +
  xlab("Number of enhancers") +
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual( values = c("#66C2A5","#969696") ) +
  theme_linedraw() + 
  theme(text = element_text(size=24), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)
text3

## % gene pairs with strong association to same enhancer
copN = data.table(table(copsMerge[peak_corr1 > peakCorrCutoff][peak_corrFDR1 < peakAssociationFDRCutoff][peak_corr2 > peakCorrCutoff][peak_corrFDR2 < peakAssociationFDRCutoff]$pairID))
noncopN = data.table(table(noncopsMerge[peak_corr1 > peakCorrCutoff][peak_corrFDR1 < peakAssociationFDRCutoff][peak_corr2 > peakCorrCutoff][peak_corrFDR2 < peakAssociationFDRCutoff]$pairID))

dt = c(nrow(copN),totCOP-nrow(copN),nrow(noncopN),totNonCOP-nrow(noncopN))
test = fisher.test(matrix(c(nrow(copN),totCOP,nrow(noncopN),totNonCOP),nrow=2))

mergedData = data.table( N = dt, group = c("COP","COP","Non-COP","Non-COP"), fill = c(3,2,1,0), sharing = c(1,0,1,0))
perc1 = mergedData[group == "COP"][sharing == 1]$N * 100.0 / (mergedData[group == "COP"][sharing == 1]$N + mergedData[group == "COP"][sharing == 0]$N)
perc2 = mergedData[group == "Non-COP"][sharing == 1]$N * 100.0 / (mergedData[group == "Non-COP"][sharing == 1]$N + mergedData[group == "Non-COP"][sharing == 0]$N)
text = paste("Fisher test OR:",round(test$estimate,1),"p-value:",format.pval(test$p.value,1,eps = -Inf))
ggplot( mergedData, aes(x = group, y = N, fill = as.factor(fill) ))  +
  geom_bar(stat = "identity", color = "black", size = 1) +
  scale_fill_manual( values = c("#d9d9d9","#969696","#ccece6","#66C2A5"), label = c("Not Sharing","Sharing", "Not Sharing", "Sharing")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "COP", y = mergedData[group=="COP"][sharing == 1]$N/2, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "COP", y = mergedData[group=="COP"][sharing == 0]$N + mergedData[group=="COP"][sharing == 1]$N/1.02, size = 8) +
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "Non-COP", y = mergedData[group=="Non-COP"][sharing == 1]$N/2, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "Non-COP", y = mergedData[group=="COP"][sharing == 0]$N +  mergedData[group=="COP"][sharing == 1]$N/1.5, size = 8) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 2.5, size = 7.5, fontface = "bold"  ) +
  ylim(c(0,max(mergedData$N)+300)) +
  xlab("group") +
  ylab("Number of gene pairs") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

