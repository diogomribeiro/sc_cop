#15-Oct-2021 Diogo Ribeiro @ UNIL
options(scipen=1)

library(data.table)
library(ggplot2)

copFile = "/scratch/dribeir1/single_cell/cop_indentification/share_seq_ma2020/1000R_binary_before_norm_1MB/final_dataset/CODer_distance_controlled_null.bed_positive"
# copFile = "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_identification/GTEx/all_tissues/Cells_EBV-transformed_lymphocytes/final_dataset/CODer_distance_controlled_null.bed_positive"

copData =  fread( copFile, stringsAsFactors = FALSE, header = T, sep="\t")
copData = unique(copData[,.(centralPhenotype,cisPheno,significant, distance, nullId)])
table(copData$significant)
copData$centralPhenotype = data.table(unlist(lapply(copData$centralPhenotype, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
copData$cisPheno = data.table(unlist(lapply(copData$cisPheno, function(x) unlist(strsplit(x,"[.]"))[1])))$V1

genes = unique(c(copData$centralPhenotype,copData$cisPheno))

groFile = "/scratch/dribeir1/single_cell/raw_input/gro_seq_gm12878/share_seq_tss_gro_seq_support_allreads.bed"

groData = fread( groFile, stringsAsFactors = FALSE, header = F, sep="\t")
groData$readSign = "+"
groData[V9 < 0]$readSign = "-"
groData = groData[V5 == readSign]
groData$V9 = abs(groData$V9)
groData = unique(groData[,.(V4,V9)])
groData = groData[V4 %in% genes]
groData$V9 = rank(groData$V9)

mergedData = merge(copData,groData, by.x = "centralPhenotype", by.y = "V4", all.x = T)
mergedData = merge(mergedData,groData, by.x = "cisPheno", by.y = "V4", all.x = T)

table(mergedData$significant)

## Correlation of reads
# mergedData[is.na(V9.x)]$V9.x = 0
# mergedData[is.na(V9.y)]$V9.y = 0
mergedData$dataset = "-1"
mergedData[significant == 1]$dataset = "COP"
mergedData[significant == 0]$dataset = "Non-COP"

mergedData = mergedData[!is.na(V9.x)][!is.na(V9.y)]
table(mergedData$significant)

correlation = cor.test(mergedData[significant == 1]$V9.x,mergedData[significant == 1]$V9.y, method = "spearman")
correlationText = paste("COP Spearman R:",round(correlation$estimate,3), "P-value:",format.pval(correlation$p.value,2), sep = " ")
correlation2 = cor.test(mergedData[significant == 0]$V9.x,mergedData[significant == 0]$V9.y, method = "spearman")
correlationText2 = paste("Non-COP Spearman R:",round(correlation2$estimate,3), "P-value:",format.pval(correlation2$p.value,2), sep = " ")
cor.test(mergedData[significant == 1]$V9.x,mergedData[significant == 1]$V9.y, method = "pearson")
cor.test(mergedData[significant == 0]$V9.x,mergedData[significant == 0]$V9.y, method = "pearson")

ggplot(mergedData[dataset == "COP"], aes(x = V9.x, y = V9.y) ) +
  geom_bin2d(binwidth = 200) +
  geom_point(alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1.05, vjust = 1.5, size = 6.5, fontface = "bold"  ) +
  ylab("Gene1 reads (rank)") +
  xlab("Gene2 reads (rank)") +
  ggtitle("COP GRO-seq support") +
  scale_fill_gradient(low = "#e5f5e0", high = "#00441b") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1,
        legend.title=element_blank(), legend.text = element_text(size = 20))

ggplot(mergedData[dataset == "Non-COP"], aes(x = V9.x, y = V9.y) ) +
  geom_bin2d(binwidth = 200) +
  geom_point(alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = correlationText2, hjust = 1.05, vjust = 1.5, size = 6.5, fontface = "bold"  ) +
  ylab("Gene1 reads (rank)") +
  xlab("Gene2 reads (rank)") +
  ggtitle("Non-COP GRO-seq support") +
  scale_fill_gradient(low = "#f0f0f0", high = "#252525") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1,
        legend.title=element_blank(), legend.text = element_text(size = 20))

