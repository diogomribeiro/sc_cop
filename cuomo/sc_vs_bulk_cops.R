#26-02-2021 # Diogo Ribeiro @ UNIL
# Script to plot distance distribution

library(data.table)
require(ggplot2)

bulk = unique(fread("../source_data/cuomo/bulk_CODer_distance_controlled_null.bed_positive", stringsAsFactors=F, header=T, sep="\t"))
sc = unique(fread("../source_data/cuomo/CODer_final_dataset_cops_merged_removedoutliers.bed", stringsAsFactors=F, header=T, sep="\t"))
scInter = unique(fread("../source_data/cuomo/COPs_5_or_more.txt", stringsAsFactors=F, header=F, sep="\t"))

bulk = unique(bulk[significant == 1])

scInter[V1 %in% unique(bulk$pairID)]
scInter[!V1 %in% unique(bulk$pairID)]

sc$donorRun = sc$dataset

o1 = length(unique(bulk[pairID %in% sc$pairID]$pairID))
o2 = length(unique(sc[pairID %in% bulk$pairID]$pairID))

dt = data.table(group = c("BulkCOPs","BulkCOPs","scCOPs","scCOPs"),
                fill = c("3: Bulk+Sc","1: Bulk only","3: Bulk+Sc","2: Sc only"),
                # fill = c("1:",1,1,0,1,0),
                Number_COPs = c(o1, length(unique(bulk$pairID))-o1, o2,length(unique(sc$pairID))-o2))

text1 = paste("Total",sum(dt[group == "BulkCOPs"]$Number_COPs))
text2 = paste("Total",sum(dt[group == "scCOPs"]$Number_COPs))

ggplot( dt, aes(x = group, y = Number_COPs, fill = as.factor(fill) ))  +
  geom_bar( stat = "identity", color = "black", size = 1, width = 0.7, alpha = 0.5) +
  geom_text(aes(label = Number_COPs, y = Number_COPs/1.3), vjust = 1.1, size = 8) +
  annotate("text", x = "BulkCOPs", y = sum(dt[group == "BulkCOPs"]$Number_COPs), label = text1, vjust = -0.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = "scCOPs", y = sum(dt[group == "scCOPs"]$Number_COPs), label = text2, vjust = -0.5, size = 7, fontface = "bold"  ) +
  xlab("Dataset") +
  ylim(c(0,max(dt$Number_COPs) + 500) ) +
  ylab("Number of COPs") +
  scale_fill_brewer(palette= "Set1") +
  theme_linedraw() + 
  theme(text = element_text(size=24), 
        legend.position = "none", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


bulk$dataset = "BulkCOPs"
sc$dataset = "scCOPs"
distDT = rbind(unique(bulk[,.(pairID, distance,dataset)]), unique(sc[,.(pairID,distance,dataset)]))

distDT$detail = "bulk"

table(distDT$dataset)

ggplot( distDT, aes(x = distance, fill = dataset ))  +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette= "Set2") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        legend.text=element_text(size=24), legend.title=element_blank(), legend.position = c(0.8,0.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

