
library(data.table)
require(ggplot2)

bulk = unique(fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/bulk/15pcs_genefilt/final_dataset/CODer_distance_controlled_null.bed_positive", stringsAsFactors=F, header=T, sep="\t"))
sc = unique(fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/scran/CODer_final_dataset_cops_merged_removedoutliers.bed", stringsAsFactors=F, header=T, sep="\t"))

# Filter for individuals used in bulk data
listDonors = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/geuvadis/yoruba/individuals_used",header = F)
sc = sc[dataset %in% listDonors$V1]

write.table(sc$pairID,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/scran/cops_37_individuals.txt",quote=F,row.names=F,col.names=F)
freq = data.table(table(sc$pairID))
repCops = freq[N >= 5]$V1
write.table(repCops,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/scran/cops_37_individuals_5_or_more.txt",quote=F,row.names=F,col.names=F)


# Filter for genes used in sc data
listGenes = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/sarkar2019/matrices_scran/list_genes_used",header = F)
bulk = unique(bulk[significant == 1])
bulk$gene1 = data.table(unlist(lapply(bulk$centralPhenotype, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
bulk$gene2 = data.table(unlist(lapply(bulk$cisPheno, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
bulk = bulk[gene1 %in% listGenes$V1][gene2 %in% listGenes$V1]

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
  geom_bar( stat = "identity", color = "black", size = 1, width = 0.6, alpha = 0.7) +
  geom_text(aes(label = Number_COPs, y = Number_COPs/1.3), vjust = 1.1, size = 8) +
  annotate("text", x = "BulkCOPs", y = sum(dt[group == "BulkCOPs"]$Number_COPs), label = text1, vjust = -0.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = "scCOPs", y = sum(dt[group == "scCOPs"]$Number_COPs), label = text2, vjust = -0.5, size = 7, fontface = "bold"  ) +
  xlab("Dataset") +
  ylim(c(0,max(dt$Number_COPs) + 500) ) +
  ylab("Number of COPs") +
  scale_fill_brewer(palette= "Set2") +
  theme_linedraw() + 
  theme(text = element_text(size=24), 
        legend.position = "none", legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

## Confusion table

testedPairs=fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/both/overlap/intersection.txt", header =F, sep="\t")
bulkInter=testedPairs[V1 %in% bulk$pairID]
scInter=testedPairs[V1 %in% sc$pairID]

tp=nrow(bulkInter[V1 %in% scInter$V1])
fp=nrow(bulkInter[!V1 %in% scInter$V1])-tp
fn=nrow(scInter[!V1 %in% bulkInter$V1])-tp
tn=nrow(testedPairs)-(fp+fn+tp)
m = matrix(c(tp,fp,fn,tn),nrow=2)
m
f = fisher.test(m)
f$p.value
f$estimate
