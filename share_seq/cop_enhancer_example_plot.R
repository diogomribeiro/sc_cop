
library(data.table)
library(ggplot2)

wantedEnhStart = 38738400
wantedEnhEnd = 38739000
rnaData = fread("../source_data/share_seq/rep3.rna.counts.ensg.filtered.coding.binary.chr21.bed.gz", sep ="\t")

rnaData = rnaData[start > wantedEnhStart - 500000][end < wantedEnhEnd + 500000]
rnaData[,1:10]

atacData = fread("../processing/rep3_atac_peak_epimap_enhancer_merge_0.5_chr21.bed")
atacData$one = data.table(unlist(lapply(atacData$V7, function(x) unlist(strsplit(x,"[,]"))[1])))$V1
atacData$two = data.table(unlist(lapply(atacData$V7, function(x) unlist(strsplit(x,"[,]"))[2])))$V1
atacData$three = data.table(unlist(lapply(atacData$V7, function(x) unlist(strsplit(x,"[,]"))[3])))$V1
atacData$cell = paste(atacData$one,atacData$two,atacData$three,sep=",")
atacData$one = NULL
atacData$two = NULL
atacData$three = NULL

atacCells = unique(atacData$cell)

geneModels = rnaData[,1:4]
geneDT = data.table()
for (g in geneModels$gene){
  print(g)
  d = rnaData[gene == g]

  d = data.table(t(d[,7:ncol(d)]))
  d$cell = colnames(rnaData)[7:ncol(rnaData)]
  
  d$one = data.table(unlist(lapply(d$cell, function(x) unlist(strsplit(x,"[,]"))[1])))$V1
  d$two = data.table(unlist(lapply(d$cell, function(x) unlist(strsplit(x,"[,]"))[2])))$V1
  d$three = data.table(unlist(lapply(d$cell, function(x) unlist(strsplit(x,"[,]"))[3])))$V1
  d$cell = paste(d$one,d$two,d$three,sep=",")
  d = d[V1 == 1][cell %in% atacData$cell]
  
  data = data.table(cell = atacCells)
  data$value = 0
  data[cell %in% d$cell]$value = 1
  data$gene = g

  geneDT = rbind(geneDT, data)
}
geneMerge = merge(geneDT, geneModels, by = "gene", all.x = T)


#####
atacData = atacData[V2 > wantedEnhStart - 500000][V3 < wantedEnhEnd + 500000]
atacData$tag = paste(atacData$V1,atacData$V2,atacData$V3,sep="_")

enhDT = data.table()
for (a in unique(atacData$tag)){
  print(a)
  d = atacData[tag == a]
  data = data.table(cell = atacCells)
  data$value = 0
  data[cell %in% d$cell]$value = 1
  data$tag = a
  enhDT = rbind(enhDT, data)
}

atacModels = unique(atacData[,.(V1,V2,V3,tag)])
enhMerge = merge(enhDT, atacModels, by = "tag", all.x = T)

# Calculate frequency of each cell to rank them
mergedData = rbind(geneMerge[value == 1][,.(cell)], enhMerge[value == 1][,.(cell)])
cellFreq = data.table(table(mergedData$cell))
cellFreq = cellFreq[order(-N)]
cellFreq$rank = seq(1,nrow(cellFreq))
colnames(cellFreq)[1] = "cell"

enhMerge = merge(enhMerge,cellFreq, by = "cell",all.x=T)
geneMerge = merge(geneMerge,cellFreq, by = "cell",all.x=T)

ggplot() +
  geom_rect(data = geneModels, aes(xmin = start, xmax = end, ymin = -1800, ymax = -500), color = "white", fill = "white") + # just for ylim
  geom_rect(data = geneModels, aes(xmin = start, xmax = end, ymin = -500, ymax = -700), color = "black", fill = "#377eb8", size = 0.5) +
  geom_rect(data = atacModels, aes(xmin = V2-1000, xmax = V3+1000, ymin = -400, ymax = -800), color = "black", fill = "#4daf4a", size = 0.5) +
  # geom_rect(xmin = gene2Start, xmax = gene2End, ymin = 1000, ymax = 1010, color = "green", fill = "green") +
  geom_point(data = geneMerge[value == 1],aes(y = rank, x = (start+end)/2 ), color = "#377eb8", alpha = 0.01) +
  geom_point(data = enhMerge[value == 1],aes(y = rank, x = (V2+V3)/2 ), color = "#4daf4a", alpha = 0.01) +
  scale_y_reverse( ) +
  xlim(c( wantedEnhStart - 380000, wantedEnhEnd + 190000)) +
  xlab("Genome position") +
  ylab("Cells ranked by activity") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18),
      axis.text.x = element_text(vjust=0.6),
      legend.text=element_text(size=20), legend.title=element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1))

length(unique(geneMerge[!is.na(rank)]$cell))
