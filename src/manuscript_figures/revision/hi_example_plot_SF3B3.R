
library(data.table)
library(ggplot2)

args=commandArgs(trailingOnly=TRUE)

pack = "16_70473140_70473140" #args[1] #

c = data.table(unlist(lapply(pack, function(x) unlist(strsplit(x,"_"))[1])))$V1
wantedEnhStart = as.numeric(data.table(unlist(lapply(pack, function(x) unlist(strsplit(x,"_"))[2])))$V1) 
wantedEnhEnd = as.numeric(data.table(unlist(lapply(pack, function(x) unlist(strsplit(x,"_"))[3])))$V1) 

# c = args[1] #16
# wantedEnhStart = as.numeric(args[2]) #70286198
# wantedEnhEnd = as.numeric(args[3]) #70286198

windowSize = 150000

paste0(c,"_",wantedEnhStart,"_",wantedEnhEnd)

# png(paste0(c,"_",wantedEnhStart,"_",wantedEnhEnd,".png"),800,600)

inFile = paste0("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/cod_analysis/geuvadis/HiC/5kb_normalised/chr",c,"_5kb.RAWobserved_KR_normalised.gz")
data = fread(inFile)
data = data[V1 > wantedEnhStart - windowSize & V2 < wantedEnhStart + windowSize]

inFile = paste0("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/share_seq/rna_seq/processing/rep3/binary_nofilt/rep3.rna.counts.ensg.filtered.coding.binary.chr",c,".bed")
rnaData = fread(inFile, sep ="\t")
rnaData = rnaData[start > wantedEnhStart - windowSize & end < wantedEnhEnd + windowSize]
geneModels = rnaData[,1:6]
geneModels$tss = geneModels$start
geneModels[strand == "-"]$tss = geneModels[strand == "-"]$end
geneModels$geneName = data.table(unlist(lapply(geneModels$info, function(x) unlist(strsplit(x,"="))[5])))$V1

inFile = paste0("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/share_seq/atac_seq/processing/rep3_atac_peak_epimap_enhancer_merge_0.5_chr",c,".bed")
atacData = fread(inFile)
atacData = atacData[V2 > wantedEnhStart - windowSize & V3 < wantedEnhEnd + windowSize]
atacData$tag = paste(atacData$V1,atacData$V2,atacData$V3,sep="_")
atacModels = unique(atacData[,.(V1,V2,V3,tag)])
colnames(atacModels) = c("chr","start","end","tag")

data$start = (data$V1 + data$V2) / 2 
data$end = (data$V2 - data$V1) 
data$contact = log10(data$V3)

ggplot() + 
  geom_point(data=data,  aes(x = start/1000, y = end/1000, fill = contact, color = contact ), size = 3.4, shape = 16) +
  geom_rect(data = geneModels, aes(xmin = start/1000, xmax = end/1000, ymin = -30, ymax = -40), color = "black", fill = "#377eb8", size = 0.5) +
  geom_text(data = geneModels[geneName != "RP11-529K1.3"][strand == "+"], aes(x = (start+end)/2/1000, y = -10, label = geneName), color = "black",size=3) +
  geom_text(data = geneModels[strand == "-"], aes(x = (start+end)/2/1000, y = -60, label = geneName), color = "black",size=3) +
  geom_rect(data = atacModels, aes(xmin = start/1000-1, xmax = end/1000+1, ymin = -25, ymax = -45), color = "black", fill = "#4daf4a", size = 0.5) +
  geom_rect(data = geneModels, aes(xmin = tss/1000, xmax = tss/1000, ymin = -20, ymax = -40), color = "black", fill = "black", size = 0.5) +
  geom_rect(data = geneModels[strand == "+"], aes(xmin = tss/1000, xmax = tss/1000+5, ymin = -20, ymax = -20), color = "black", fill = "black", size = 0.5) +
  geom_rect(data = geneModels[strand == "-"], aes(xmin = tss/1000-5, xmax = tss/1000, ymin = -20, ymax = -20), color = "black", fill = "black", size = 0.5) +
  scale_color_gradient2(low = "white", midpoint = 1.8, high = "#e41a1c") +
  scale_fill_gradient2(low = "white", midpoint = 1.8, high = "#e41a1c") +
  ylim(c(-60000/1000,max(data$end)/1000)) +
  # xlim(c(wantedEnhStart-windowSize, wantedEnhEnd + windowSize)) +
  ylab("Distance from bin (kb)") +
  xlab("Position (kb)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# dev.off()



