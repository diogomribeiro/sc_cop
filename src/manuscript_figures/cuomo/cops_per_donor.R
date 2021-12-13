#04-Oct-2021 Diogo Ribeiro @ UNIL
# COPs across donor-run

options(scipen=1)

library(data.table)
library(ggplot2)

inFile = "../data/cuomo2021_sc_cops_final_dataset.bed.gz"
medataFile = "../data/cuomo2021_CellSampleMetadata_reformat.txt"

data = fread(inFile, header = T, sep = "\t")
metadata = fread(medataFile, header = T, sep = "\t")

summary(data$corr)

data$originalDataset = data$dataset
# summarize donor-run into donor
data$dataset = data.table(unlist(lapply(data$dataset, function(x) unlist(strsplit(x,"[-]"))[1])))$V1

dt = unique(data[,.(pairID,dataset)])

length(unique(dt$pairID))


## COPs per donor
d = data.table(table(dt$dataset))

## COPs per donor vs sample size
# sampleSize = data.table(table(metadata$`donor_id-run_id`))
sampleSize = data.table(table(metadata$donor))
summary(sampleSize$N)
m = merge(d,sampleSize, by = "V1")
colnames(m) = c("donor","cops","cells")
summary(m$cops)
sd(m$cops)

correlation = cor.test(m$cops, m$cells, method = "spearman")
correlationText = paste("Spearman R:",round(correlation$estimate,2), "P-value:",format.pval(correlation$p.value,2, eps = -Inf), sep = " ")
meanText = paste("Mean # COPs: ", round(mean(m$cops),1), "Mean # cells:",round(mean(m$cells),1) )
ggplot( m, aes(x = cells, y = cops ))  +
  geom_point( ) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = meanText, hjust = 1.05, vjust = 3.5, size = 7, fontface = "bold"  ) +
  xlab("Number of cells") +
  ylab("Number of COPs") +
  ylim(c(0,700)) +
  theme_linedraw() + 
  theme(text = element_text(size=24), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)


## COP replication
dt = unique(dt[,.(pairID,dataset)])
dd = data.table(table(dt$pairID))
dd$perc = dd$N * 100 / length(unique(dt$dataset))
ggplot( dd, aes(x = perc ))  +
  geom_histogram( binwidth = 2, color = "black", size = 1) +
  xlab("% individuals") +
  ylab("# COPs (log10-scale)") +
  scale_y_log10() +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        axis.text.x = element_text(vjust=0.6),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


d = data.table(c("1 individual","2 to 5","5 or more"), c(nrow(dd[N==1]),nrow(dd[N>1][N<5]),nrow(dd[N>=5]) ) )
# write.table(unique(dd[N>=5]$V1),"COPs_5_or_more.txt", quote = F,row.names=F,col.names=F)

text = paste("Total COPs:", nrow(dd))
ggplot( d, aes(x = V1, y = V2, fill = V1))  +
  geom_bar( stat = "identity", color = "black", size = 1, alpha = 0.8) +
  geom_text(aes(label = V2), size = 8, vjust = 2) +
  ggtitle("COP replication") +
  xlab("Number of individuals") +
  ylab("Number of COPs") +
  scale_fill_brewer(palette = "Set1") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        axis.text.x = element_text(angle = 20, vjust=0.6),  legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

length(unique(data$pairID))


data$run = data.table(unlist(lapply(data$originalDataset, function(x) unlist(strsplit(x,"[-]"))[2])))$V1
data = unique(data[,.(pairID,dataset,run)])
perDonorRun = data.table(table(data$dataset, data$run))
colnames(perDonorRun) = c("donor","run","dr_cops")

perDonorRun = unique(perDonorRun[dr_cops>0])
sum(perDonorRun$dr_cops)

perDonorRun[dr_cops > 500]

perDonorRun = merge(perDonorRun,m, by ="donor")

ggplot( perDonorRun, aes(x = reorder(donor,dr_cops, sum), y = dr_cops, fill = run))  +
  geom_bar( stat = "identity", color = "black", size = 1, alpha = 0.8) +
  xlab("Individuals") +
  ylab("# COPs") +
  coord_flip() +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        axis.text.y = element_text(size = 12), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
