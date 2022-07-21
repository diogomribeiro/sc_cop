#04-Oct-2021 Diogo Ribeiro @ UNIL
# COPs across donor-run

options(scipen=1)

library(data.table)
library(ggplot2)

inFile = "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/scran/CODer_final_dataset_cops_merged_removedoutliers.bed"

data = fread(inFile, header = T, sep = "\t")

# Filter for individuals used in bulk data
listDonors = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/raw_input/geuvadis/yoruba/individuals_used",header = F)
data = data[dataset %in% listDonors$V1]


summary(data$corr)

data$originalDataset = data$dataset
# summarize donor-run into donor
data$dataset = data.table(unlist(lapply(data$dataset, function(x) unlist(strsplit(x,"[-]"))[1])))$V1

dt = unique(data[,.(pairID,dataset)])

length(unique(dt$pairID))

medataFile = "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/raw_input/sarkar2019/sarkar2019_metadata.txt"
metadata = fread(medataFile, header = T, sep = "\t")

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
  # ggtitle("Number COPs per number cells") +
  xlab("Number of cells") +
  ylab("Number of COPs") +
  # coord_flip() +
  # xlim(c(0,400)) +
  ylim(c(0,350)) +
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
  geom_histogram( binwidth = 4, color = "black", size = 1) +
  # ggtitle("COP replication") +
  xlab("% individuals") +
  ylab("# COPs") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        axis.text.x = element_text(vjust=0.6),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

d = data.table(c("1 individual","2 to 5","5 or more"), c(nrow(dd[N==1]),nrow(dd[N>1][N<5]),nrow(dd[N>=5]) ) )
# write.table(unique(dd[N>=5]$V1),"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/manuscript/revision/sarkar2019/scran/COPs_5_or_more.txt", quote = F,row.names=F,col.names=F)

text = paste("Total COPs:", nrow(dd))
ggplot( d, aes(x = V1, y = V2, fill = V1))  +
  geom_bar( stat = "identity", color = "black", size = 1, alpha = 0.8) +
  geom_text(aes(label = V2), size = 8, vjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 5.5, size = 9, fontface = "bold"  ) +
  ggtitle("COP replication") +
  xlab("Number of individuals") +
  ylab("Number of COPs") +
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=30), 
        axis.text.x = element_text(angle = 20, vjust=0.6),  legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

length(unique(data$pairID))
