#04-Oct-2021 Diogo Ribeiro @ UNIL
# COPs across donor-run

options(scipen=1)

library(data.table)
library(ggplot2)

inFile = "/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/G1/CODer_final_dataset_cops_merged_removedoutliers.bed"
# inFile = "/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/S/CODer_final_dataset_cops_merged_removedoutliers.bed"
# inFile = "/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/G2M/CODer_final_dataset_cops_merged_removedoutliers.bed"
data = fread(inFile, header = T, sep = "\t")

### IMPORTANT: update this!
cyclePhase = "G1"

medataFile = "/scratch/dribeir1/single_cell/raw_input/cuomo2021/sc_rna_seq/regress_out_cycle_ribo/cell_metadata_cycle.txt"
metadata = fread(medataFile, header = T, sep = "\t")

data$originalDataset = data$dataset
data$dataset = data.table(unlist(lapply(data$dataset, function(x) unlist(strsplit(x,"[-]"))[1])))$V1
dt = unique(data[,.(pairID,dataset)])

table(unique(data[,.(pairID,corrSign)])$corrSign)

# Total COPs
length(unique(dt$pairID))
# COPs per donor
d = data.table(table(dt$dataset))

## COPs per donor vs sample size
# sampleSize = data.table(table(metadata$`donor_id-run_id`))
sampleSize = data.table(table(metadata[Phase == cyclePhase]$donor))
summary(sampleSize$N)
m = merge(d,sampleSize, by = "V1")
colnames(m) = c("donor","cops","cells")
summary(m$cops)

correlation = cor.test(m$cops, m$cells, method = "spearman")
correlationText = paste("Spearman R:",round(correlation$estimate,2), "P-value:",format.pval(correlation$p.value,2), sep = " ")
meanText = paste("Mean # COPs: ", round(mean(m$cops),1), "Mean # cells:",round(mean(m$cells),1) )
ggplot( m, aes(x = cells, y = cops ))  +
  geom_point( ) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = meanText, hjust = 1.05, vjust = 3.5, size = 7, fontface = "bold"  ) +
  ggtitle(paste(cyclePhase,"phase")) +
  xlab("Number of cells") +
  ylab("Number of COPs") +
  # coord_flip() +
  # xlim(c(0,400)) +
  # ylim(c(0,700)) +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio=1)
