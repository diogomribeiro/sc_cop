# Plotting results from fisher exact test and its odds ratio on groups

library(data.table)
require(ggplot2)
library(grid)
library(gridExtra)

inputFile = "../data/functional_enrichment/cuomo2021_functional_enrichment_cell_cycle.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

colnames(dataset) = c("ExternalList","Size1","ItemGroup","Size2","Overlap","OddsRatio","Pvalue")

background = 265392

dataset$theSum = background - (dataset$Size1 + dataset$Size2) - dataset$Overlap
dataset$theFirst = dataset$Size1 - dataset$Overlap
dataset$theSecond = dataset$Size2 - dataset$Overlap

dataset$odds = apply(dataset, 1, function(x) fisher.test(matrix(c(as.numeric(x[5]),as.numeric(x[9]),as.numeric(x[10]),as.numeric(x[8]) ),nrow = 2), conf.level = 0.95)$estimate )
dataset$confmin = apply(dataset, 1, function(x) fisher.test(matrix(c(as.numeric(x[5]),as.numeric(x[9]),as.numeric(x[10]),as.numeric(x[8]) ),nrow = 2), conf.level = 0.95)$conf.int[1] )
dataset$confmax = apply(dataset, 1, function(x) fisher.test(matrix(c(as.numeric(x[5]),as.numeric(x[9]),as.numeric(x[10]),as.numeric(x[8]) ),nrow = 2), conf.level = 0.95)$conf.int[2] )
dataset$pval = apply(dataset, 1, function(x) fisher.test(matrix(c(as.numeric(x[5]),as.numeric(x[9]),as.numeric(x[10]),as.numeric(x[8]) ),nrow = 2), conf.level = 0.95)$p.value )

# number of digits after comma
options("scipen"=100, "digits"=2)

# Python float point limit is 2.225e-308
dataset[Pvalue == 0]$Pvalue = 2.225e-308
dataset$log10pval = -log10(dataset$Pvalue)

dataset = dataset[ExternalList != "paralog_genes"]
dataset[ItemGroup == "1: all phases"]$ItemGroup = "all phases"
dataset[ItemGroup == "2: G1 phase"]$ItemGroup = "G1 phase"
dataset[ItemGroup == "3: S phase"]$ItemGroup = "S phase"
dataset[ItemGroup == "4: G2M phase"]$ItemGroup = "G2M phase"
dataset[ExternalList == "same_complex"]$ExternalList = "Same complex"
dataset[ExternalList == "same_pathway"]$ExternalList = "Same pathway"
dataset[ExternalList == "conserved_cops"]$ExternalList = "Conserved COPs"

# Odds ratio plot
g1 = ggplot( ) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(data = dataset[ItemGroup == "all phases"], aes(x = confmin, xend = confmax, y = ExternalList, yend = ExternalList, color = ItemGroup), position = position_nudge(y = 0.3), size = 3 ) +
  geom_point(data = dataset[ItemGroup == "all phases"], aes(x = odds, y = ExternalList, fill = ItemGroup), position = position_nudge(y = 0.3), size = 6, shape = 21) + 
  geom_segment(data = dataset[ItemGroup == "G1 phase"], aes(x = confmin, xend = confmax, y = ExternalList, yend = ExternalList, color = ItemGroup), position = position_nudge(y = 0.1), size = 3 ) +
  geom_point(data = dataset[ItemGroup == "G1 phase"], aes(x = odds, y = ExternalList, fill = ItemGroup), position = position_nudge(y = 0.1), size = 6, shape = 24) + 
  geom_segment(data = dataset[ItemGroup == "S phase"], aes(x = confmin, xend = confmax, y = ExternalList, yend = ExternalList, color = ItemGroup), position = position_nudge(y = -0.1), size = 3 ) +
  geom_point(data = dataset[ItemGroup == "S phase"], aes(x = odds, y = ExternalList, fill = ItemGroup), position = position_nudge(y = -0.1), size = 6, shape = 24) +
  geom_segment(data = dataset[ItemGroup == "G2M phase"], aes(x = confmin, xend = confmax, y = ExternalList, yend = ExternalList, color = ItemGroup), position = position_nudge(y = -0.3), size = 3 ) +
  geom_point(data = dataset[ItemGroup == "G2M phase"], aes(x = odds, y = ExternalList, fill = ItemGroup), position = position_nudge(y = -0.3), size = 6, shape = 24) +
  xlab("Odds ratio (log-scale)") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  ylab("") +
  scale_x_log10() +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.position = c(0.72,0.92),
        legend.text=element_text(size=16), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


dataset$percOverlap = dataset$Overlap * 100 / dataset$Size2
g2 = ggplot() +
  geom_bar(data = dataset[ItemGroup == "all phases"], stat = "identity", aes(x = percOverlap, y = ExternalList, fill = ItemGroup), position = position_nudge(y = 0.25), color = "black",  width = 0.15) +
  geom_bar(data = dataset[ItemGroup == "G1 phase"], stat = "identity", aes(x = percOverlap, y = ExternalList, fill = ItemGroup), position = position_nudge(y = 0.10), color = "black",  width = 0.15) +
  geom_bar(data = dataset[ItemGroup == "S phase"], stat = "identity", aes(x = percOverlap, y = ExternalList, fill = ItemGroup), position = position_nudge(y = -0.05), color = "black",  width = 0.15) +
  geom_bar(data = dataset[ItemGroup == "G2M phase"], stat = "identity", aes(x = percOverlap, y = ExternalList, fill = ItemGroup), position = position_nudge(y = -0.20), color = "black",  width = 0.15) +
  scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#999999","#a65628","#f781bf")) +
  xlab("% overlap") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none", legend.key = element_rect(size = 1), legend.key.size = unit(1, 'lines'), legend.title = element_text(size = 12), legend.text = element_text(size = 10), 
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
        axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0.2,0.8,0.2,0), "cm"))

lay <- rbind(c(1,1,1,1,2))
grid.arrange(g1,g2, layout_matrix = lay)

summary(dataset$odds)

summary(dataset$log10pval)
