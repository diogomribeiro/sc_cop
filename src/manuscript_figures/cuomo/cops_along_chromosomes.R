#21-May-2020 Diogo Ribeiro @ UNIL
# Plotting COPs along chromosomes

# pdf("cops_along_chromosomes.pdf",16,8)

library(data.table)
require(ggplot2)
library(grid)
library(gridExtra)

resFile = paste("../data/cuomo2021_raw_gene_models.bed.gz",sep="")
copFile = paste("../data/cuomo2021_sc_cops_union_dataset.bed.gz",sep="")

resData = fread( resFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
copData = fread( copFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

resData$chr = resData$`#chr`
copData$chr = copData$`#chr`

# get only central phenotype info
centralData = unique(resData[,.(chr,centralStart,centralEnd,centralPhenotype)])
centralData = centralData[order(chr)]

# adding gene order (sequencial gene IDs)
geneNumbers = c()
numberGenesPerChr = data.table(table(centralData$chr))$N
for (i in seq(1:length(numberGenesPerChr))){
  geneOrder = seq(1:numberGenesPerChr[i])
  geneNumbers = c(geneNumbers, geneOrder)
}
centralData$geneOrder = geneNumbers
# splitting centralData into genes in COPs or not
centralData$inCops = 0
listGenesInCops = unique(c(copData$centralPhenotype, copData$cisPheno))
centralData[centralData$centralPhenotype %in% listGenesInCops]$inCops = 1

### Proportion of genes in COPs for each chromosome
dt = data.table(table(centralData$chr,centralData$inCops))
dt2 = data.table(table(centralData$chr))
propCOP = merge(dt, dt2, by = "V1")
propCOP$perc = propCOP$N.x * 100.0 / propCOP$N.y
colnames(propCOP) = c("chr","inCops","value","total","percentage")

summary(propCOP[inCops == 1]$percentage)

### Plot of COPs along chromosomes (gene index)
ggplot( ) +
  geom_rect( data = centralData, aes(xmin = geneOrder-0.5, xmax = geneOrder+0.5, ymin = chr-0.22, ymax = chr+0.22, fill = as.factor(inCops) ) ) +
  geom_hline(aes(yintercept = centralData$chr-0.4), color = "black", size = 0.5 ) +
  geom_hline(aes(yintercept = centralData$chr+0.4), color = "black", size = 0.5 ) +
  geom_text(data = propCOP[inCops == 1], aes(x = 2200, y = chr, label = paste(round(percentage,1),"%", sep = "")), size = 4.5, fontface = "bold", hjust = 0.7) +
  scale_fill_manual(values=c("black", "#e41a1c"), labels=c("Non-COP", "COP")) + 
  xlab("Gene index") +
  ylab("Chromosomes") +
  scale_y_discrete(limits = seq(1,max(centralData$chr),1)) +
  # scale_fill_manual(values = c("black","#00CED1")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")


# dev.off()
