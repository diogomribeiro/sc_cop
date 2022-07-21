
library(data.table)
library(ggplot2)

inFile = "../source_data/cuomo/CODer_final_dataset_cops_merged_removedoutliers.bed.gz"

data = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")
data = unique(data[,.(pairID,dataset)])

allCOPs = unique(data$pairID)
  
data$donor = data.table(unlist(lapply(data$dataset, function(x) unlist(strsplit(x,"[-]"))[1])))$V1
data$experiment = data.table(unlist(lapply(data$dataset, function(x) unlist(strsplit(x,"[-]"))[2])))$V1

## Filter for at least 1 other experiment/donor
d = data.table(table(unique(data[,.(donor,experiment)])$donor))
wantedDonors = d[N > 1]$V1

finalDT = data.table()
for (dn in unique(wantedDonors)){
  dt = data[donor == dn]
  for (exp in unique(dt$experiment)){
    expCOPs = unique(dt[experiment == exp]$pairID)
    otherexpCOPs = unique(dt[experiment != exp]$pairID)
    real = sum(expCOPs %in% otherexpCOPs) * 100 / length(expCOPs)

    randDT = data.table()  
    for (i in seq(100)){
      randomCOPs = sample(allCOPs, length(unique(dt[experiment != exp]$pairID)))
      rand = sum(expCOPs %in% randomCOPs) * 100 / length(expCOPs)
      randDT = rbind(randDT, data.table(rand))
    }
    
    finalDT = rbind(finalDT, data.table(donor = dn, experiment = exp, observed = real, expected = mean(randDT$rand) ) )
    
  }
}


finalDT

meltedDT = melt(finalDT)

p = wilcox.test(finalDT$observed, finalDT$expected)
text = paste("Observed mean:", round(mean(finalDT$observed),1), "\nExpected mean: ", round(mean(finalDT$expected),1), "\nWilcoxon p-value", format.pval(p$p.value,1))
ggplot( meltedDT, aes(x = variable, y = value, color = variable)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(height = 0, width = 0.2) +
  annotate("text", x = Inf, y = Inf, label = text, vjust = 1.1, hjust = 1, size = 7, fontface = "bold"  ) +
  scale_fill_brewer( palette = "Set2") +
  xlab("Background group") + 
  ylab("% COPs replicated") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1, 
        legend.title=element_blank(), legend.position = "None", legend.text = element_text(size = 22))


