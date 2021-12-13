#13-May-2021 Diogo Ribeiro @ UNIL
# Exploring overlap of proteomics data with COPs and non-COPs

library(data.table)
library(ggplot2)

options(scipen=1)
args = commandArgs(trailingOnly=TRUE)

##########
# Input parameters
##########

proteinFile = args[1] #"mirauta_2020_elife/elife-57390-fig1-data3-v3.xls"
proteinData = fread( proteinFile, stringsAsFactors = FALSE, header = T, sep="\t")

wantedGenes = proteinData$ensembl_gene_id

nCOPCutoff = 10 #minimum number of COPs of donor-run for analysis 

listDonorRun = fread( "../data/cuomo2021_list_individual_experiment.tsv", stringsAsFactors = FALSE, header = F, sep="\t")
baseFolder = args[2] # all_runs/

##########
# Process protein data (sum intensities for several isoforms of same gene)
##########
ncol(proteinData)

# Remove HPSI0314i-bubh_3 (control)
# Remove Strontium
dropCols = grep("Strontium", colnames(proteinData))
proteinData[, (dropCols) := NULL]
dropCols = grep("bubh", colnames(proteinData))
proteinData[, (dropCols) := NULL]
ncol(proteinData)

geneCol = grep("ensembl_gene_id", colnames(proteinData))
wantedCols = grep("HPSI", colnames(proteinData))
d = data.table(ensembl_gene_id = proteinData$ensembl_gene_id)
geneLabels = unique(d[order(ensembl_gene_id)])

finalDT = data.table()
for (c in wantedCols){
  co = c(c,geneCol)
  dt = proteinData[,..co]
  print(colnames(dt)[1])
  res = data.table(aggregate(dt[,1], list(intensity = dt$ensembl_gene_id), sum ))
  finalDT = cbind(finalDT, data.table(res[,2]))
  
}
finalDT = cbind(finalDT,geneLabels)

geneNames = finalDT$ensembl_gene_id

##########
# List of donors in common
##########
listDonorRun$donor = data.table(unlist(lapply(listDonorRun$V1, function(x) unlist(strsplit(x,"[-]"))[1])))$V1

proteomicsDonors = data.table(colnames(proteinData))

proteomicsDonors$donor = data.table(unlist(lapply(proteomicsDonors$V1, function(x) unlist(strsplit(x,"[-]"))[2])))$V1
proteomicsDonors$donor = data.table(unlist(lapply(proteomicsDonors$donor, function(x) unlist(strsplit(x,"[@]"))[1])))$V1

donorMerge = merge(listDonorRun,proteomicsDonors, by = "donor", all.x = T)

donorMerge = donorMerge[!is.na(V1.y)]
length(unique(donorMerge[is.na(V1.y)]$donor))

##########
# Loop each donor
##########

countCOPs = c()
results = data.table()
for (dn in unique(donorMerge$donor)){ # for each donor
    # Get list of runs per donor               
  folder = listDonorRun[donor == dn]$V1
  
  for (f in folder){ # for each donor-run
    print(f)
    
    # Get COP/non-COP data
    copFile = paste(baseFolder, f, "/final_dataset/CODer_distance_controlled_null.bed_positive",sep="")
    copData = fread( copFile, stringsAsFactors = FALSE, header = T, sep="\t")

    if(nrow(copData)/2 > nCOPCutoff){
      # print(nrow(copData))
      
      countCOPs = c(countCOPs, unique(copData[significant == 1]$pairID))
      
      # # # Control: shuffling gene pairs
      # copData$centralPhenotype = copData[sample(nrow(copData))]$centralPhenotype
      # copData$cisPheno = copData[sample(nrow(copData))]$cisPheno
      
      # Get protein data
      protCol = unique(donorMerge[donor == dn]$V1.y)
      proteinData = finalDT[,get(protCol),"ensembl_gene_id"]
      
      colnames(proteinData) = c("centralPhenotype","intensity_1")
      mergedData = merge(copData, proteinData, by = "centralPhenotype", all.x =T)
      colnames(proteinData) = c("cisPheno","intensity_2")
      mergedData = merge(mergedData, proteinData, by = "cisPheno", all.x =T)
      
      mergedData = mergedData[centralPhenotype %in% geneNames][cisPheno %in% geneNames]
      
      # mergedData[is.na(intensity_1)]$intensity_1 = 0
      # mergedData[is.na(intensity_2)]$intensity_2 = 0
      
      correlation = cor.test(mergedData[significant == 1]$intensity_1,mergedData[significant == 1]$intensity_2, method = "spearman")
      print(paste("COP Spearman R:",round(correlation$estimate,2), "P-value:",format.pval(correlation$p.value,2), sep = " "))
      correlation2 = cor.test(mergedData[significant == 0]$intensity_1,mergedData[significant == 0]$intensity_2, method = "spearman")
      print(paste("Non-COP Spearman R:",round(correlation2$estimate,2), "P-value:",format.pval(correlation2$p.value,2), sep = " "))
      
      results = rbind(results, data.table(donor_run = f, sample_size = nrow(copData)/2,
                                          cop_corr_r = correlation$estimate, cop_corr_pval = correlation$p.value, 
                                          noncop_corr_r = correlation2$estimate, noncop_corr_pval = correlation2$p.value ))
      
    }else{ print("No enough COPs for analysis") }
    
  }
}

results

length(unique(countCOPs))

meltedResults = melt(results, id.vars = c("donor_run", "sample_size"), measure.vars = c("cop_corr_r","noncop_corr_r"))

# summarize donor-run into donor
meltedResults$donor = data.table(unlist(lapply(meltedResults$donor_run, function(x) unlist(strsplit(x,"[-]"))[1])))$V1
length(unique(meltedResults$donor))

text1 = paste("COP mean correlation:",round(mean(meltedResults[variable == "cop_corr_r"]$value),2))
text2 = paste("Non-COP mean correlation",round(mean(meltedResults[variable == "noncop_corr_r"]$value),3))
t = wilcox.test(meltedResults[variable == "cop_corr_r"]$value, meltedResults[variable == "noncop_corr_r"]$value)
text = paste("Wilcoxon test P-value", format.pval(t$p.value,2) )
meltedResults[variable == "cop_corr_r"]$variable = "COP"
meltedResults[variable == "noncop_corr_r"]$variable = "Non-COP"
ggplot(meltedResults, aes(x = variable, y = value, fill = variable) ) +
  geom_abline(intercept = 0, slope = 0, color = "grey", linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.3, size = 1) +
  geom_jitter(width = 0.3, alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = text1, hjust = 1.05, vjust = 1.5, size = 6, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = text2, hjust = 1.05, vjust = 3.5, size = 6, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 5.5, size = 6, fontface = "bold"  ) +
  ylab("Correlation") +
  ylim(c(-1,1)) +
  xlab("Gene pair group") +
  # ggtitle("COP proteomics support across donors") +
  scale_fill_manual( values = c("#66C2A5","#969696")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)
