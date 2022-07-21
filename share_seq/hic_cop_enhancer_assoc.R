#09-Aug-2021 Diogo Ribeiro @ UNIL
# Support for gene-enhancer interactions using Hi-C 

library(ggplot2)
library(data.table)

data = fread( "../source_data/share_seq/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts_5kb.gz", stringsAsFactors = FALSE, header = T, sep="\t")
#data = fread( "../source_data/share_seq/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts_5kb.gz", stringsAsFactors = FALSE, header = T, sep="\t")
#data = fread( "../source_data/share_seq/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts_5kb.gz", stringsAsFactors = FALSE, header = T, sep="\t")

copFile = "../source_data/share_seq/CODer_distance_controlled_null.bed_positive"

peakAssociationFDRCutoff = 0.05
peakCorrCutoff = 0.05
wantedQuantile = 0.75

data$fdr = p.adjust(data$corrPval, method = "BH")

data[is.na(normalised_contact)]$normalised_contact = 0
data$normalised_contact = log2(data$normalised_contact+1)

# accounting for distance
data$tss = data$centralStart
data[centralStrand == "-"]$tss = data[centralStrand == "-"]$centralEnd
data$dist = abs(data$tss - data$cisStart)
data$normalised_contact = residuals(lm(data$normalised_contact ~ data$dist ))

data$significant = "no interaction"
data[fdr < peakAssociationFDRCutoff][corrSign == "+"][corr > peakCorrCutoff]$significant = "interaction"
data$real = as.factor(data$real)
data[real == 1]$real = "real"
data[real == 0]$real = "control"
data$group = paste(data$significant,data$real)
table(data$group)
quantileContact = quantile(summary(data[group != "no interaction control"][group != "interaction control"]$normalised_contact), probs = wantedQuantile)

### Read gene-gene co-expression
geneGeneData = fread( copFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

geneGeneData$centralPhenotype = data.table(unlist(lapply(geneGeneData$centralPhenotype, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
geneGeneData$centralPhenotype = data.table(unlist(lapply(geneGeneData$centralPhenotype, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
geneGeneData$cisPheno = data.table(unlist(lapply(geneGeneData$cisPheno, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
geneGeneData$cisPheno = data.table(unlist(lapply(geneGeneData$cisPheno, function(x) unlist(strsplit(x,"[.]"))[1])))$V1

table(geneGeneData$significant)
geneGeneData = unique(geneGeneData[,.(pairID,centralPhenotype,cisPheno,distance,significant,nullId)])

## Interaction data
enhancerGeneData = unique(data[,.(centralPhenotype,cisPheno,normalised_contact,group)])

## Merge peak-gene and gene-gene data
colnames(enhancerGeneData) = c("centralPhenotype","enhancer","hic1","group1")
m1 = merge(geneGeneData,enhancerGeneData, by = "centralPhenotype", allow.cartesian = T, all.x = T)
colnames(enhancerGeneData) = c("cisPheno","enhancer","hic2","group2")
mergedData = merge(m1,enhancerGeneData, by = c("cisPheno","enhancer"), allow.cartesian = T, all.x = T)
mergedData$significant = as.factor(mergedData$significant)
# mergedData = mergedData[group1 == group2]

## Q1: For the same COP-enhancer pair, how many times the real data has higher contacts than control (distance-matched)?
d1 = mergedData[group1 == "interaction real"][group2 == "interaction real"][,.(pairID,enhancer,hic1,hic2,group1,significant)]
d2 = mergedData[group1 == "interaction control"][group2 == "interaction control"][,.(pairID,enhancer,hic1,hic2,group1,significant)]
colnames(d2) = c("pairID_control","enhancer_control","hic1_control","hic2_control","group1_control","significant_control")
d3 = cbind(d1,d2)

d3$mean1 = (d3$hic1 + d3$hic2) / 2
d3$mean2 = (d3$hic1_control + d3$hic2_control) / 2
d3[mean1 > mean2]

text = paste("% cases real > control: ",round(nrow(d3[significant == 1][mean1 > mean2]) * 100 / nrow(d3[significant == 1]),1), "%", sep = "")
d3$contact_diff = d3$mean1 - d3$mean2
ggplot(d3[significant == 1], aes(x = contact_diff) ) +
  geom_histogram() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  annotate(geom = "text", label = text, y = Inf, x = Inf, hjust = 1, vjust = 2, size = 6) +
  xlab("Real-control Hi-C contact difference") +
  ylab("# COP-enhancer combinations") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
  )

## Q2: # gene pairs with Hi-C support (for at least 1 enhancer)
d5 = mergedData[group1 == "interaction real" & group2 == "interaction real"]

d = unique(mergedData[,.(pairID, significant)])
table(d$significant)

d5$pass = 0
d5[hic1 > quantileContact][hic2 > quantileContact]$pass = 1

v1 = length(unique(d5[significant == 1 & pass == 1]$pairID))
v2 = length(unique(d5[significant == 0 & pass == 1]$pairID))
d6 = data.table( group = c("COP", "COP","Non-COP", "Non-COP"), 
                 hic_support = c(1,0,1,0),
                 N = c( v1,  length(unique(mergedData[significant == 1]$pairID)) - v1, v2, length(unique(mergedData[significant == 0]$pairID)) - v2 ),
                 fill = c(3,2,1,0)                 )


perc1 = d6[group == "COP"][hic_support == 1]$N * 100.0 / (d6[group == "COP"][hic_support == 1]$N + d6[group == "COP"][hic_support == 0]$N)
perc2 = d6[group == "Non-COP"][hic_support == 1]$N * 100.0 / (d6[group == "Non-COP"][hic_support == 1]$N + d6[group == "Non-COP"][hic_support == 0]$N)
ggplot(d6, aes(x = group, y = N, fill = as.factor(fill) ) ) +
  geom_bar(stat = "identity", color = "black", size = 1) + 
  scale_fill_manual( values = c("#d9d9d9","#969696","#fb9a99","#e31a1c"), label = c("Not supported","Supported", "Not supported", "Supported")) +
  annotate(geom = "text", label = paste0(round(perc1,1),"%"), x = "COP", y = d6[group=="COP"][hic_support == 1]$N/2, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc1,1),"%"), x = "COP", y = d6[group=="COP"][hic_support == 0]$N + d6[group=="COP"][hic_support == 1]$N/1.8, size = 8) +
  annotate(geom = "text", label = paste0(round(perc2,1),"%"), x = "Non-COP", y = d6[group=="Non-COP"][hic_support == 1]$N/2, size = 8) +
  annotate(geom = "text", label = paste0(round(100-perc2,1),"%"), x = "Non-COP", y = d6[group=="COP"][hic_support == 0]$N +  d6[group=="COP"][hic_support == 1]$N / 5, size = 8) +
  ylab("Number of gene pairs") + 
  xlab("group") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
