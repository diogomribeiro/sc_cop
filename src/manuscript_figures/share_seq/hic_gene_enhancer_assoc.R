#09-Aug-2021 Diogo Ribeiro @ UNIL
# Support for gene-enhancer interactions using Hi-C 

library(ggplot2)
library(data.table)

# data = fread( "/scratch/dribeir1/single_cell/multi_omics/hic/coex_F0.5_enhancer_gene/coex_peak_mergenh_F_0.5.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")
data = fread( "/scratch/dribeir1/single_cell/multi_omics/hic/coex_F0.5_enhancer_gene_coding_nofilt/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")
# data = fread( "/scratch/dribeir1/single_cell/multi_omics/hic/coex_F0.5_enhancer_gene_coding/coex_peak_mergenh_F_0.5_coding.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")
# data = fread( "/scratch/dribeir1/single_cell/multi_omics/hic/coex_F0.5_enhancer_gene_rep2/coex_peak_rep2_mergenh_F_0.5.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")

data = fread( "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/multi_omics/hic/25kb_normalised/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")
# data = fread( "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/multi_omics/hic/10kb_normalised/coex_peak_mergenh_F_0.5_coding_nofilt.bed_contacts", stringsAsFactors = FALSE, header = T, sep="\t")


data$fdr = p.adjust(data$corrPval, method = "BH")

peakAssociationFDRCutoff = 0.05
peakCorrCutoff = 0.05

data[is.na(normalised_contact)]$normalised_contact = 0
# data = data[!is.na(normalised_contact)]
data$normalised_contact = log2(data$normalised_contact+1)

# accounting for distance
data$tss = data$centralStart
data[centralStrand == "-"]$tss = data[centralStrand == "-"]$centralEnd
data$dist = abs(data$tss - data$cisStart)
cor.test(data$dist,data$normalised_contact)
data$res = residuals(lm(data$normalised_contact ~ data$dist ))

p = cor.test(data[real == 1]$corr,data[real == 1]$res, method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value))
ggplot(data[real == 1], aes(x = corr, y = res) ) +
  geom_bin_2d() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 8, fontface = "bold"  ) +
  xlab("Gene-enhancer association") +
  ylab("Hi-C contact (log distance-scaled)") +
  # scale_fill_gradient(low = "#e5f5f9", high = "#00441b") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )

p = cor.test(data[real == 0]$corr,data[real == 0]$res, method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value))
ggplot(data[real == 0], aes(x = corr, y = res) ) +
  geom_bin_2d() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 8, fontface = "bold"  ) +
  xlab("Gene-enhancer association") +
  ylab("Hi-C contact (log distance-scaled)") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )

# ## Distance vs contacts
# p = cor.test(data[real == 1]$dist,data[real == 1]$normalised_contact, method = "spearman")
# text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value))
# ggplot(data[real == 1], aes(x = dist, y = normalised_contact) ) +
#   geom_bin_2d() +
#   geom_smooth(method = "lm") +
#   annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
#   xlim(c(0,1000000)) +
#   xlab("Distance") +
#   ylab("Hi-C contact (log-scaled)") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1),
#   )


data$significant = "no interaction"
data[fdr < peakAssociationFDRCutoff][corrSign == "+"][corr > peakCorrCutoff]$significant = "interaction"

data$real = as.factor(data$real)
data[real == 1]$real = "real"
data[real == 0]$real = "control"
data$group = paste(data$significant,data$real)

table(data$group)

wilcox.test(data[group == "interaction real"]$res, data[group == "interaction control"]$res)
meansDF = data.table(aggregate(data$res, list(data$group), mean))
colnames(meansDF) = c("group","mean")
ggplot(data[group != "no interaction control"], aes(x = group, y = res, fill = group) ) +
  geom_violin() + 
  geom_boxplot(width = 0.2) +
  geom_text(data = meansDF[group != "no interaction control"], aes(x = group, y = mean, label = round(mean,1)), size = 7, fontface = "bold", color = "#525252", nudge_x = 0.35) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1),
  )

#######
# Matching each real-control case
#######
realData = data[group == "interaction real"][order(match_id)][,.(centralPhenotype,cisPheno,match_id, normalised_contact, res)]
controlData = data[group == "interaction control"][order(match_id)][,.(centralPhenotype,cisPheno, match_id, normalised_contact, res)]
sum(realData$match_id == controlData$match_id) == nrow(realData)
colnames(realData) = c("gene_real","enhancer_real","match_id_real","contact_real","res_real")
colnames(controlData) = c("gene_control","enhancer_control","match_id_control","contact_control","res_control")
mergedData = cbind(realData,controlData)
mergedData[res_real > res_control]
mergedData$contact_diff = mergedData$res_real - mergedData$res_control

text = paste("% cases real > control: ",round(nrow(mergedData[res_real > res_control]) * 100 / nrow(mergedData),1), "%", sep = "")

ggplot(mergedData, aes(x = contact_diff) ) +
  geom_histogram() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  annotate(geom = "text", label = text, y = Inf, x = Inf, hjust = 1, vjust = 2, size = 7) +
  ylab("# gene-enhancer combinations") +
  xlab("Contact difference") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=28), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


