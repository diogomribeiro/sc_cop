
library(ggplot2)
library(data.table)

## EPIMAP

data = fread( "../source_data/share_seq/coex_epimap_intersect.bed.gz", stringsAsFactors = FALSE, header = F, sep="\t")
colnames(data) = c("chr_coex","start_coex","end_coex","gene_coex","corr","corr_sign","corr_pval","chr_epi","start_epi","end_epi","gene_epi","epi_score")
data[corr_sign == "-"]$corr = -data[corr_sign == "-"]$corr

data$corr_norm = qnorm((rank(data$corr))/(length(data$corr)+1))
data$epi_score_norm = qnorm((rank(data$epi_score))/(length(data$epi_score)+1))

p = cor.test(data$corr_norm,data$epi_score_norm,method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value,2))
ggplot(data, aes(x = corr_norm, y = epi_score_norm) ) +
  geom_bin_2d( bins = 25) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  xlab("SHARE-seq coex correlation") +
  ylab("EpiMap score") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


## ABC model

data = fread( "../source_data/share_seq/coex_abc_intersect.bed", stringsAsFactors = FALSE, header = F, sep="\t")
colnames(data) = c("chr_coex","start_coex","end_coex","gene_coex","corr","corr_sign","corr_pval","chr_abc","start_abc","end_abc","gene_abc","abc_score")
data[corr_sign == "-"]$corr = -data[corr_sign == "-"]$corr

data$corr_norm = qnorm((rank(data$corr))/(length(data$corr)+1))
data$abc_score_norm = qnorm((rank(data$abc_score))/(length(data$abc_score)+1))

p = cor.test(data$corr_norm,data$abc_score_norm,method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value,2))
ggplot(data, aes(x = corr_norm, y = abc_score_norm) ) +
  geom_bin_2d() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  xlab("SHARE-seq coex correlation") +
  ylab("ABC model score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


## Epi vs ABC

data = fread( "../source_data/share_seq/epi_abc_intersect.bed", stringsAsFactors = FALSE, header = F, sep="\t")
colnames(data) = c("chr_epi","start_epi","end_epi","gene_epi","epi_score","chr_abc","start_abc","end_abc","gene_abc","abc_score")

data$epi_score_norm = qnorm((rank(data$epi_score))/(length(data$epi_score)+1))
data$abc_score_norm = qnorm((rank(data$abc_score))/(length(data$abc_score)+1))

p = cor.test(data$epi_score_norm,data$abc_score_norm,method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", format.pval(p$p.value,2))
ggplot(data, aes(x = epi_score_norm, y = abc_score_norm) ) +
  geom_bin_2d() +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 7, fontface = "bold"  ) +
  xlab("EpiMap score") +
  ylab("ABC model score") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )

