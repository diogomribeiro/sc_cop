#07-10-2021 # Diogo Ribeiro @ UNIL
# Script to remove samples with too high number of COPs

library(data.table)
require(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("All arguments need to be provided", call.=FALSE)
}

inFile = args[1]
medataFile = args[2]
outFile = args[3]
iqrMulti = as.numeric(args[4])

data = fread(inFile, header = T, sep = "\t")

metadata = fread(medataFile, header = T, sep = "\t")

perDonorRun = data.table(table(data$dataset))
colnames(perDonorRun) = c("donor-run","dr_cops")

cells = data.table(table(metadata$`donor_id-run_id`))

perDonorRun = merge(perDonorRun,cells, by.x ="donor-run", by.y = "V1")


## Outlier removal
perDonorRun$ratio = perDonorRun$dr_cops/perDonorRun$N
iqr = IQR(perDonorRun$ratio)
quan3 = quantile(perDonorRun$ratio, 0.75)
ggplot( perDonorRun, aes(x = "x", y = ratio ))  +
  geom_boxplot( ) +
  geom_text(aes(x = "x", label = `donor-run`), position=position_jitter(height = 0.5, width = 0.1)) +
  geom_hline(yintercept = quantile(perDonorRun$ratio, 0.75) + iqr*iqrMulti, color = "red") +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


toRemove = perDonorRun[ratio > quantile(perDonorRun$ratio, 0.75) + iqr*iqrMulti]$`donor-run`

write.table(data[!dataset %in% toRemove],outFile, quote=F, sep="\t", row.names=F, col.names = T)
