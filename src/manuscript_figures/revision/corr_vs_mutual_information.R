
library(data.table)
library(ggplot2)

miData = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/cop_indentification/share_seq_ma2020/mutual_information/1000R/final_dataset/CODer_raw_results.bed")
miData = miData[!is.na(distance)]
miData = unique(miData[,.(pairID,corr)])

corrData = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/cop_indentification/share_seq_ma2020/1000R_binary_before_norm_1MB/final_dataset/CODer_raw_results.bed")
corrData = corrData[!is.na(distance)]
corrData = unique(corrData[,.(pairID,corr)])

mergedData = merge(miData,corrData, by = "pairID")

hist(corrData$corr, breaks = 50)
hist(miData$corr, breaks = 50)

mergedData$corr.x.norm = normalize(mergedData$corr.x)
mergedData$corr.y.norm = normalize(mergedData$corr.y)

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

ggplot(mergedData,aes(x = corr.x, y = corr.y)) +
  geom_point(alpha = 0.1, shape = 21) +
  xlab("Mutual Information") +
  ylab("Spearman correlation") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=16), panel.grid.minor.y = element_line(size = 0.03, color = "black"), 
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),  panel.grid.major.y = element_line(size = 0.03, color = "black"),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1
  )


plot(mergedData$corr.x,mergedData$corr.y, xlab = "Mutual Information", ylab = "Spearman Correlation")
cor.test(mergedData$corr.x,mergedData$corr.y)
cor.test(mergedData$adjustedPval.x,mergedData$adjustedPval.y)
plot(mergedData$adjustedPval.x,mergedData$adjustedPval.y, xlab = "Mutual Information", ylab = "Spearman Correlation")

# Get the mean value relative to 0.2 correlation

summary(mergedData[corr.y > 0.199 & corr.y < 0.201]$corr.x)


discovery = mergedData[corr.y > 0.2]
replication = mergedData[corr.x > min(mergedData[corr.y > 0.2]$corr.x)]

discovery[pairID %in% replication$pairID]
replication[pairID %in% discovery$pairID]

# e.g. using correlation cutoff of 0.2 we have 14 significant cases, 13 of these are found with MI cutoff 0.01625.
# MI finds 2 other associations (at correlation 0.19).

mergedData[adjustedPval.y <= 0.000999][adjustedPval.x > 0.000999]

#########################

miData = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/cop_indentification/share_seq_ma2020/mutual_information/1000R/final_dataset/CODer_distance_controlled_null.bed_positive")
miData = miData[significant == 1]
miData = unique(miData[,.(pairID,corr,adjustedPval)])

corrData = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/cop_indentification/share_seq_ma2020/1000R_binary_before_norm_1MB/final_dataset/CODer_distance_controlled_null.bed_positive")
corrData = corrData[significant == 1]
corrData = unique(corrData[,.(pairID,corr,adjustedPval)])

miData[pairID %in% corrData$pairID]
corrData[pairID %in% miData$pairID]

m = merge(miData,corrData,all.x = T, all.y = T, by = "pairID")

missing = m[is.na(corr.y)]$pairID

summary(corrData[pairID %in% missing])


# 2420 out of 2589 (93.5%) correlation COPs are found with MI, with the same amount of COP discoveries between the two methods
2420/2589

