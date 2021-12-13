#27-Jan-2021 # Diogo Ribeiro @ UNIL
# Script to filter gene pairs by several conditions

library(data.table)

# number of digits after comma
options("scipen"=100, "digits"=2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("All arguments need to be provided", call.=FALSE)
}

########
# Load/process data
########

inFile = args[1] # CODer_raw_results.bed
corrCutoff = as.numeric(args[2])
excludedGenePairsFile = args[3]
outFile = args[4]

mergedData = fread( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
excludeData = fread( excludedGenePairsFile, stringsAsFactors = FALSE, header = F, sep="\t")

currentNrow = nrow(mergedData)
currentNrowNonredundant = length(unique(mergedData$pairID))

paste("Initial number of gene pairs:", currentNrow)
paste("Initial number of non-redundant gene pairs:", currentNrowNonredundant)

#######
# Remove extremely high correlation 
#######
nas = mergedData[is.na(corr)]
mergedData = mergedData[corr < corrCutoff]
mergedData = rbind(mergedData,nas)

paste("Remove correlation > ", corrCutoff,": ", currentNrow - nrow(mergedData), " pairs removed", sep="") 
currentNrow = nrow(mergedData)

#######
# Remove gene pairs
#######
mergedData = mergedData[!pairID %in% excludeData$V1]

paste("Remove pairs from list: ", currentNrow - nrow(mergedData), "pairs removed")
paste("Remove non-redundant pairs from list: ", currentNrowNonredundant - length(unique(mergedData$pairID)), "pairs removed")
currentNrow = nrow(mergedData)
currentNrowNonredundant = length(unique(mergedData$pairID))
#######

#######                                                                                                            
# Remove non-COP one-way gene pairs
#######
# Because of these filters, some gene pairs were removed in only one direction
# e.g. gene1-gene2 but not gene2-gene1. Which causes problems in subsequent code
toRemove = data.table(table(mergedData$pairID))
mergedData = mergedData[!pairID %in% toRemove[N != 2]$V1]

paste("Remove pairs for consistency: ", currentNrow - nrow(mergedData), "pairs removed")
paste("Remove pairs for consistency: ", currentNrowNonredundant - length(unique(mergedData$pairID)), "pairs removed")

########
# Write file
########
paste("Final dataset: ", nrow(mergedData), "gene pairs") 
paste("Final dataset: ", length(unique(mergedData$pairID)), "non-redudant gene pairs") 
write.table(mergedData, outFile, row.names = FALSE, quote = FALSE, sep ="\t")

