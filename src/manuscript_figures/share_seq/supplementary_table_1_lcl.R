library(data.table)

inFile = "zcat ../data/shareseq_sc_cops_control.bed.gz"

data = fread(inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
data = data[significant == 1]

data$chr = data$`#chr`

data$gene1_tss = data$centralStart
data[centralStrand == "-"]$gene1_tss = data[centralStrand == "-"]$centralEnd

data$gene2_tss = data$cisStart
data[cisStrand == "-"]$gene2_tss = data[cisStrand == "-"]$cisEnd

data$gene1_name = data.table(unlist(lapply(data$centralInfo, function(x) unlist(strsplit(x,"[=]"))[5])))$V1
data$gene2_name = data.table(unlist(lapply(data$cisInfo, function(x) unlist(strsplit(x,"[=]"))[5])))$V1

data$gene1 = data$centralPhenotype
data$gene2 = data$cisPheno

data$corr = round(data$corr,2)
finalDT = data[,.(chr,gene1,gene1_name,gene1_tss,gene2,gene2_name,gene2_tss,corr)]

