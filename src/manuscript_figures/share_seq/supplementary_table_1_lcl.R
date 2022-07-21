library(data.table)

inFile = "/scratch/dribeir1/single_cell/cop_indentification/share_seq_ma2020/1000R_binary_before_norm_1MB/final_dataset/CODer_distance_controlled_null.bed_positive"

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

write.table(finalDT, "/home/dribeiro/Downloads/supplementary_table_1_lcl.tsv",sep = "\t",quote = F, col.names=T, row.names=F)

