
library(ggplot2)
library(data.table)

data = fread( "zcat /home/dribeiro/git/sc_cop/data/git/sc_cop/data/share_seq_gene_enhancer_associations.tsv.gz", stringsAsFactors = FALSE, header = T, sep="\t")
