
library(gprofiler2)
library(data.table)
options("scipen"=100, "digits"=2)

scGenes = fread("../data/functional_enrichment/cuomo2021_sc_genes.txt", header = F)
bulkGenes = fread("../data/functional_enrichment/cuomo2021_bulk_genes.txt", header = F)
backgroundGenes = fread("../data/functional_enrichment/cuomo2021_background_genes.txt", header = F)

gostres <- gost(query = list(scGenes$V1,bulkGenes$V1),organism = "hsapiens", domain_scope = "custom", custom_bg = backgroundGenes$V1, sources = c("GO:BP","KEGG","REAC"), multi_query = T)

res = data.table(gostres$result)

res = res[,.(term_id, p_values, significant, source, term_name)]

res$p_val1 = data.table(unlist(lapply(res$p_values, function(x) unlist(x)[1])))$V1
res$p_val2 = data.table(unlist(lapply(res$p_values, function(x) unlist(x)[2])))$V1

res$sign1 = data.table(unlist(lapply(res$significant, function(x) unlist(x)[1])))$V1
res$sign2 = data.table(unlist(lapply(res$significant, function(x) unlist(x)[2])))$V1

res[sign2 == TRUE][sign1 == TRUE]
res[sign1 == TRUE][sign2 == FALSE]

scGeneRes = res[sign1 == TRUE][source == "GO:BP"][,.(term_id, p_val1)]
bulkGeneRes = res[sign2 == TRUE][source == "GO:BP"][,.(term_id, p_val2)]

scRes = data.table(res[sign1 == TRUE])
scRes = scRes[,.(term_id,source,term_name,p_val1)]
scRes$p_val1 = as.numeric(scRes$p_val1)
scRes$log10_pval = -log10(scRes$p_val1)
scRes$p_val1 = format.pval(scRes$p_val1,2)
scRes = scRes[order(-log10_pval)]

res[term_id == "GO:0006412"]
res[term_id == "GO:0043933"]
res[term_name == "Ribosome"]
scRes[term_name == "Ribosome"]
res[term_name == "Cell Cycle"]
scRes[term_name == "Cell Cycle"]
res[term_name == "Cellular responses to stress"]
scRes[term_name == "Cellular responses to stress"]


bulkRes = data.table(res[sign2 == TRUE])
bulkRes = bulkRes[,.(term_id,source,term_name,p_val2)]
bulkRes$p_val2 = as.numeric(bulkRes$p_val2)
bulkRes$log10_pval = -log10(bulkRes$p_val2)
bulkRes$p_val2 = format.pval(bulkRes$p_val2,2)
bulkRes = bulkRes[order(-log10_pval)]

# write.table(scRes,"gprofiler/sc_results.tsv", quote=F, row.names=F, sep = "\t")
# write.table(bulkRes,"gprofiler/bulk_results.tsv", quote=F, row.names=F, sep = "\t")

# write.table(scGeneRes,"gprofiler/sc_gene_BP_results.tsv", quote=F, sep = "\t",row.names=F, col.names=F)
# write.table(bulkGeneRes,"gprofiler/bulk_gene_BP_results.tsv", quote=F, sep = "\t",row.names=F, col.names=F)
