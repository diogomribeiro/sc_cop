
library(gprofiler2)
library(data.table)
options("scipen"=100, "digits"=2)

g1Genes = fread("../data/functional_enrichment/cuomo2021_g1_genes.txt", header = F)
sGenes = fread("../data/functional_enrichment/cuomo2021_s_genes.txt", header = F)
g2mGenes = fread("../data/functional_enrichment/cuomo2021_g2m_genes.txt", header = F)
backgroundGenes = fread("../data/functional_enrichment/cuomo2021_cycle_background_genes.txt", header = F)

gostres <- gost(query = list(g1Genes$V1,sGenes$V1,g2mGenes$V1),organism = "hsapiens", domain_scope = "custom", custom_bg = backgroundGenes$V1, sources = c("GO:BP","KEGG","REAC"), multi_query = T)

res = data.table(gostres$result)

res = res[,.(term_id, p_values, significant, source, term_name)]

res$p_val1 = data.table(unlist(lapply(res$p_values, function(x) unlist(x)[1])))$V1
res$p_val2 = data.table(unlist(lapply(res$p_values, function(x) unlist(x)[2])))$V1
res$p_val3 = data.table(unlist(lapply(res$p_values, function(x) unlist(x)[3])))$V1

res$sign1 = data.table(unlist(lapply(res$significant, function(x) unlist(x)[1])))$V1
res$sign2 = data.table(unlist(lapply(res$significant, function(x) unlist(x)[2])))$V1
res$sign3 = data.table(unlist(lapply(res$significant, function(x) unlist(x)[3])))$V1

res[sign1 == TRUE][sign2 == TRUE][sign3 == TRUE]

res[term_name == "Ribosome"]
res[term_name == "Cellular responses to stress"]

g1Res = data.table(res[sign1 == TRUE])
g1Res = g1Res[,.(term_id,source,term_name,p_val1)]
g1Res$p_val1 = as.numeric(g1Res$p_val1)
g1Res$log10_pval = -log10(g1Res$p_val1)
g1Res$p_val1 = format.pval(g1Res$p_val1,2)
g1Res = g1Res[order(-log10_pval)]

sRes = data.table(res[sign2 == TRUE])
sRes = sRes[,.(term_id,source,term_name,p_val2)]
sRes$p_val2 = as.numeric(sRes$p_val2)
sRes$log10_pval = -log10(sRes$p_val2)
sRes$p_val2 = format.pval(sRes$p_val2,2)
sRes = sRes[order(-log10_pval)]

g2mRes = data.table(res[sign3 == TRUE])
g2mRes = g2mRes[,.(term_id,source,term_name,p_val3)]
g2mRes$p_val3 = as.numeric(g2mRes$p_val3)
g2mRes$log10_pval = -log10(g2mRes$p_val3)
g2mRes$p_val3 = format.pval(g2mRes$p_val3,2)
g2mRes = g2mRes[order(-log10_pval)]


# write.table(g1Res,"gprofiler/g1_results.tsv", quote=F, row.names=F, sep = "\t")
# write.table(sRes,"gprofiler/s_results.tsv", quote=F, row.names=F, sep = "\t")
# write.table(g2mRes,"gprofiler/g2m_results.tsv", quote=F, row.names=F, sep = "\t")
