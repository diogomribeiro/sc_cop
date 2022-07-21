
library(gprofiler2)
library(data.table)
options("scipen"=100, "digits"=2)

g1Genes = fread("/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/G1_genes.txt", header = F)
sGenes = fread("/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/S_genes.txt", header = F)
g2mGenes = fread("/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/G2M_genes.txt", header = F)
backgroundGenes = fread("/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/background_genes.txt", header = F)

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


write.table(g1Res,"/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/g1_results.tsv", quote=F, row.names=F, sep = "\t")
write.table(sRes,"/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/s_results.tsv", quote=F, row.names=F, sep = "\t")
write.table(g2mRes,"/scratch/dribeir1/single_cell/cop_indentification/cuomo2021/sc_rna_seq/per_donor_per_experiment/cell_cycle/functional_enrichment/gprofiler/g2m_results.tsv", quote=F, row.names=F, sep = "\t")
