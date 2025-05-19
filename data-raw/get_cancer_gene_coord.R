# Run from refset development directory
library(devtools)
library(data.table)

load_all()

# Downloaded from http://www.network-cancer-genes.org/
# See 10.1186/s13059-022-02607-z
cg = fread('inst/extdata/NCG_cancerdrivers_annotation_supporting_evidence.tsv.gz')
cg = unique(cg[symbol %in% ces.refset.hg38$gene_names, .(symbol, NCG_oncogene, NCG_tsg)])
all_cancer_genes = cg$symbol


# Protein-coding genes again
# Leaving out chrY. (Note some genes have X and Y listings in the transcript table.
transcripts = ces.refset.hg38$transcripts[gene_name %in% ces.refset.hg38$gene_names & chr != 'Y']

# Didn't lose any cancer genes by dropping chrY.
stopifnot(all(all_cancer_genes %in% transcripts$gene_name))

pure_oncogenes = cg[NCG_oncogene == 1, symbol]
pure_tsg = cg[NCG_tsg == 1, symbol]
stopifnot(cg[NCG_oncogene == 1 & NCG_tsg == 1, .N] == 0)
transcripts[, cancer_anno := fcase(gene_name %in% pure_oncogenes, 'oncogene',
                                   gene_name %in% pure_tsg, 'TSG',
                                   gene_name %in% all_cancer_genes, 'other',
                                   default = 'noncancer')]


gene_id_counts = transcripts[, .(cancer_anno = cancer_anno[1], num_id = uniqueN(gene_id)), by = 'gene_name']
multi_id_genes = gene_id_counts[num_id > 1, gene_name]


# Gene names with multiple gene IDs are all non-cancer. For simplicity, we'll just take gene ID from the MANE entry.
# We lose two genes that have no MANE entry; fine.
stopifnot(sum(multi_id_genes %in% all_cancer_genes) == 0)
transcripts = transcripts[is_mane == T | ! gene_name %in% multi_id_genes]


general_coord = transcripts[, .(cancer_anno = cancer_anno[1], chr = chr[1], start = min(start), end = max(end)), 
                            by = 'gene_name']

mane_info = transcripts[is_mane == T, .(mane_start = min(start), mane_end = max(end)), 
                        by = 'gene_name']

gene_coord = merge.data.table(general_coord, mane_info, by = 'gene_name', all.x = T)
gene_coord[, chr_factor := factor(chr, levels = c(1:22, 'X'))]
gene_coord = gene_coord[order(chr_factor, start, end)]
gene_coord$chr_factor = NULL
setnames(gene_coord, 'gene_name', 'gene')
saveRDS(gene_coord, 'inst/refset/gene_coord_cancer_anno.rds')

