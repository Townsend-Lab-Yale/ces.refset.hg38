library(data.table)
library(rtracklayer)
library(cancereffectsizeR)

# Require cancereffectsizeR >= 2.3.1 for latest version of refset due to small changes in create_refset/build_RefCDS
stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('2.3.1'))

# External data source:
# Gencode GRCh38 "basic" GTF, release 38

# Loading Gencode GTF, release 38 (latest release as of 09/15/21)
gen = as.data.table(import("~/reference/gencode/gencode.v38.basic.annotation.gtf.gz"))

# Subset to protein-coding since we're building CDS annotations
# In the future, would be nice to include additional annotations somewhere
gen = gen[gene_type == 'protein_coding']
gen[, support := transcript_support_level]
gen[is.na(support) | support == "NA", support := '6'] # yeah, both "NA" and NAs
gen[, support := as.numeric(support)]


# For each gene, prioritize CCDS (consensus CDS) with highest transcript support, and then
# for genes with no such transcripts, accept the most-supported non-CCDS transcripts. When
# there are multiple transcripts tied for most supported, keep them all.
just_ccds = gen[type == 'CDS' & ! is.na(ccdsid), .SD[support == min(support)], by = "gene_id"]
missing_genes = setdiff(gen$gene_id, just_ccds$gene_id)
from_missing = gen[missing_genes, on = 'gene_id'][type == 'CDS', .SD[support == min(support)], by = "gene_id"]
best_by_gene = rbind(just_ccds, from_missing)

sum(width(reduce(makeGRangesFromDataFrame(best_by_gene)))) # 33.1 Mb (similar to old RefCDS)

## A few transcripts will be tossed by build_RefCDS, but most will be good
# uniqueN(best_by_gene$gene_id) # 19,955
# uniqueN(best_by_gene$transcript_id) # 30,654


# chr17:7675994 (TP53) is not at an "essential splice site" as defined by Martincorena et al.,
# but mutations there have experimentally confirmed effects on splicing.
# We will grab all TP53 transcripts that include this site and manually specify the splice status.
# (It turns out that for all such transcripts, 7675994 is at a CDS start position.)
## sort(best_by_gene[gene_name == "TP53"][, abs(7675994 - as.numeric(end))])[1:10]
## sort(best_by_gene[gene_name == "TP53"][, abs(7675994 - as.numeric(start))])[1:15]
custom_splice_list = best_by_gene[gene_name == "TP53" & start == 7675994, .(protein_id, start)]
custom_splice_list = setNames(as.list(custom_splice_list$start), custom_splice_list$protein_id)

refcds = build_RefCDS(gtf = best_by_gene, genome = 'hg38', use_all_transcripts = TRUE, 
                      additional_essential_splice_pos = custom_splice_list)

# For a default exome, we'll take all exon regions of protein-coding genes (with no transcript support filtering)
# This gives 80 Mb due to UTRs, as opposed to 34.7 Mb if just CDS regions were used.
exome = reduce(makeGRangesFromDataFrame(gen[type == 'exon']))

dir.create('tmp_ref') # will move data into refset package inst/refset directory

create_refset(output_dir = "tmp_ref/", refcds_output = refcds, species_name = 'human',
              genome_build_name = 'hg38', BSgenome_name = 'hg38', supported_chr = c(1:22, 'X', 'Y'), 
              default_exome = exome, exome_interval_padding = 0)

# Finish by copying into inst/refset