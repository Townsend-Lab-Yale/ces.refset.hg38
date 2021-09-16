library(data.table)
library(rtracklayer)


# External data sources:
# Gencode GRCh38 "basic" GTF, release 38

# Loading Gencode GTF, release 38 (latest release as of 09/15/21)
gen = as.data.table(import("~/reference/gencode/gencode.v38.basic.annotation.gtf.gz"))

# Subset to protein-coding since we're building CDS annotations
# In the future, would be nice to include additional annotations somewhere
gen = gen[gene_type == 'protein_coding']
gen[, support := transcript_support_level]
gen[is.na(support) | support == "NA", support := '6'] # yeah, both "NA" and NAs
gen[, support := as.numeric(support)]


# For each gene, take the CDS regions for only its most-supported transcripts
best_by_gene = gen[type == 'CDS', .SD[support == min(support)], by = "gene_id"]
sum(width(reduce(makeGRangesFromDataFrame(best_by_gene)))) # 32.9 Mb (basically same as old RefCDS)

## A few transcripts will be tossed by build_RefCDS, but most will be good
# uniqueN(best_by_gene$gene_id) # 19,955
# uniqueN(best_by_gene$transcript_id) # 33,865


# chr17:7675994 (TP53) is not at an "essential splice site" as defined by Martincorena et al.,
# but mutations there have experimentally confirmed effects on splicing.
# We will grab all TP53 transcripts that include this site and manually specify the splice status.
# (It turns out that for all such transcripts, 7675994 is at a CDS start position.)
## sort(best_by_gene[gene_name == "TP53"][, abs(7675994 - as.numeric(end))])[1:10]
## sort(best_by_gene[gene_name == "TP53"][, abs(7675994 - as.numeric(start))])[1:15]
custom_splice_list = best_by_gene[gene_name == "TP53" & start == 7675994, .(protein_id, start)]
custom_splice_list = setNames(as.list(custom_splice_list$start), custom_splice_list$protein_id)

refcds = build_RefCDS(gtf = best_by_gene, genome = 'hg38', use_all_transcripts = T, 
                      additional_essential_splice_pos = custom_splice_list)

# For a default exome, we'll take all coding regions (with no transcript support filtering)
exome = reduce(makeGRangesFromDataFrame(gen[type == 'CDS']))

# 100bp padding would bring in most commonly callable UTR regions, and it would increase
# the total size of the regions from 34.7 Mb to 73.1 Mb. Instead, we'll leave padding at
# 0, since users can keep out-of-coverage data anyway. preload_maf() can help a user
# decide whether uncovered calls are mostly valid calls from wide exome capture kits, or
# identify likely sequencing/calling error.
dir.create('tmp_ref') # will move data into refset package inst/refset directory

export.bed(exome, 'hg38_gencode38_cds_regions.bed')
create_refset(output_dir = "tmp_ref/", refcds_output = refcds, species_name = 'human',
              genome_build_name = 'hg38', BSgenome_name = 'hg38', supported_chr = c(1:22, 'X', 'Y'), 
              default_exome_bed = 'hg38_gencode38_cds_regions.bed', exome_interval_padding = 0)
