library(data.table)
library(rtracklayer)
library(cancereffectsizeR)

# Require cancereffectsizeR >= 2.4.0 for latest version of refset due to fixes in create_refset/build_RefCDS
stopifnot(packageVersion('cancereffectsizeR') >= as.package_version('2.4.0'))

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

# remove the handful of chrM records
best_by_gene = best_by_gene[seqnames != 'chrM']

sum(width(reduce(makeGRangesFromDataFrame(best_by_gene)))) # 33.1 Mb (similar to old RefCDS)

## A few transcripts will be tossed by build_RefCDS, but most will be good
# uniqueN(best_by_gene$gene_id) # 19,942
# uniqueN(best_by_gene$transcript_id) # 30,641


# A handful of well-studied genes like TP53 have lots of transcripts. For simplicity,
# we'll find transcripts that do not have any exons not already incorporated in longer transcripts,
# and we'll filter these out.
transcripts_by_gene = split(best_by_gene, by = 'gene_id')
best_by_gene[, exon_number := as.numeric(exon_number)]
exons_per_transcript = best_by_gene[, .(num_exons = max(exon_number)), by = 'transcript_id']

best_by_gene[exons_per_transcript, num_exons := num_exons, on = 'transcript_id']
best_by_gene = best_by_gene[order(gene_id, -num_exons)]

find_novel_transcripts = function(gene_dt) {
  used_exons = character()
  unneeded_transcripts = character()
  gene_dt = setDT(gene_dt[, .(transcript_id, exon_id)], key = 'transcript_id')
  setkey(gene_dt, 'transcript_id')
  unique_transcripts = unique(gene_dt$transcript_id)
  for (transcript in unique_transcripts) {
    current_exons = gene_dt[transcript, exon_id]
    new_exons = setdiff(current_exons, used_exons)
    if(length(new_exons) == 0) {
      unneeded_transcripts = c(unneeded_transcripts, transcript)
    }
    used_exons = c(used_exons, new_exons)
  }
  return(setdiff(unique_transcripts, unneeded_transcripts))
}

# Takes a minute to run
useful_transcripts_by_gene = best_by_gene[, .(transcript_id = find_novel_transcripts(.SD)), by = 'gene_id']

best_by_gene = best_by_gene[useful_transcripts_by_gene$transcript_id, on = 'transcript_id']
# Same number of genes (19,942), but now down to 29,475 transcripts.
# uniqueN(best_by_gene$gene_id)
# uniqueN(best_by_gene$transcript_id)


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
exome = reduce(makeGRangesFromDataFrame(gen[type == 'exon' & seqnames != 'chrM']))

dir.create('tmp_ref') # will move data into refset package inst/refset directory

create_refset(output_dir = "tmp_ref/", refcds_output = refcds, species_name = 'human',
              genome_build_name = 'hg38', BSgenome_name = 'hg38', supported_chr = c(1:22, 'X', 'Y'), 
              default_exome = exome, exome_interval_padding = 0)

# Finish by copying into inst/refset