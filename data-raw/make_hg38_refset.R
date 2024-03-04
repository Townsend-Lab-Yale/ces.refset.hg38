library(data.table)
library(cancereffectsizeR) # v2.9.0
library(GenomicRanges)

# Fun fact: rtracklayer's GTF parser doesn't handle tags correctly.
gtf_file = 'untracked/gencode.v45.basic.annotation.gtf.gz' # local Gencode GTF
gtf = as.data.table(rtracklayer::import(gtf_file))

# Read in file again to get tags          
gtf$tag = NULL
tags = fread(gtf_file, skip = 'chr') # first data line starts with a chromosome name

# 9th column (V9) has tag entries (along with other stuff)
stopifnot(identical(tags$V4, gtf$start)) # Verify records are in order 

gtf[, tag := tags$V9]

gtf[, support := transcript_support_level]
gtf[is.na(support) | support == "NA", support := '6'] # yeah, both "NA" and NAs
gtf[, support := as.numeric(support)]
gtf[, exon_number := as.numeric(exon_number)]
gtf[, seqnames := gsub('^chr', '', as.character(seqnames))]
gtf = gtf[seqnames %in% c(1:22, 'X', 'Y')]


gtf[, is_mane := tag %like% 'tag \\"MANE_Select']
gtf[, is_mane_plus := tag %like% 'tag \\"MANE_Plus_Clinical']
# gtf[, is_basic := tag %like% 'tag \\"basic\\"'] # all records are basic since we're now using the basic GTF
gtf[, is_seleno := tag %like% 'tag \\"seleno']
gtf[, is_non_ATG_start := tag %like% 'tag \\"non_ATG_start']
gtf[, is_NMD_exception := tag %like% 'tag \\"NMD_exception']

# Capture various splice annotations (all "non_canonical" tags relate to splicing)
gtf[, noncanonical_splice := tag %like% 'tag \"non_canonical']
gtf[, NAGNAG_splice_site := tag %like% 'tag \"NAGNAG_splice_site']

cds = gtf[transcript_type == 'protein_coding' & type == 'CDS']

uniqueN(cds$transcript_id) # 63,727
sum(width(reduce(makeGRangesFromDataFrame(cds)))) # 34,760,548

cds[, .(uniqueN(gene_id), uniqueN(transcript_id))] # 19,681 genes with 63,727 transcripts

# chr17:7675994 (TP53) is not at an "essential splice site" as defined by Martincorena et al.,
# but mutations there have experimentally confirmed effects on splicing.
# We will grab all TP53 transcripts that include this site and manually specify the splice status.
# (It turns out that for all such transcripts, 7675994 is at a CDS start position.)
## sort(cds[gene_name == "TP53"][, abs(7675994 - as.numeric(end))])[1:10]
## sort(cds[gene_name == "TP53"][, abs(7675994 - as.numeric(start))])[1:15]
custom_splice_list = cds[gene_name == "TP53" & start == 7675994, .(protein_id, start)]
custom_splice_list = setNames(as.list(custom_splice_list$start), custom_splice_list$protein_id)


refcds_anno = build_RefCDS(gtf = cds, genome = 'hg38', use_all_transcripts = TRUE, 
                           additional_essential_splice_pos = custom_splice_list, cores = 4)

# There will be a warning about some of the protein IDs in the custom_splice_list not being used; that's okay.
refcds_dndscv = build_RefCDS(gtf = cds, genome = 'hg38', use_all_transcripts = FALSE, 
                             additional_essential_splice_pos = custom_splice_list, cores = 4)

# For a default exome, we'll take all exon regions of protein-coding genes (with no transcript support filtering)
# This gives 81 Mb due to UTRs, as opposed to ~35 Mb if just CDS regions were used.
exome = reduce(makeGRangesFromDataFrame(gtf[type == 'exon' & transcript_type == 'protein_coding']))

# We will store some information about coding and noncoding transcripts
transcript_info = gtf[(transcript_type == 'protein_coding' & type %in% c('UTR', 'transcript', 'CDS')) |
                        (transcript_type != 'protein_coding' & type == 'transcript'), 
                      .(chr = seqnames, start, end, strand, type, gene_id, gene_type, gene_name, transcript_id,
                        transcript_type, exon_number = as.integer(exon_number), protein_id, is_mane, is_mane_plus, non_ATG_start = is_non_ATG_start,
                        NMD_exception = is_NMD_exception, seleno = is_seleno, noncanonical_splice, NAGNAG_splice_site)]

# Verify these columns are transcript-level, then combine into one field
transcript_level_cols = c('non_ATG_start', 'NMD_exception', 'noncanonical_splice', 'NAGNAG_splice_site', 'seleno')
for(col in transcript_level_cols) {
  stopifnot(uniqueN(transcript_info[, .SD, .SDcols = c('transcript_id', col)]) == uniqueN(transcript_info$transcript_id))
}

transcript_melted = unique(melt(transcript_info, id.vars = 'transcript_id', measure.vars = transcript_level_cols)[value == TRUE])
transcript_melted[, variable := paste0('<', variable, '>')]
transcript_melted[, transcript_tags := paste(variable, collapse = ''), by = 'transcript_id']
to_merge = unique(transcript_melted[, .(transcript_id, transcript_tags)])

transcript_info[to_merge, transcript_tags := transcript_tags, on = 'transcript_id']
transcript_info = transcript_info[, .SD, .SDcols = setdiff(names(transcript_info), transcript_level_cols)]

# Note that makeGRangesFromDataFrame(transcript_info[transcript_type == 'protein_coding' & type %in% c('UTR', 'CDS')])
# reproduces the default exome.

dir.create('tmp_ref') # will move data into refset package inst/refset directory

create_refset(output_dir = "tmp_ref/", refcds_anno = refcds_anno, refcds_dndscv = refcds_dndscv,
              species_name = 'human', genome_build_name = 'hg38', BSgenome_name = 'hg38', 
              supported_chr = c(1:22, 'X', 'Y'), default_exome = exome, exome_interval_padding = 0,
              transcripts = transcript_info)









