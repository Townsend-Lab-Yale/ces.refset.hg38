# Generate site counts for each COSMIC DBS type. Each dinucleotide in the sequence can have 9 
# different DBS mutations, so the total counts will be 9x the CDS length, which then get
# normalized to proportions.

# Run from dev refset directory

RefCDS = ces.refset.hg38$RefCDS

# Get the reference dinucleotides for the COSMIC DBS classes, in order
dn_names = sub(">.*", "", cosmic_dbs_classes)
cosmic_dbs_dinuc = unique(dn_names)
reverse_dn_names = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(dn_names)))

process_cds = function(entry) {
  dn_counts = integer(length(cosmic_dbs_classes)) # number of distinct dinucleotide mutations (78)
  # for each transcripts, need to consider each exon and 1 base up/downstream for trinucleotide context
  intervals = entry$intervals_cds # two-column matrix with start/stop coordinates of each cds
  cds_lengths = intervals[,2] - intervals[,1] 
  
  # if transcript is on negative strand, flip exon order
  if (entry$strand == -1) {
    cds_lengths = rev(cds_lengths)
  }
  start = 1
  for (cds_length in cds_lengths) {
    end = start + cds_length
    # xscat and subseq are much more efficient than plain concatenation and subsetting
    seq = Biostrings::xscat(Biostrings::subseq(entry$seq_cds1up, start = start, width = 1), 
                            Biostrings::subseq(entry$seq_cds, start = start, end = end), 
                            Biostrings::subseq(entry$seq_cds1down, start = end, width = 1))
    
    # get trinucleotide counts for cds
    # this function returns full 64-cell table, including counts of 0 (otherwise subsetting below wouldn't work)
    curr_dn_counts = Biostrings::dinucleotideFrequency(seq)
    
    
    cosmic_dinuc_counts = curr_dn_counts[cosmic_dbs_dinuc]
    other_dinuc_names = setdiff(names(curr_dn_counts), cosmic_dbs_dinuc)
    other_dinuc_counts = curr_dn_counts[other_dinuc_names]
    other_names_reversed = as.character(reverseComplement(DNAStringSet(other_dinuc_names)))
    cosmic_dinuc_counts[other_names_reversed] = cosmic_dinuc_counts[other_names_reversed] + other_dinuc_counts
    
    dn_counts = dn_counts + cosmic_dinuc_counts[dn_names]
    
    # next cds sequence starts with next base in the seq_cds sequence
    start = end + 1
  }
  
  # normalize and add trinuc comp to environment
  comp = dn_counts / sum(dn_counts)
  return(unname(comp))
}

dbs_comp = pbapply::pblapply(RefCDS, process_cds)
names(dbs_comp) = names(RefCDS)
dbs_comp = list2env(dbs_comp, parent = emptyenv())
saveRDS(dbs_comp, 'inst/refset/cds_dbs_exposure.rds')
