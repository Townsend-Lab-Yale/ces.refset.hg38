# Run from refset dev directory.
library(data.table)
library(devtools)

load_all()

# Load blacklist, as downloaded from UCSC Table Browser.
bl = rtracklayer::import('inst/extdata/encBlacklistv2_hg38.bed.gz')
ranges = clean_granges_for_cesa(refset_env = ces.refset.hg38, gr = bl, reduce_sort_strip = FALSE)
saveRDS(ranges, 'inst/refset/ENCODE_blacklist_v2.rds')
