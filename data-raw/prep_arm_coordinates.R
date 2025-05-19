# Run from refset dev directory.
library(data.table)
library(devtools)

load_all()

# Get arm starts/ends from UCSC cytoband info.
cytoband_info = fread('inst/extdata/cytobands_hg38.txt')
setnames(cytoband_info, c('#chrom', 'chromStart', 'chromEnd'), c('chr', 'start', 'end'))
cytoband_info[, chr := sub('^chr', '', chr)]
cytoband_info[, arm := substr(name, 1, 1)]
cytoband_info = cytoband_info[chr %in% c(1:22, 'X', 'Y')]
stopifnot(identical(unique(cytoband_info$arm), c('p', 'q')))

# Text data from UCSC table browser is 0-based, half-open. Convert to 1-based closed.
cytoband_info[, start := start + 1]

# Get arm boundaries and centromeric regions (note centromeric regions include p and q segments)
coordinates = cytoband_info[, .(p_start = min(start), p_end = .SD[arm == 'p', max(end)],
                                q_start = .SD[arm == 'q', min(start)], q_end = max(end),
                                cen_start = .SD[arm == 'p' & gieStain == 'acen', min(start)],
                                cen_end = .SD[arm == 'q' & gieStain == 'acen', max(end)]), by = 'chr']
stopifnot(coordinates[, all(q_end > q_start & cen_end > q_start & 
                              cen_end > cen_start & p_end > cen_start & p_start == 1)])

saveRDS(coordinates, 'inst/refset/arm_coordinates.rds')



