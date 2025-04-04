library(devtools)

# Run from ces.refet.hg38 dev directory
load_all()

library(BSgenome.Hsapiens.UCSC.hg38)
bsg = BSgenome::getBSgenome(ces.refset.hg38::ces.refset.hg38$genome)

ucsc_info = getFromNamespace(".UCSC_cached_chrom_info", "GenomeInfoDb")[["hg38"]]
ucsc_info = GenomeInfoDb:::.add_ensembl_column(ucsc_info, "hg38")

# Figure out name by inspecting output of names(getFromNamespace(".NCBI_cached_chrom_info", "GenomeInfoDb")) and choosing latest version
curr_ncbi_name = "GCF_000001405.40"
ncbi_info = getFromNamespace(".NCBI_cached_chrom_info", "GenomeInfoDb")[[curr_ncbi_name]]

cached_chromInfo = list()
cached_chromInfo[['UCSC']] = list(name = "hg38", value = ucsc_info)
cached_chromInfo[['NCBI']] = list(name = "GCF_000001405.40", value = ncbi_info)
saveRDS(cached_chromInfo, 'inst/refset/cached_chromInfo.rds')
