library(devtools)

# Run from ces.refet.hg38 dev directory
load_all()

ucsc_info = getFromNamespace(".UCSC_cached_chrom_info", "GenomeInfoDb")[["hg38"]]
ucsc_info = GenomeInfoDb:::.add_ensembl_column(ucsc_info, "hg38")

# Figure out name by inspecting output of getFromNamespace(".NCBI_cached_chrom_info", "GenomeInfoDb")
ncbi_info = getFromNamespace(".NCBI_cached_chrom_info", "GenomeInfoDb")[["GCF_000001405.39"]]

cached_chromInfo = list()
cached_chromInfo[['UCSC']] = list(name = "hg38", value = ucsc_info)
cached_chromInfo[['NCBI']] = list(name = "GCF_000001405.39", value = ncbi_info)
saveRDS(cached_chromInfo, 'inst/refset/cached_chromInfo.rds')
