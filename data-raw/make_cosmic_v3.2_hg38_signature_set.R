library(data.table)
library(cancereffectsizeR)

# Run from ces.refset.hg38 package directory
# Signature definitions downloaded from COSMIC on 09/15/21.
cosmic = fread("inst/extdata/COSMIC_v3.2_SBS_GRCh38.txt")
metadata = fread("inst/extdata/COSMIC_v3.2_signature_metadata.txt")

# column names will be deconstructSigs-style trinuc mutations
dS_muts = cosmic$Type

# drop non-signature columns
cosmic = cosmic[, .SD, .SDcols = patterns("SBS")]
sig_names = colnames(cosmic)

cosmic_df = as.data.frame(t(cosmic))
rownames(cosmic_df) = sig_names
colnames(cosmic_df) = dS_muts

# put columns in canonical order (the order used by deconstructSigs, originally)
deconstructSigs_trinuc_string = getFromNamespace("deconstructSigs_trinuc_string", "cancereffectsizeR")
cosmic_df = cosmic_df[, deconstructSigs_trinuc_string]
signature_set = list(name = "COSMIC v3.2", signatures = cosmic_df, meta = metadata)

# trigger an error if this signature set isn't valid
validate_signature_set(signature_set)

# save in hg38 reference data collection
out_path = "inst/refset/signatures/COSMIC_v3.2_signatures.rds"
saveRDS(signature_set, out_path)


