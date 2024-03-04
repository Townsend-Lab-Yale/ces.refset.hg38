library(data.table)
library(cancereffectsizeR)

# Run from ces.refset.hg38 package directory
# Signature definitions downloaded from COSMIC on 09/15/21.
cosmic = fread("inst/extdata/COSMIC_v3.4_SBS_GRCh38.txt")

old_meta = fread("inst/extdata/COSMIC_v3.2_signature_metadata.txt")

current_signatures = setdiff(names(cosmic), 'Type')
sbs_etiologies = cosmic_signature_info() # From cancereffectsizeR's data

stopifnot(all(current_signatures %in% sbs_etiologies$name))
new_meta = sbs_etiologies[current_signatures, on = 'name']

new_signatures = setdiff(new_meta$name, old_meta$Signature)
# new_meta[new_signatures, on = 'name']
## New signatures consist of 22 getting split to 22a,b; 40 getting split to 40a,b,c; and 95-99. None of the new
## signatures are hypermutation signatures. SBS95 is a new artifact signature. All other information can be merged in from v3.2.

new_meta[old_meta, let(Likely_Artifact = Likely_Artifact, Exome_Min = Exome_Min, 
                       Genome_Min = Genome_Min), on = c(name = 'Signature')]
new_meta[name %in% new_signatures, Likely_Artifact := FALSE]
new_meta[name == 'SBS95', Likely_Artifact := TRUE]


# column names will be deconstructSigs-style trinuc mutations
dS_muts = cosmic$Type

# drop non-signature columns
cosmic = cosmic[, .SD, .SDcols = patterns("SBS")]
sig_names = names(cosmic)

cosmic_df = as.data.frame(t(cosmic))
rownames(cosmic_df) = sig_names
colnames(cosmic_df) = dS_muts

# put columns in canonical order (the order used by deconstructSigs, originally)
deconstructSigs_trinuc_string = getFromNamespace("deconstructSigs_trinuc_string", "cancereffectsizeR")
cosmic_df = cosmic_df[, deconstructSigs_trinuc_string]
signature_set = list(name = "COSMIC v3.4", signatures = cosmic_df, meta = new_meta)

# trigger an error if this signature set isn't valid
validate_signature_set(signature_set)

# save in hg38 reference data collection
out_path = "inst/refset/signatures/COSMIC_v3.4_signatures.rds"
saveRDS(signature_set, out_path)


