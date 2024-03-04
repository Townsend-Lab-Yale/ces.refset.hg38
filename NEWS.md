# ces.refset.hg38 1.3.0
* Annotations are now derived from the Gencode v45 basic annotation GTF.
* More information about coding transcripts is now available, including which are canonical (MANE) transcripts.
* Coding and noncoding regions of protein-coding and noncoding transcripts are viewable in a transcript table.
* Added COSMIC SBS signature definitions v3.4.

# ces.refset.hg38 1.2.2
* A tweak to ensure package can load without an internet connection.

# ces.refset.hg38 1.2.1
* Added BSgenome.Hsapiens.UCSC.hg38 dependency to DESCRIPTION.

# ces.refset.hg38 1.2.0
* Corrected gene trinculeotide content proportions. Incorrect values were caused by a bug in cancereffectsizeR::create_refset() that was fixed in cancereffectsizeR v2.4.0.
* Removed some extraneous transcripts from RefCDS (transcript where all exons were already included in other, longer transcripts).

# ces.refset.hg38 1.1.0
* Modest change to transcript set selection to prioritize consensus coding regions (CCDS) and reduce the number of relatively redundant transcripts per gene. (All transcripts still taken from Gencode basic gene annotations, release 38.)
* Default exome expanded to cover all exon regions of protein-coding genes (as opposed to just CDS regions).

# ces.refset.hg38 1.0.0
* Initial release of hg38 reference data package.
* Uses Gencode "basic" gene annotations (release 38) and COSMIC v3.2 signature definitions.
* Regional (gene) mutation rate covariates data are taken from ces.refset.hg19 (updated covariates may come in a future release).
* Additional genomic site annotations from RepeatMasker (via UCSC Genome Browser), COSMIC, and gnomAD (via Annovar).


