# # # Global package variables
# #
signature.annotations = read.delim("data/signature_annotations.tsv")
PCAWG.sigs.subset = read.delim("data/PCAWG_sigProfiler_sigs_subset.tsv")

usethis::use_data(signature.annotations, PCAWG.sigs.subset, overwrite = T)
