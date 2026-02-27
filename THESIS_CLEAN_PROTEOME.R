# ----------------------------------------------------------
# Provenance / AI assistance statement:
# ----------------------------------------------------------
# This script was developed by the author with assistance from
# Microsoft 365 Copilot for code streamlining, naming clarity,
# and minor robustness improvements (e.g., type standardisation,
# safer subsetting, and improved annotation spacing). All analytical
# decisions (candidate gene list, contrasts, thresholds, and
# interpretation) were made by the author, and results were validated
# by the author.

library(seqinr)

# read FASTA
seqs <- read.fasta("Lolium_2.6.1_V3_all_transcripts_PROT.fasta",
                   seqtype = "AA", as.string = TRUE)

# remove dots from sequences only
seqs_clean <- lapply(seqs, function(s) gsub("\\.", "", s))

# write cleaned FASTA
write.fasta(seqs_clean,
            names(seqs),
            file.out = "lp_proteome.fasta")
