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
