

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


#This section produces the fcData file, before turning it into a counts file (with 0s), then without 0s
#loads the R ouput
fcData <- read.delim("ryegrass_exon_counts.txt", comment.char="#")

#renames the columns to the sample name, not the extended version
colnames(fcData) <- gsub(".*mapping\\.|\\.bam$", "", colnames(fcData))

#removes the first few columns with other data
counts_unfiltered = fcData[,7:76]
gene_loci = fcData[,1:6]
rownames(counts_unfiltered) <- fcData$Geneid

# current column names
old_names <- colnames(counts_unfiltered)

# identify samples from timepoint 1_5
t15 <- grepl("_1_5$", old_names)

# swap C <-> N only for those samples
new_names <- old_names
new_names[t15] <- gsub("^C", "TEMP", new_names[t15])   # C → TEMP
new_names[t15] <- gsub("^N", "C", new_names[t15])      # N → C
new_names[t15] <- gsub("^TEMP", "N", new_names[t15])   # TEMP → N

# apply the corrected names
colnames(counts_unfiltered) <- new_names

#removes rows with all 0s, as there are a lot of them Rows removed: 26945 - makes the final counts file
filter_nonzero <- function(mat) {
  keep <- rowSums(mat) > 0
  removed <- sum(!keep)
  cat("Rows removed:", removed, "\n")
  mat[keep, ]
}
counts <- filter_nonzero(counts_unfiltered)

#makes a list of the sample names
samples <- colnames(counts)

## TO PUT TREATMENT AND TIMEPOINT AS METADATA FACTORS
# Split on "_", into e.g. part1=C1, part2=5, part3=17
parts <- do.call(rbind, strsplit(samples, "_"))

# Extract treatment and replicate from part1 (C1 > C + 1)
treatment <- substr(parts[,1], 1, 1)
replicate <- substr(parts[,1], 2, nchar(parts[,1]))

# Build metadata frame
coldata <- data.frame(
  sample = samples,
  treatment = factor(treatment, levels = c("C", "N")),
  day = factor(parts[,2]),
  time = factor(parts[,3]),
  replicate = factor(replicate),
  row.names = samples
)

# Reorder timepoints properly: NB+1 → NB+5 → NB+17, day 1 2 3
coldata$time <- factor(coldata$time, levels = c("1", "5", "17"))
coldata$day <- factor(coldata$day,
                      levels = c("1", "2", "3"),
                      labels = c("Day_1", "Day_10", "Day_30"))

#to confirm it worked
head(coldata)

#CPM calculation
# Ensure counts is a numeric matrix
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

# Library sizes (total reads per sample)
lib_sizes <- colSums(counts)

# CPM calculation
cpm <- sweep(counts, 2, lib_sizes / 1e6, FUN = "/")
log_cpm <- log2(cpm + 1)

# ----------------------------------------------------------
# citations() (at the bottom, as requested)
# ----------------------------------------------------------
citations <- function() {
  pkgs <- c("dplyr", "DESeq2", "purrr", "ggplot2", "tidyverse")
  cat("\n================ PACKAGE CITATIONS ================\n")
  for (p in pkgs) {
    cat("\n--- ", p, " ---\n", sep = "")
    print(citation(p))
  }
}

citations()
