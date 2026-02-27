
############################################################
# Pairwise LD Scan (EM haplotype inference + permutation test)
# Tabular-only output (no plots / no haplotype matrix printing)
#
# Loci included:
#   CO, FT3, FTL9, VRN1, VRN2a, VRN2b, PHYB, PPD1
#
# Provenance / AI assistance statement:
# This script was developed by the author with assistance from
# Microsoft 365 Copilot for code formatting, naming clarity,
# and minor robustness/cleanliness improvements (e.g., simplifying
# to tabular-only output and adding a BH-corrected significance flag).
# All analytical decisions (LD statistic, permutation strategy,
# multiple testing method, and interpretation) were made by the author.
# Results should be validated by the author on the target dataset.
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(haplo.stats)
})

############################################################
# WORKFLOW SUMMARY
#
# Input:
#   - genotypes.csv with columns {LOCUS_H1, LOCUS_H2} for each locus.
#
# Steps:
#   1) Read CSV and aggressively clean allele strings.
#   2) For each locus pair:
#        - integer-encode alleles (0 = missing for haplo.em)
#        - infer haplotype frequencies via EM (haplo.em)
#        - compute LD evidence:
#            * G = likelihood-ratio statistic vs independence
#            * V = Cramér’s V effect size
#        - permutation p-value:
#            * shuffle locus2 genotypes across individuals B times
#   3) Adjust permutation p-values across pairs with BH FDR (q_BH).
#   4) Add linked flag: linked = TRUE if q_BH < 0.05 else FALSE.
#
# Output:
#   - Printed summary table (and optional CSV) with:
#       Locus1, Locus2, G, V, p_perm, q_BH, linked
############################################################

# ==========================================================
# USER SETTINGS
# ==========================================================
infile           <- "genotypes.csv"   # <--- update path if needed
B_permutations   <- 5000              # increase for more precise p-values
rng_seed         <- 42                # reproducibility

save_results_csv <- TRUE
out_csv          <- "pairwise_ld_results_with_PHYB_PPD1.csv"

# Loci to include (must have *_H1 and *_H2 columns in the CSV)
loci <- c("CO","FT3","FTL9","VRN1","VRN2a","VRN2b","PHYB","PPD1")

# ==========================================================
# 1) READ & CLEAN DATA
# ==========================================================
df_raw <- read_csv(
  infile,
  col_types = cols(.default = col_character()),
  na = c("", "NA", "N/A", ".", "na", "n/a")
)

# Strip UTF-8 BOM if present
names(df_raw) <- gsub("\uFEFF", "", names(df_raw))

# Check required columns
required_cols <- as.vector(rbind(paste0(loci, "_H1"), paste0(loci, "_H2")))
missing_cols  <- setdiff(required_cols, names(df_raw))
if (length(missing_cols)) {
  stop("Missing columns in CSV: ", paste(missing_cols, collapse = ", "))
}

# Aggressive allele cleaner:
# - trims whitespace (including NBSP), removes internal whitespace
# - keeps only letters A–Z
# - uppercases
# - blank -> NA
clean_allele <- function(x) {
  x <- gsub("\u00A0", " ", x, fixed = TRUE)  # NBSP -> space
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("[^A-Za-z]", "", x)
  x <- toupper(x)
  x[x == ""] <- NA
  x
}

df <- df_raw
for (g in loci) {
  df[[paste0(g, "_H1")]] <- clean_allele(df[[paste0(g, "_H1")]])
  df[[paste0(g, "_H2")]] <- clean_allele(df[[paste0(g, "_H2")]])
}

# Optional: brief allele summary (comment out if you want absolutely minimal output)
cat("Locus summary after cleaning:\n")
for (g in loci) {
  a1 <- df[[paste0(g, "_H1")]]
  a2 <- df[[paste0(g, "_H2")]]
  alleles <- sort(unique(na.omit(c(a1, a2))))
  cat(sprintf("  %-5s alleles={%s}\n", g, paste(alleles, collapse=",")))
}

# ==========================================================
# 2) UTILITIES
# ==========================================================

# Integer-encode alleles per locus with shared levels across H1/H2
# (0 = missing for haplo.em)
encode_locus <- function(h1, h2) {
  lvls <- sort(unique(na.omit(c(h1, h2))))
  a1 <- match(h1, lvls)
  a2 <- match(h2, lvls)
  a1[is.na(a1)] <- 0L
  a2[is.na(a2)] <- 0L
  list(a1 = as.integer(a1), a2 = as.integer(a2), levels = lvls)
}

# EM-based LD for a pair of loci + permutation p-value
# Returns tabular essentials only (no plotting, no matrix output).
ld_two_loci_em <- function(df, locus1, locus2, B = 2000, seed = 42) {
  
  enc1 <- encode_locus(df[[paste0(locus1, "_H1")]], df[[paste0(locus1, "_H2")]])
  enc2 <- encode_locus(df[[paste0(locus2, "_H1")]], df[[paste0(locus2, "_H2")]])
  
  if (length(enc1$levels) == 0 || length(enc2$levels) == 0) {
    return(list(Locus1 = locus1, Locus2 = locus2,
                G = NA_real_, V = NA_real_, p_perm = NA_real_,
                note = "all-missing"))
  }
  
  # Monomorphic locus -> trivial LD
  if (length(enc1$levels) == 1 || length(enc2$levels) == 1) {
    return(list(Locus1 = locus1, Locus2 = locus2,
                G = 0, V = 0, p_perm = 1,
                note = "monomorphic"))
  }
  
  # EM haplotype inference (unphased)
  geno <- cbind(enc1$a1, enc1$a2, enc2$a1, enc2$a2)
  em <- haplo.em(geno = geno, locus.label = c(locus1, locus2), miss.val = 0)
  
  # Rebuild haplotype frequency matrix from EM output
  r <- length(enc1$levels)
  c <- length(enc2$levels)
  hap_mat <- matrix(0, nrow = r, ncol = c)
  
  ai <- as.integer(em$haplotype[[1]])
  bi <- as.integer(em$haplotype[[2]])
  
  for (i in seq_along(em$hap.prob)) {
    if (!is.na(ai[i]) && !is.na(bi[i]) && ai[i] >= 1 && ai[i] <= r && bi[i] >= 1 && bi[i] <= c) {
      hap_mat[ai[i], bi[i]] <- hap_mat[ai[i], bi[i]] + em$hap.prob[i]
    }
  }
  hap_mat <- hap_mat / sum(hap_mat)
  
  # Likelihood-ratio G statistic vs independence + Cramér's V
  N_gam <- 2 * nrow(df)
  eps <- 1e-12
  
  P  <- hap_mat + eps; P <- P / sum(P)
  pi <- rowSums(P)
  pj <- colSums(P)
  P0 <- outer(pi, pj) + eps
  
  G <- 2 * N_gam * sum(P * log(P / P0))
  V <- sqrt((G / N_gam) / min(nrow(P) - 1, ncol(P) - 1))
  
  # Permutation p-value: shuffle locus2 across individuals
  set.seed(seed)
  
  perm_once <- function() {
    ix <- sample.int(nrow(df))
    em_p <- haplo.em(
      geno = cbind(enc1$a1, enc1$a2, enc2$a1[ix], enc2$a2[ix]),
      locus.label = c(locus1, locus2),
      miss.val = 0
    )
    
    ai_p <- as.integer(em_p$haplotype[[1]])
    bi_p <- as.integer(em_p$haplotype[[2]])
    hm_p <- matrix(0, nrow = r, ncol = c)
    
    for (k in seq_along(em_p$hap.prob)) {
      if (!is.na(ai_p[k]) && !is.na(bi_p[k]) && ai_p[k] >= 1 && ai_p[k] <= r && bi_p[k] >= 1 && bi_p[k] <= c) {
        hm_p[ai_p[k], bi_p[k]] <- hm_p[ai_p[k], bi_p[k]] + em_p$hap.prob[k]
      }
    }
    hm_p <- hm_p / sum(hm_p)
    
    Pp  <- hm_p + eps; Pp <- Pp / sum(Pp)
    pi_ <- rowSums(Pp)
    pj_ <- colSums(Pp)
    P0p <- outer(pi_, pj_) + eps
    
    2 * N_gam * sum(Pp * log(Pp / P0p))
  }
  
  G_perm <- replicate(B, perm_once())
  p_perm <- (sum(G_perm >= G) + 1) / (B + 1)
  
  list(Locus1 = locus1, Locus2 = locus2, G = G, V = V, p_perm = p_perm, note = NA_character_)
}

# ==========================================================
# 3) RUN ALL PAIRS + SUMMARY TABLE (BH-corrected linked flag)
# ==========================================================
pairs <- t(combn(loci, 2))
cat("\n--- Running all pairs (", nrow(pairs), " tests) ---\n", sep = "")

res_list <- apply(
  pairs, 1,
  function(x) ld_two_loci_em(df, x[1], x[2], B = B_permutations, seed = rng_seed)
)

res_tbl <- bind_rows(lapply(res_list, as_tibble)) %>%
  mutate(
    q_BH   = p.adjust(p_perm, method = "BH"),
    linked = !is.na(q_BH) & q_BH < 0.05
  ) %>%
  arrange(q_BH, p_perm) %>%
  select(Locus1, Locus2, G, V, p_perm, q_BH, linked)

cat("\n--- LD summary table (sorted by q_BH, then p_perm) ---\n")
print(res_tbl, n = nrow(res_tbl))

if (save_results_csv) {
  write_csv(res_tbl, out_csv)
  cat("\nSaved results to:", out_csv, "\n")
}

# ==========================================================
# citations() and session info (reproducibility)
# ==========================================================
citations <- function() {
  pkgs <- c("readr", "dplyr", "tibble", "haplo.stats")
  cat("\n================ PACKAGE CITATIONS ================\n")
  for (p in pkgs) {
    cat("\n--- ", p, " ---\n", sep = "")
    print(citation(p))
  }
}

citations()

cat("\n================ sessionInfo() ================\n")
print(sessionInfo())
