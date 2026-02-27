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

# ============================================
# LD test for PHYB × VRN2b (EM + permutation)
# ============================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(haplo.stats)
  library(ggplot2)
})

# ---- Settings ----
infile <- "genotypes.csv"   # <- change if needed
B_permutations <- 5000      # increase for more precise p-values
rng_seed <- 42
save_plot_png <- TRUE

# ---- 1) Read & clean data ----
df_raw <- read_csv(
  infile,
  col_types = cols(.default = col_character()),
  na = c("", "NA", "N/A", ".", "na", "n/a")
)
# Strip BOM if present
names(df_raw) <- gsub("\uFEFF", "", names(df_raw))

# Check required columns
req <- c("VRN2b_H1","VRN2b_H2","PHYB_H1","PHYB_H2")
missing_cols <- setdiff(req, names(df_raw))
if (length(missing_cols))
  stop("Missing columns in CSV: ", paste(missing_cols, collapse=", "))

# Aggressive allele cleaner: keep only letters, trim spaces/NBSP/punct, uppercase, blanks->NA
clean_allele <- function(x) {
  x <- gsub("\u00A0", " ", x, fixed = TRUE)  # NBSP -> space
  x <- trimws(x)
  x <- gsub("\\s+", "", x)                    # remove internal whitespace
  x <- gsub("[^A-Za-z]", "", x)               # keep only A–Z
  x <- toupper(x)
  x[x == ""] <- NA
  x
}

df <- df_raw %>%
  mutate(
    VRN2b_H1 = clean_allele(VRN2b_H1),
    VRN2b_H2 = clean_allele(VRN2b_H2),
    PHYB_H1  = clean_allele(PHYB_H1),
    PHYB_H2  = clean_allele(PHYB_H2)
  )

cat("Alleles after cleaning:\n")
cat("  VRN2b: ", paste(sort(unique(na.omit(c(df$VRN2b_H1, df$VRN2b_H2)))), collapse=","), "\n", sep="")
cat("  PHYB : ", paste(sort(unique(na.omit(c(df$PHYB_H1,  df$PHYB_H2)))),  collapse=","), "\n", sep="")

# ---- 2) Utilities ----
# Integer-encode per locus with shared levels across H1/H2 (0 = missing for haplo.em)
encode_locus <- function(h1, h2) {
  lvls <- sort(unique(na.omit(c(h1, h2))))
  a1 <- match(h1, lvls); a2 <- match(h2, lvls)
  a1[is.na(a1)] <- 0L; a2[is.na(a2)] <- 0L
  list(a1 = as.integer(a1), a2 = as.integer(a2), levels = lvls)
}

# Core: EM-based LD for a pair + permutation p-value
ld_two_loci_em <- function(df, locus1, locus2, B = 2000, seed = 42) {
  # locus1/locus2 must be "PHYB" and "VRN2b" in this script
  enc1 <- encode_locus(df[[paste0(locus1,"_H1")]], df[[paste0(locus1,"_H2")]])
  enc2 <- encode_locus(df[[paste0(locus2,"_H1")]], df[[paste0(locus2,"_H2")]])
  
  if (length(enc1$levels) == 0 || length(enc2$levels) == 0)
    stop("All-missing locus at ", locus1, " or ", locus2, " after cleaning.")
  
  # If monomorphic at either locus, LD is trivial
  if (length(enc1$levels) == 1 || length(enc2$levels) == 1) {
    hm <- matrix(1, nrow = length(enc1$levels), ncol = length(enc2$levels),
                 dimnames = list(enc1$levels, enc2$levels))
    storage.mode(hm) <- "double"
    return(list(locus1=locus1, locus2=locus2, hap_mat=hm, G=0, V=0, p_perm=1, note="monomorphic"))
  }
  
  # EM haplotype inference (unphased)
  geno <- cbind(enc1$a1, enc1$a2, enc2$a1, enc2$a2)
  em <- haplo.em(geno = geno, locus.label = c(locus1, locus2), miss.val = 0)
  
  # Rebuild haplotype matrix from em$haplotype indices and our level lists
  r <- length(enc1$levels); c <- length(enc2$levels)
  hap_mat <- matrix(0, nrow = r, ncol = c)
  ai <- as.integer(em$haplotype[[1]])
  bi <- as.integer(em$haplotype[[2]])
  for (i in seq_along(em$hap.prob)) {
    if (!is.na(ai[i]) && !is.na(bi[i]) && ai[i]>=1 && ai[i]<=r && bi[i]>=1 && bi[i]<=c) {
      hap_mat[ai[i], bi[i]] <- hap_mat[ai[i], bi[i]] + em$hap.prob[i]
    }
  }
  hap_mat <- hap_mat / sum(hap_mat)
  rownames(hap_mat) <- enc1$levels
  colnames(hap_mat) <- enc2$levels
  storage.mode(hap_mat) <- "double"
  
  # Likelihood-ratio G statistic vs independence + Cramér's V
  N_gam <- 2 * nrow(df)       # total "gametes"
  eps <- 1e-12
  P <- hap_mat + eps; P <- P / sum(P)
  pi <- rowSums(P); pj <- colSums(P)
  P0 <- outer(pi, pj) + eps
  G <- 2 * N_gam * sum(P * log(P / P0))
  V <- sqrt((G / N_gam) / min(nrow(P)-1, ncol(P)-1))
  
  # Permutation p-value (shuffle second locus across individuals)
  set.seed(seed)
  perm_once <- function() {
    ix <- sample.int(nrow(df))
    em_p <- haplo.em(
      geno = cbind(enc1$a1, enc1$a2, enc2$a1[ix], enc2$a2[ix]),
      locus.label = c(locus1, locus2), miss.val = 0
    )
    ai_p <- as.integer(em_p$haplotype[[1]])
    bi_p <- as.integer(em_p$haplotype[[2]])
    hm_p <- matrix(0, nrow = r, ncol = c)
    for (k in seq_along(em_p$hap.prob)) {
      if (!is.na(ai_p[k]) && !is.na(bi_p[k]) && ai_p[k]>=1 && ai_p[k]<=r && bi_p[k]>=1 && bi_p[k]<=c) {
        hm_p[ai_p[k], bi_p[k]] <- hm_p[ai_p[k], bi_p[k]] + em_p$hap.prob[k]
      }
    }
    hm_p <- hm_p / sum(hm_p)
    Pp <- hm_p + eps; Pp <- Pp / sum(Pp)
    pi_ <- rowSums(Pp); pj_ <- colSums(Pp)
    P0p <- outer(pi_, pj_) + eps
    2 * N_gam * sum(Pp * log(Pp / P0p))
  }
  G_perm <- replicate(B, perm_once())
  p_perm <- (sum(G_perm >= G) + 1) / (B + 1)
  
  list(locus1=locus1, locus2=locus2, hap_mat=hap_mat, G=G, V=V, p_perm=p_perm, note=NA)
}

# ---- 3) Run PHYB × VRN2b ----
cat("\n--- PHYB × VRN2b ---\n")
res <- ld_two_loci_em(df, "PHYB", "VRN2b", B = B_permutations, seed = rng_seed)

cat("G =", round(res$G, 2),
    "| Cramer's V =", round(res$V, 3),
    "| permutation p =", signif(res$p_perm, 4), "\n")
cat("\nInferred haplotype frequencies (PHYB alleles as rows × VRN2b alleles as cols):\n")
print(round(res$hap_mat, 4))

# ---- 4) Plot heatmap ----
hm <- res$hap_mat; storage.mode(hm) <- "double"
p <- as.data.frame(hm) %>%
  rownames_to_column("PHYB") %>%
  pivot_longer(-PHYB, names_to = "VRN2b", values_to = "Freq") %>%
  ggplot(aes(VRN2b, PHYB, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq, accuracy = 0.1)),
            color = "white", size = 3) +
  scale_fill_viridis_c(name = "Freq") +
  labs(title = "Inferred haplotype frequencies: PHYB × VRN2b",
       x = "VRN2b allele", y = "PHYB allele") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
if (save_plot_png) {
  ggsave("haplo_heatmap_PHYB_x_VRN2b.png", plot = p, width = 7, height = 5, dpi = 300)
  cat("Saved heatmap to: haplo_heatmap_PHYB_x_VRN2b.png\n")
}

# ---- 5) (Optional) Backup: genotype×genotype association (Monte-Carlo χ²) ----
do_geno_test <- FALSE
if (do_geno_test) {
  gtab <- df %>%
    mutate(
      PHYB_geno  = paste0("PHYB_",  pmin(PHYB_H1,  PHYB_H2),  pmax(PHYB_H1,  PHYB_H2)),
      VRN2b_geno = paste0("VRN2b_", pmin(VRN2b_H1, VRN2b_H2), pmax(VRN2b_H1, VRN2b_H2))
    ) %>%
    count(PHYB_geno, VRN2b_geno) %>%
    pivot_wider(names_from = VRN2b_geno, values_from = n, values_fill = 0)
  
  mat <- as.matrix(column_to_rownames(gtab, "PHYB_geno"))
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]
  
  set.seed(rng_seed)
  chi_mc <- chisq.test(mat, simulate.p.value = TRUE, B = 100000)
  n <- sum(mat)
  Vc <- sqrt(as.numeric(chi_mc$statistic) / (n * min(nrow(mat)-1, ncol(mat)-1)))
  
  cat("\n[Genotype×Genotype] Monte-Carlo χ² (PHYB×VRN2b):\n")
  print(chi_mc)
  cat("Cramer's V =", round(Vc, 3), "\n")
}
