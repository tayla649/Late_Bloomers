
############################################################
# GLOBAL DESeq2 ANALYSIS (Discovery / Genome-wide)
#
# Provenance / AI assistance statement:
# This script was developed by the author with assistance from
# Microsoft 365 Copilot for code streamlining, naming clarity,
# and minor robustness improvements (e.g., type standardisation,
# safer subsetting, and plot styling). All analytical decisions
# (model design, contrasts, thresholds, and interpretation)
# were made by the author, and results were validated by the author.
############################################################

# Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DESeq2)
  library(Biostrings)
  library(grid)  # for unit()
})

# ==========================================================
# 0) Standardise metadata types (robust subsetting / ordering)
# ==========================================================
# Ensures day/time behave consistently whether numeric/character/factor
coldata <- coldata %>%
  mutate(
    day  = factor(as.character(day),  levels = c("Day_1", "Day_10", "Day_30")),
    time = factor(as.character(time), levels = c("1", "5", "17"))
  )

global_days  <- levels(coldata$day)   # Day_1, Day_10, Day_30
global_times <- levels(coldata$time)  # 1, 5, 17

# ==========================================================
# 1) DESeq2 per Day × Timepoint subset (design: ~ treatment)
# ==========================================================
global_results_list <- list()

for (d in global_days) {
  for (t in global_times) {
    message("GLOBAL: Processing ", d, " | time ", t)
    
    # Subset samples for this day/time
    global_subset_samples <- rownames(
      coldata[coldata$day == d & coldata$time == t, , drop = FALSE]
    )
    
    # Skip if too few samples to run a meaningful model
    if (length(global_subset_samples) < 2) {
      message("GLOBAL: Skipping ", d, " ", t, " — too few samples")
      next
    }
    
    global_counts_sub  <- counts[, global_subset_samples, drop = FALSE]
    global_coldata_sub <- coldata[global_subset_samples, , drop = FALSE]
    
    # Skip if treatment has only one level (needs both C and N)
    if (length(unique(global_coldata_sub$treatment)) < 2) {
      message("GLOBAL: Skipping ", d, " ", t, " — missing C or N samples")
      next
    }
    
    # Create DESeq2 object for this subset (treatment-only model)
    global_dds_sub <- DESeqDataSetFromMatrix(
      countData = global_counts_sub,
      colData   = global_coldata_sub,
      design    = ~ treatment
    )
    
    # Remove genes with zero counts across all samples in this subset
    global_dds_sub <- global_dds_sub[rowSums(counts(global_dds_sub)) > 0, ]
    
    # Fit DESeq2 model
    global_dds_sub <- DESeq(global_dds_sub, quiet = TRUE)
    
    # Extract contrast: N vs C (N relative to C)
    global_res_sub <- results(global_dds_sub, contrast = c("treatment", "N", "C"))
    
    # Store results + identifiers
    global_res_df <- as.data.frame(global_res_sub)
    global_res_df$gene <- rownames(global_res_df)
    global_res_df$day  <- d
    global_res_df$time <- as.character(t)
    
    global_key <- paste0(d, "_", t)
    global_results_list[[global_key]] <- global_res_df
  }
}

# ==========================================================
# 2) Combine results across all Day × Time subsets
# ==========================================================
global_all_tp_DE <- bind_rows(global_results_list) %>%
  mutate(
    day  = factor(as.character(day),  levels = global_days),
    time = factor(as.character(time), levels = c("1", "5", "17"))
  )

# ==========================================================
# 3) Apply global DEG filter (genome-wide per subset)
# ==========================================================
global_filtered_DE <- global_all_tp_DE %>%
  filter(
    !is.na(padj),
    padj < 0.01,
    abs(log2FoldChange) > 1
  )

# ==========================================================
# 4A) Consistent DEGs: same timepoint across ALL days
# ==========================================================
global_DE_same_time <- global_filtered_DE %>%
  group_by(time, gene) %>%
  summarise(n_days = n_distinct(day), .groups = "drop") %>%
  filter(n_days == length(global_days))

global_DE_same_time_full <- global_filtered_DE %>%
  inner_join(global_DE_same_time, by = c("time", "gene"))

# ==========================================================
# 4B) Consistent DEGs: same day across ALL timepoints
# ==========================================================
global_DE_same_day <- global_filtered_DE %>%
  group_by(day, gene) %>%
  summarise(n_times = n_distinct(time), .groups = "drop") %>%
  filter(n_times == length(global_times))

global_DE_same_day_full <- global_filtered_DE %>%
  inner_join(global_DE_same_day, by = c("day", "gene"))

############################################################
# OUTPUT OBJECTS (GLOBAL)
#
# global_all_tp_DE:
#   All genes; DESeq2 results for every Day × Time subset.
#
# global_filtered_DE:
#   Genes passing padj < 0.01 AND |log2FC| > 1 in any subset.
#
# global_DE_same_time_full:
#   Genes passing the DEG filter at the same timepoint in all days,
#   with per-day result rows retained.
#
# global_DE_same_day_full:
#   Genes passing the DEG filter on the same day across all timepoints,
#   with per-timepoint result rows retained.
############################################################

# ==========================================================
# 5) PLOT HELPERS — a priori look & feel (global significance)
# ==========================================================

# --- Mapping provided by you (same-time V3 -> short) ---
same_time_id_to_short <- c(
  "V3.Lp_chr1_0G14956" = "ASL",
  "V3.Lp_chr2_0G19994" = "FLZ2",
  "V3.Lp_chr2_0G694"   = "PPD1",
  "V3.Lp_chr3_0G8870"  = "GI",
  "V3.Lp_chr5_0G1614"  = "CPRF2",
  "V3.Lp_chr5_0G17244" = "PRR95"
)

# Helper: safe filename
safe_filename <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

# Helper (vectorised): get short labels with fallback to V3 IDs
short_label_vec <- function(v3_ids, which_set = c("same_time", "same_day")) {
  which_set <- match.arg(which_set)
  v3_ids <- as.character(v3_ids)
  out <- rep(NA_character_, length(v3_ids))
  
  if (which_set == "same_time") {
    mapped <- unname(same_time_id_to_short[v3_ids])
    out <- mapped
  }
  # If you later add same_day_id_to_short, map here similarly.
  
  # Fallback to V3 ID where mapping is missing
  out[is.na(out)] <- v3_ids
  out
}

# Helper (scalar) for plot titles/labels
short_label_scalar <- function(v3_id, which_set = c("same_time", "same_day")) {
  which_set <- match.arg(which_set)
  lbl <- short_label_vec(v3_id, which_set = which_set)
  lbl[[1]]
}

# Global significance lookup for stars (padj/log2FC from global_all_tp_DE)
global_sig_lookup <- function(gene_id) {
  stopifnot(exists("global_all_tp_DE"))
  global_all_tp_DE %>%
    filter(gene == gene_id) %>%
    mutate(
      day  = factor(as.character(day),  levels = c("Day_1", "Day_10", "Day_30")),
      time = factor(as.character(time), levels = c("1", "5", "17")),
      sig  = !is.na(padj) & padj < 0.01 & abs(log2FoldChange) > 1
    )
}

# Reskinned global plot (matches a priori style)
plot_global_reskinned <- function(
    gene_id,
    out_dir,
    label_set = c("same_time", "same_day"),
    save_png = TRUE,
    save_svg = TRUE,
    png_w = 2.85, png_h = 3.14,
    svg_w = 3.145, svg_h = 2.87,
    dpi = 300
) {
  label_set <- match.arg(label_set)
  
  # Preconditions
  if (!exists("log_cpm")) stop("log_cpm not found. Create log_cpm (log2(CPM+1)) before plotting.")
  if (!gene_id %in% rownames(log_cpm)) {
    message("Skipping gene not in log_cpm: ", gene_id)
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Determine label (short if known for the set; fallback to V3)
  short_name <- short_label_scalar(gene_id, which_set = label_set)
  
  # Data for plotting
  expr <- log_cpm[gene_id, ]
  df <- data.frame(
    sample = names(expr),
    expression = as.numeric(expr)
  ) %>%
    mutate(
      day       = coldata[sample, "day"],
      time      = coldata[sample, "time"],
      treatment = coldata[sample, "treatment"],
      replicate = coldata[sample, "replicate"]
    )
  
  mean_df <- df %>%
    group_by(day, time, treatment) %>%
    summarise(mean_expr = mean(expression), .groups = "drop") %>%
    mutate(
      time = factor(as.character(time), levels = c("1", "5", "17")),
      day  = factor(as.character(day),  levels = c("Day_1", "Day_10", "Day_30"))
    )
  
  # Stars from global significance
  de_info <- global_sig_lookup(gene_id) %>%
    mutate(ypos = Inf)
  
  # In-panel day labels (bottom-right, behind points)
  label_df <- mean_df %>%
    distinct(day) %>%
    mutate(
      label = gsub("Day_", "Day ", as.character(day)),
      x = Inf, y = -Inf
    )
  
  # Facet layer compatible across ggplot2 versions
  facet_layer <- NULL
  if (utils::packageVersion("ggplot2") >= "3.4.0") {
    facet_layer <- facet_wrap(
      ~ day, ncol = 1,
      scales = "fixed",
      labeller = labeller(day = function(x) gsub("Day_", "Day ", x)),
      strip.position = "top",
      axes = "margins"
    )
  } else {
    facet_layer <- facet_grid(
      day ~ .,
      labeller = labeller(day = function(x) gsub("Day_", "Day ", x))
    )
  }
  
  # Y label (italic(short) expression (log2CPM)); fallback prints V3 if unmapped
  ylab_expr <- parse(text = sprintf("italic('%s')~' expression ('*log[2]*'CPM)'", short_name))
  
  base_size <- 11 * 1
  
  p <- ggplot() +
    # In-panel day labels
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.6,
      color = "black",
      alpha = 0.15,
      size = 6,
      fontface = "bold"
    ) +
    geom_point(
      data = df,
      aes(x = time, y = expression, color = treatment),
      size = 0.5,
      alpha = 0.8,
      position = position_jitter(width = 0.05)
    ) +
    geom_line(
      data = mean_df,
      aes(x = time, y = mean_expr, color = treatment, group = treatment),
      linewidth = 0.5
    ) +
    geom_text(
      data = de_info,
      aes(x = time, y = ypos, label = ifelse(sig, "*", "")),
      vjust = 1,
      position = position_nudge(x = 0.1),
      color = "black",
      size = 5,
      fontface = "bold"
    ) +
    facet_layer +
    theme_bw(base_size = base_size) +
    labs(
      x = "Timepoint",
      y = ylab_expr
    ) +
    scale_x_discrete(labels = c("1" = "NB+1h", "5" = "NB+5h", "17" = "NB+17h")) +
    coord_cartesian(ylim = c(0, NA)) +
    theme(
      legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = unit(1.2, "lines"),
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.8),
      panel.spacing.y = unit(0.15, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  print(p)
  
  # Save (filenames use short when available)
  base_file <- safe_filename(short_name)
  if (save_png) {
    ggsave(
      filename = file.path(out_dir, paste0(base_file, "_expression_plot.png")),
      plot = p, width = png_w, height = png_h, dpi = dpi, bg = "white", limitsize = FALSE
    )
  }
  if (save_svg) {
    ggsave(
      filename = file.path(out_dir, paste0(base_file, "_expression_plot.svg")),
      plot = p, width = svg_w, height = svg_h, bg = "white", limitsize = FALSE
    )
  }
  
  invisible(p)
}

# ==========================================================
# 5A) PLOTS: SAME-TIME consistent DE genes (reskinned)
# ==========================================================
if (!exists("log_cpm")) {
  stop("log_cpm not found. Create log_cpm (log2(CPM+1)) before plotting.")
}
if (!exists("global_DE_same_time") || nrow(global_DE_same_time) == 0) {
  stop("global_DE_same_time not found or empty. Run the global consistency section first.")
}

dir_same_time <- "global_same_time_gene_plots_reskinned"
dir.create(dir_same_time, showWarnings = FALSE, recursive = TRUE)

unique_same_time_genes <- unique(global_DE_same_time$gene)
message("Plotting (reskinned) ", length(unique_same_time_genes), " same-time genes into: ", dir_same_time)

for (gid in unique_same_time_genes) {
  plot_global_reskinned(
    gene_id   = gid,
    out_dir   = dir_same_time,
    label_set = "same_time"   # use your short-name mapping when present
  )
}

# Membership with short names (for convenience)
same_time_membership <- global_DE_same_time %>%
  group_by(gene) %>%
  summarise(timepoints = paste(sort(unique(time)), collapse = ","), .groups = "drop") %>%
  mutate(shortname = short_label_vec(gene, "same_time")) %>%
  select(v3_id = gene, shortname, timepoints)

write.csv(
  same_time_membership,
  file = file.path(dir_same_time, "same_time_gene_membership_with_shortnames.csv"),
  row.names = FALSE
)

# ==========================================================
# 5B) PLOTS: SAME-DAY consistent DE genes (reskinned)
# ==========================================================
if (!exists("global_DE_same_day") || nrow(global_DE_same_day) == 0) {
  stop("global_DE_same_day not found or empty. Run the global consistency section first.")
}

dir_same_day <- "global_same_day_gene_plots_reskinned"
dir.create(dir_same_day, showWarnings = FALSE, recursive = TRUE)

unique_same_day_genes <- unique(global_DE_same_day$gene)
message("Plotting (reskinned) ", length(unique_same_day_genes), " same-day genes into: ", dir_same_day)

for (gid in unique_same_day_genes) {
  plot_global_reskinned(
    gene_id   = gid,
    out_dir   = dir_same_day,
    label_set = "same_day"   # falls back to V3 unless you add a same-day map
  )
}

# Membership with short names (falls back to V3 in shortname column if unknown)
same_day_membership <- global_DE_same_day %>%
  group_by(gene) %>%
  summarise(days = paste(sort(unique(day)), collapse = ","), .groups = "drop") %>%
  mutate(shortname = short_label_vec(gene, "same_day")) %>%
  select(v3_id = gene, shortname, days)

write.csv(
  same_day_membership,
  file = file.path(dir_same_day, "same_day_gene_membership_with_shortnames.csv"),
  row.names = FALSE
)

# ==========================================================
# 6) FASTA EXPORT: sequences for consistent-gene sets
# ==========================================================
# Exports protein sequences for genes in:
# - global_DE_same_time_full (same-time consistent)
# - global_DE_same_day_full  (same-day consistent)
#
# Matching strategy:
# - Uses the DESeq2 "gene" IDs directly.
# - FASTA headers are lightly cleaned by:
#     (i) keeping only the first token before whitespace
#     (ii) removing a trailing isoform suffix like ".1", ".2"
# - After cleaning, sequences are extracted by exact ID matching.
# ==========================================================

# 1) IDs to retrieve (using DESeq2 gene IDs)
same_time_lp <- unique(global_DE_same_time_full$gene)
same_day_lp  <- unique(global_DE_same_day_full$gene)

# 2) Load L. perenne proteome (AA FASTA)
#    Ensure "lp_proteome.fasta" is present in the working directory
lp_proteome <- readAAStringSet("lp_proteome.fasta")

# 3) Clean FASTA headers for matching
names(lp_proteome) <- sub("\\s.*$", "", names(lp_proteome))    # keep first token
names(lp_proteome) <- sub("\\.[0-9]+$", "", names(lp_proteome)) # drop ".1" etc

# 4) Extract sequences (exact match after header cleaning)
same_time_lp_seqs <- lp_proteome[names(lp_proteome) %in% same_time_lp]
same_day_lp_seqs  <- lp_proteome[names(lp_proteome) %in% same_day_lp]

# 5) Report missing IDs
missing_time <- setdiff(same_time_lp, names(same_time_lp_seqs))
missing_day  <- setdiff(same_day_lp,  names(same_day_lp_seqs))

if (length(missing_time) > 0) {
  cat("\nWARNING: same-time genes missing from proteome (after header cleaning):\n",
      paste(missing_time, collapse = ", "), "\n")
}
if (length(missing_day) > 0) {
  cat("\nWARNING: same-day genes missing from proteome (after header cleaning):\n",
      paste(missing_day, collapse = ", "), "\n")
}

# 6) Write FASTA outputs
writeXStringSet(same_time_lp_seqs, "same_time_consistent_Lp.fa")
writeXStringSet(same_day_lp_seqs,  "same_day_consistent_Lp.fa")

cat("\nSaved:\n  same_time_consistent_Lp.fa (", length(same_time_lp_seqs), " seqs)\n",
    "  same_day_consistent_Lp.fa  (", length(same_day_lp_seqs),  " seqs)\n", sep = "")

# 7) Optional inspection prints
cat("\n================ SAME-TIME AMINO ACID SEQUENCES =================\n")
print(same_time_lp_seqs)

cat("\n================ SAME-DAY AMINO ACID SEQUENCES ==================\n")
print(same_day_lp_seqs)

cat("\n================ PREVIEW (first 60 aa) ==========================\n")
for (nm in names(same_time_lp_seqs)) {
  s <- as.character(same_time_lp_seqs[[nm]])
  cat("\n>", nm, "\n", substr(s, 1, 60), "...\n", sep = "")
}

# ==========================================================
# citations() (package citation output)
# ==========================================================
citations <- function() {
  pkgs <- c("dplyr", "tidyr", "ggplot2", "DESeq2", "Biostrings")
  cat("\n================ PACKAGE CITATIONS ================\n")
  for (p in pkgs) {
    cat("\n--- ", p, " ---\n", sep = "")
    print(citation(p))
  }
}

# Call at end if desired
# citations()
