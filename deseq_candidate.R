
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

library(dplyr)
library(DESeq2)
library(purrr)
library(ggplot2)
library(tidyverse)

gene_map <- c(
  FTL9       = "V3.Lp_chr1_0G17602",
  ELF3       = "V3.Lp_chr1_0G19888",
  PPD1       = "V3.Lp_chr2_0G694",
  GI         = "V3.Lp_chr3_0G8870",
  LUX        = "V3.Lp_chr3_0G26392",
  VRN1       = "V3.Lp_chr4_0G5312",
  PHYC       = "V3.Lp_chr4_0G5334",
  ELF4       = "V3.Lp_chr4_0G14614",
  LNK1_like  = "V3.Lp_chr4_0G18178",
  COL13_like = "V3.Lp_chr4_0G20406",
  PHYB       = "V3.Lp_chr4_0G21912",
  VRN2b      = "V3.Lp_chr4_0G28974",
  VRN2a      = "V3.Lp_chr4_0G29098",
  FT3        = "V3.Lp_chr7_0.1G1118",
  CO         = "V3.Lp_chr7_0.1G2870",
  CO9        = "V3.Lp_chr1_0G11066",
  CO2        = "V3.Lp_chr6_0G15034"
)

days  <- levels(coldata$day)     # Day_1, Day_10, Day_30
times <- levels(coldata$time)    # 1, 5, 17

all_tp_results <- list()

for (d in days) {
  for (t in times) {
    
    message("Processing: ", d, "  time ", t)
    
    # Subset samples for this day/time
    subset_samples <- rownames(coldata[coldata$day == d &
                                         coldata$time == t, ])
    
    counts_sub  <- counts[, subset_samples, drop = FALSE]
    coldata_sub <- coldata[subset_samples, ]
    
    # Must have both C and N
    if (length(unique(coldata_sub$treatment)) < 2) {
      message("Skipping ", d, " ", t, " — missing C or N samples")
      next
    }
    
    # Build DESeq2 object
    dds_sub <- DESeqDataSetFromMatrix(
      countData = counts_sub,
      colData   = coldata_sub,
      design    = ~ treatment
    )
    
    dds_sub <- DESeq(dds_sub)
    
    # Extract N vs C
    res_sub <- results(dds_sub, contrast = c("treatment", "N", "C"))
    res_df  <- as.data.frame(res_sub)
    
    res_df$gene <- rownames(res_df)
    res_df$day  <- d
    res_df$time <- t
    
    key <- paste0(d, "_", t)
    all_tp_results[[key]] <- res_df
  }
}

# Combine all day × time results
all_tp_DE <- bind_rows(all_tp_results)

# Restrict to your gene panel + panel-level FDR
DESeq2_small <- all_tp_DE %>%
  filter(gene %in% gene_map) %>%
  mutate(
    FDR_small = p.adjust(pvalue, method = "BH"),
    logFC     = log2FoldChange,
    bio_sig   = abs(logFC) >= 1
  )

# ----------------------------------------------------------
# N vs C STATS TABLE (replaces former within-C t-tests)
# ----------------------------------------------------------
save_nvc_stats_table <- function(short_name, gene_id, out_dir = "gene_expression") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  stats_tbl <- DESeq2_small %>%
    filter(gene == gene_id) %>%
    mutate(
      stat_sig = FDR_small < 0.05,
      bio_sig  = abs(logFC) >= 1,
      sig      = stat_sig & bio_sig
    ) %>%
    select(
      day, time,
      log2FoldChange, lfcSE, stat,
      pvalue, padj,
      FDR_small, sig
    ) %>%
    arrange(
      match(day,  c("Day_1", "Day_10", "Day_30")),
      match(time, c("1", "5", "17"))
    )
  
  out_path <- file.path(out_dir, paste0(short_name, "_stats_table.csv"))
  write.csv(stats_tbl, out_path, row.names = FALSE)
  
  cat("\n================ N vs C STATS TABLE ================\n")
  print(stats_tbl)
  
  return(stats_tbl)
}


plot_gene <- function(short_name,
                      out_dir = "gene_expression",
                      save_png = TRUE,
                      save_svg = TRUE,
                      width = 6.5,
                      height = 8,
                      dpi = 300) {
  
  gene_id <- gene_map[[short_name]]
  if (is.null(gene_id)) stop("Gene not found in gene_map")
  if (!gene_id %in% rownames(log_cpm)) stop("Gene ID not found in rownames(log_cpm): ", gene_id)
  
  # Make sure output directory exists
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  expr <- log_cpm[gene_id, ]
  
  df <- data.frame(
    sample = names(expr),
    expression = as.numeric(expr)
  ) %>%
    mutate(
      day       = coldata[sample, "day"],
      time      = coldata[sample, "time"],
      treatment = coldata[sample, "treatment"],
      replicate = coldata[sample, "replicate"],
      day_time  = paste0(day, "_", time),
      rep_label = paste0(treatment, replicate)
    )
  
  # ----------------------------------------------------------
  # SAVE N vs C STATS TABLE (DESeq2 only)
  # ----------------------------------------------------------
  nvc_stats <- save_nvc_stats_table(short_name, gene_id, out_dir = out_dir)
  # ----------------------------------------------------------
  
  mean_df <- df %>%
    group_by(day, time, treatment) %>%
    summarise(mean_expr = mean(expression), .groups = "drop") %>%
    mutate(
      time = factor(time, levels = c("1", "5", "17")),
      day  = factor(day,  levels = c("Day_1", "Day_10", "Day_30"))
    )
  
  df_wide <- df %>%
    select(day_time, rep_label, expression) %>%
    pivot_wider(names_from = rep_label, values_from = expression) %>%
    mutate(
      C_mean = rowMeans(select(., starts_with("C")), na.rm = TRUE),
      N_mean = rowMeans(select(., starts_with("N")), na.rm = TRUE)
    )
  
  # --- DE significance lookup (small-set FDR) ---
  ymax_gene <- max(df$expression, na.rm = TRUE) * 1.1
  
  de_info <- DESeq2_small %>%
    filter(gene == gene_id) %>%
    mutate(
      stat_sig = FDR_small < 0.05,
      bio_sig  = abs(logFC) >= 1,
      sig      = stat_sig & bio_sig,
      ypos     = ymax_gene
    )
  
  final_table <- df_wide %>%
    separate(day_time, into = c("prefix", "day_num", "time"), sep = "_", remove = FALSE) %>%
    mutate(
      day  = paste0("Day_", day_num),
      time = time
    ) %>%
    mutate(
      day  = factor(day,  levels = c("Day_1", "Day_10", "Day_30")),
      time = factor(time, levels = c("1", "5", "17"))
    ) %>%
    left_join(
      de_info %>% select(day, time, logFC),
      by = c("day", "time")
    ) %>%
    arrange(
      match(day,  c("Day_1", "Day_10", "Day_30")),
      match(time, c("1", "5", "17"))
    ) %>%
    select(
      day_time,
      logFC,
      `C (mean)` = C_mean,
      `N (mean)` = N_mean,
      C1, C2, C3, C4,
      N1, N2, N3, N4
    )
  
  # ----------------------------------------------------------
  # Labels (title removed as requested)
  # ----------------------------------------------------------
  ylab_expr <- parse(
    text = sprintf("italic('%s')~' expression ('*log[2]*'CPM)'", short_name)
  )
  
  # ----------------------------------------------------------
  # Facet labels: Day_1 -> Day 1 (etc.)
  # ----------------------------------------------------------
  day_labeller <- labeller(day = function(x) gsub("Day_", "Day ", x))
  
  # ----------------------------------------------------------
  # Glue panels: x-axis only at bottom when supported
  # ----------------------------------------------------------
  facet_layer <- NULL
  
  # ggplot2 changed these facet args across versions.
  # Use only args that exist for the detected version.
  if (utils::packageVersion("ggplot2") >= "3.4.0") {
    facet_layer <- facet_wrap(
      ~ day, ncol = 1,
      scales = "fixed",
      labeller = day_labeller,
      strip.position = "top",
      axes = "margins"       # OK in >= 3.4.0
      # axis_labels removed for compatibility
    )
  } else {
    # Fallback that won't error on older ggplot2
    facet_layer <- facet_grid(
      day ~ .,
      labeller = day_labeller
    )
  }
  
  # ----------------------------------------------------------
  # In-panel day label data (bottom-right, behind points)
  # ----------------------------------------------------------
  label_df <- mean_df %>%
    distinct(day) %>%
    mutate(
      label = gsub("Day_", "Day ", as.character(day)),
      x = Inf,   # right edge
      y = -Inf   # bottom edge (will be nudged into the panel)
    )
  
  # ----------------------------------------------------------
  # Plot (mean triangles removed; mean lines kept)
  # ----------------------------------------------------------
  base_size <- 11 * 1  # all text ~1.5× bigger
  
  p <- ggplot() +
    # In-panel day labels (drawn first so they sit behind data)
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,   # nudge inside from right
      vjust = -0.6,   # nudge up from bottom
      color = "black",
      alpha = 0.15,   # semi-transparent
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
      aes(x = time, y = Inf, label = ifelse(sig, "*", "")),
      vjust = 1,      # adjust 1.1–1.6 depending on your text size
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
      # title removed
    ) +
    scale_x_discrete(labels = c("1" = "NB+1h", "5" = "NB+5h", "17" = "NB+17h")) +
    coord_cartesian(ylim = c(0, NA)) +
    theme(
      # Legend: bottom + bigger (keys)
      legend.position = "none",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(1.2, "lines"),
      
      # Keep clear separators between panels
      panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.8),
      panel.spacing.y = grid::unit(0.15, "lines"),
      
      # Hide facet strip labels since we draw labels inside panels
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    guides(color = guide_legend(title = "Treatment", nrow = 1))
  
  print(p)
  
  # -------------------
  # SAVE PLOTS
  # -------------------
  if (save_png) {
    png_path <- file.path(out_dir, paste0(short_name, "_expression_plot.png"))
    ggsave(
      filename = png_path,
      plot = p,
      width = 2.85,
      height = 3.14,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  if (save_svg) {
    svg_path <- file.path(out_dir, paste0(short_name, "_expression_plot.svg"))
    ggsave(
      filename = svg_path,
      plot = p,
      width = 3.145,
      height = 2.87,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  return(final_table)
}

plot_gene("VRN2b")
plot_gene("VRN2a")
plot_gene("VRN1")
plot_gene("FT3")
plot_gene("FTL9")
plot_gene("ELF3")
plot_gene("ELF4")
plot_gene("GI")
plot_gene("LUX")
plot_gene("LNK1_like")
plot_gene("COL13_like")
plot_gene("PHYC")
plot_gene("PHYB")
plot_gene("PPD1")
plot_gene("CO")
plot_gene("CO9")
plot_gene("CO2")

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
