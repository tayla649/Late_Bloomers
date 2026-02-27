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

# =========================
# 1) Libraries
# =========================
# install.packages(c("readr","dplyr","tidyr","ggplot2","emmeans","rcompanion","svglite"))
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(rcompanion)
library(svglite)

# =========================
# 2) Read CSV (values are ± days from pooled 2022/23 mean)
# =========================
dat <- readr::read_csv("aberviris_preliminary_diff_from_meanHD.csv",
                       show_col_types = FALSE)

# Optional sanity check: ensure required cols exist
req <- c("PlantID","2022_HD","2023_HD","2024_HD","2024_EXT","2024_NB")
stopifnot(all(req %in% names(dat)))

# =========================
# 3) Baseline per plant = mean(2022_HD, 2023_HD), ignoring NAs
# =========================
dat <- dat %>%
  dplyr::mutate(
    baseline_22_23 = rowMeans(dplyr::across(c(`2022_HD`, `2023_HD`)), na.rm = TRUE)
  )

# =========================
# 4) Re-reference 2024 to the plant-specific baseline (Control/EXT/NB)
# =========================
dat <- dat %>%
  dplyr::mutate(
    `2024_HD_rel`  = dplyr::if_else(!is.na(`2024_HD`)  & !is.na(baseline_22_23),
                                    `2024_HD`  - baseline_22_23, NA_real_),
    `2024_EXT_rel` = dplyr::if_else(!is.na(`2024_EXT`) & !is.na(baseline_22_23),
                                    `2024_EXT` - baseline_22_23, NA_real_),
    `2024_NB_rel`  = dplyr::if_else(!is.na(`2024_NB`)  & !is.na(baseline_22_23),
                                    `2024_NB`  - baseline_22_23, NA_real_)
  )

# =========================
# 5) Reshape for plotting/analysis (ONLY the *_rel columns)
#    (use fully qualified tidyr::pivot_longer, and dplyr::select)
# =========================
plot_dat <- dat %>%
  dplyr::select(PlantID, `2024_HD_rel`, `2024_EXT_rel`, `2024_NB_rel`) %>%
  tidyr::pivot_longer(
    cols = c(`2024_HD_rel`, `2024_EXT_rel`, `2024_NB_rel`),
    names_to = "Treatment", values_to = "DeltaRel"
  ) %>%
  dplyr::filter(!is.na(DeltaRel)) %>%
  dplyr::mutate(
    Treatment = factor(
      Treatment,
      levels = c("2024_HD_rel", "2024_EXT_rel", "2024_NB_rel"),
      labels = c("Control", "Extended Daylength", "Night Break")
    )
  )

# =========================
# 6) Summary stats per treatment (IQR, median, mean±SE, n)
# =========================
summary_iqr <- plot_dat %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    mean   = mean(DeltaRel),
    se     = sd(DeltaRel)/sqrt(n),
    median = median(DeltaRel),
    q1     = quantile(DeltaRel, 0.25, names = FALSE),
    q3     = quantile(DeltaRel, 0.75, names = FALSE),
    iqr    = IQR(DeltaRel),
    .groups = "drop"
  )
print(summary_iqr)

# =========================
# 7) ANOVA + Tukey HSD and compact-letter display (no multcompView)
# =========================
fit <- lm(DeltaRel ~ Treatment, data = plot_dat)

# Tukey pairwise via emmeans
emm <- emmeans::emmeans(fit, ~ Treatment)
tk  <- summary(pairs(emm, adjust = "tukey"))

# Convert pair labels to "A-B" for cldList
tukey_tbl <- tk %>%
  as.data.frame() %>%
  dplyr::transmute(
    Comparison = gsub(" - ", "-", contrast),
    P.adj      = p.value
  )

# Letters (shared letter = not significantly different at alpha = 0.05)
letters_tbl <- rcompanion::cldList(P.adj ~ Comparison, data = tukey_tbl, threshold = 0.05) %>%
  dplyr::rename(Treatment = Group, .group = Letter) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = levels(plot_dat$Treatment)))

# y-positions above boxes for the letters
y_pos <- plot_dat %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(y = max(DeltaRel, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(y = y + 0.06 * diff(range(plot_dat$DeltaRel, na.rm = TRUE)))

letters_df <- dplyr::left_join(y_pos, letters_tbl, by = "Treatment")

# =========================
# 8) Plain boxplot + jitter + Tukey letters (A4-third layout)
# =========================
set.seed(5)
rng  <- range(plot_dat$DeltaRel); lim <- max(abs(rng))*1.1
n_df <- plot_dat %>% dplyr::count(Treatment, name = "n")

p_box <- ggplot(plot_dat, aes(Treatment, DeltaRel)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "white", color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = seq(-100, 100, by = 10), color = "grey90", linewidth = 0.25) +
  geom_boxplot(width = 0.45, color = "grey20", fill = "#DCE9F7",
               outlier.shape = 16, outlier.size = 1.6, outlier.alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.12, height = 0),
             size = 1.4, alpha = 0.6, color = "#2C7FB8") +
  geom_text(data = n_df, aes(label = paste0("n=", n), y = lim*1.02),
            vjust = -0.4, size = 3.2, color = "grey20") +
  # Tukey letters
  geom_text(data = letters_df, aes(Treatment, y, label = .group),
            vjust = 0, size = 3.4, fontface = "bold") +
  coord_cartesian(ylim = c(-lim, lim*1.10), clip = "off") +
  labs(
    title = "2024 heading vs each plant’s 2022/23 baseline",
    subtitle = "Box = IQR with median; letters = Tukey HSD groups (α = 0.05); 0 = plant baseline",
    x = NULL,
    y = "Δ days from baseline"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title    = element_text(face = "bold", size = 10.5, margin = margin(b = 2)),
    plot.subtitle = element_text(color = "grey30", size = 8.5, margin = margin(b = 6)),
    axis.title.y  = element_text(size = 9, margin = margin(r = 6)),
    axis.text.x   = element_text(size = 8.5, margin = margin(t = 2)),
    axis.text.y   = element_text(size = 8.5),
    plot.margin   = margin(2, 2, 2, 2)
  )

print(p_box)

# =========================
# 9) Save A4‑third SVG
# =========================
ggsave("aberviris_2024_vs_baseline_boxplot_tukey_A4third.svg", p_box,
       width = 16.0, height = 8.2, units = "cm", dpi = 400)
