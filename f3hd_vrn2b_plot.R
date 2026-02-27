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

# ============================================================
# VRN2b_F3: Heading date vs VRN2b allele (A/C), coloured by Family
# Input file: VRN2b_F3.csv
# Columns required: Allele (e.g., "27-Nov"), HD (VRN2b_A/VRN2b_C/na), Family
# ============================================================

# ---- Packages ----
pkgs <- c("tidyverse", "lubridate")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(tidyverse)
library(lubridate)

# ---- Load ----
df <- read_csv("VRN2b_F3.csv", show_col_types = FALSE)

# ---- Clean + Parse ----
dummy_year <- 2025
dec1 <- as.Date(paste0(dummy_year, "-12-01"))

dat <- df %>%
  transmute(
    Family = str_trim(Family),
    # In your file the date string is in column "Allele"
    date_str = str_trim(Allele),
    # In your file the allele label is in column "HD"
    allele = str_trim(HD)
  ) %>%
  mutate(
    allele = na_if(allele, ""),
    allele = if_else(allele == "na", NA_character_, allele),
    allele = recode(allele, "VRN2b_A" = "A", "VRN2b_C" = "C"),
    allele = factor(allele, levels = c("A", "C")),
    parsed = parse_date_time(date_str, orders = c("d-b", "d-B"),
                             locale = "C", quiet = TRUE),
    HD_date = as.Date(update(parsed, year = dummy_year))
  ) %>%
  filter(!is.na(HD_date), !is.na(allele), !is.na(Family))

# ---- Axis breaks (5-day majors; daily minors; label 1-Dec) ----
y_min <- floor_date(min(dat$HD_date), unit = "day")
y_max <- ceiling_date(max(dat$HD_date), unit = "day")

major_breaks <- seq(y_min, y_max, by = "5 days")
minor_breaks <- seq(y_min, y_max, by = "1 day")
major_breaks <- sort(unique(c(major_breaks, dec1)))  # ensure 1-Dec appears

# ---- Plot ----
set.seed(123)

p <- ggplot(dat, aes(x = allele, y = HD_date)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.45,
    colour = "grey20",
    fill = "grey70"
  ) +
  geom_jitter(
    aes(colour = Family),
    width = 0.18, height = 0,
    size = 2, alpha = 0.85
  ) +
  geom_hline(yintercept = dec1, colour = "red3", linewidth = 0.9) +
  scale_y_date(
    breaks = major_breaks,
    minor_breaks = minor_breaks,
    date_labels = "%d-%b",
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Heading Time by VRN2b Allele (F3)",
    subtitle = "Two columns (A vs C). Points coloured by Family. Red line = 1-Dec.",
    x = "VRN2b allele",
    y = "2025 Heading Date",
    colour = "Family"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey30", linewidth = 0.7),
    panel.grid.minor.y = element_line(colour = "grey80", linewidth = 0.3),
    axis.text.x = element_text(hjust = 0.5),
    plot.margin = margin(10, 15, 20, 15)
  )

print(p)

# ---- Save ----
ggsave(
  filename = "VRN2b_F3_heading_A_vs_C_coloured_by_family.png",
  plot = p,
  width = 7.5, height = 5.2, dpi = 300
)
