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

library(svglite)
library(tidyverse)
library(lubridate)

df <- read.csv("nb_hd.csv")

df_long <- df %>%
  pivot_longer(cols = c(n_dates, c_dates),
               names_to = "treatment",
               values_to = "date_str") %>%
  mutate(
    treatment = recode(treatment,
                       n_dates = "Night-Break",
                       c_dates = "Control"),
    # parse with a base year
    date0 = dmy(paste(date_str, "2025")),
    # if it's January, treat as next year (common for seasonal trials)
    date  = if_else(month(date0) == 1, date0 + years(1), date0)
  )

ggplot(df_long, aes(x = treatment, y = date)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.25) +
  geom_jitter(aes(colour = Replicate),
              width = 0.12, height = 0,
              size = 3, alpha = 0.9) +
  scale_y_date(date_labels = "%d-%b", date_breaks = "1 week") +
  labs(x = NULL, y = "Heading date", colour = "Replicate") +
  theme_minimal(base_size = 12)


# A4 with 2.5 cm margins → printable width ≈ 16.0 cm; use ~1/3 page height ≈ 8.2 cm
ggsave("nb_hd_boxplot.svg", p,
       width = 16.0, height = 8.2, units = "cm", dpi = 400)
