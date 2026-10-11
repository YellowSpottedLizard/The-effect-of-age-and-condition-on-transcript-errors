# ============================================================================
# Figure S4 — Mismatch rate by Age x genotype (single panel)
#  
# ============================================================================



library(here)
library(tidyverse)
library(ggbeeswarm)

## ---- paths ---------------------------------------------------------------- #
RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here("Fig S4_genic_unfiltered")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILTER <- "unfiltered"



#   unfiltered, nuclear CDS positions only:
#     numerator   = cds_mutations  (all retained mismatches at CDS positions)
#     denominator = depth_sum_cds
NUM <- "cds_mutations"
DEN <- "depth_sum_cds"

## ---- metadata ------------------------------------------------------------- #
meta <- read_tsv(
  META_FILE, skip = 1,
  col_names = c("row_index", "sample_id", "Age", "Mutation_Status", "Deletion"),
  show_col_types = FALSE
) %>%
  mutate(
    Age             = factor(Age, c("Y", "O"), c("Young", "Old")),
    Mutation_Status = factor(Mutation_Status, c("wt", "mu"), c("Wild type", "Deletion")),
    Deletion        = factor(Deletion)
  ) %>%
  select(-row_index)

## ---- loader --------------------------------------------------------------- #
# num / den are CHARACTER VECTORS of column names; multiple entries are summed
# row-wise.
load_rate <- function(filter, num, den) {
  f <- file.path(RESULTS_DIR, paste0(filter, ".summary.tsv"))
  if (!file.exists(f)) {
    avail <- sub("\\.summary\\.tsv$", "",
                 list.files(RESULTS_DIR, pattern = "\\.summary\\.tsv$"))
    stop("No summary.tsv for filter '", filter, "'.\nAvailable:\n  ",
         paste(avail, collapse = ", "), call. = FALSE)
  }
  d <- read_tsv(f, show_col_types = FALSE)
  
  miss <- setdiff(c(num, den), names(d))
  if (length(miss)) {
    stop("Filter '", filter, "' has no column(s): ", paste(miss, collapse = ", "),
         "\nAvailable columns:\n  ", paste(names(d), collapse = ", "), call. = FALSE)
  }
  
  d %>%
    left_join(meta, by = "sample_id") %>%
    filter(!is.na(Age)) %>%
    mutate(
      numerator   = rowSums(across(all_of(num)), na.rm = TRUE),
      denominator = rowSums(across(all_of(den)), na.rm = TRUE),
      rate_per_Mb = numerator / denominator * 1e6
    )
}

dat <- load_rate(FILTER, NUM, DEN)

## ---- y window ------------------------------------------------------------- #
# Single panel, so the window is just the data range plus a little headroom.
panel_range <- function(df) {
  m <- df %>%
    group_by(Age, Mutation_Status) %>%
    summarise(mean = mean(rate_per_Mb),
              sem  = if (n() > 1) sd(rate_per_Mb) / sqrt(n()) else 0,
              .groups = "drop")
  range(c(df$rate_per_Mb, m$mean - m$sem, m$mean + m$sem), na.rm = TRUE)
}

rng  <- panel_range(dat)
SPAN <- diff(rng) * 1.12
mid  <- mean(rng)
lo   <- mid - SPAN / 2
hi   <- mid + SPAN / 2
if (lo < 0) { hi <- hi - lo; lo <- 0 }               # slide up, never clip below zero
lims <- c(lo, hi)

TICK <- signif(SPAN / 4, 1)

## ---- panel drawing -------------------------------------------------------- #
# x = genotype, colour/dodge/lines = age: one line per age group across genotypes
pal <- c("Young" = "#33BB36", "Old" = "#9B77BB")
W   <- 0.6

panel <- function(df, title, ylim, show_ylab = TRUE) {
  means <- df %>%
    group_by(Age, Mutation_Status) %>%
    summarise(mean = mean(rate_per_Mb),
              sem  = if (n() > 1) sd(rate_per_Mb) / sqrt(n()) else 0,
              .groups = "drop")
  
  ggplot(df, aes(Mutation_Status, rate_per_Mb, fill = Age, group = Age)) +
    geom_beeswarm(dodge.width = W, cex = 1.3, size = 2.5, alpha = 1,
                  priority = "density", shape = 21, colour = "white", stroke = 0.3) +
    geom_line(data = means, inherit.aes = FALSE,
              aes(Mutation_Status, mean, group = Age),
              position = position_dodge(W), linewidth = 0.7, colour = "black") +
    geom_errorbar(data = means, inherit.aes = FALSE,
                  aes(Mutation_Status, ymin = mean - sem, ymax = mean + sem, group = Age),
                  position = position_dodge(W), width = 0.14,
                  linewidth = 0.5, colour = "black") +
    geom_point(data = means, inherit.aes = FALSE,
               aes(Mutation_Status, mean, group = Age),
               position = position_dodge(W), size = 2.8, colour = "black") +
    scale_fill_manual(values = pal, name = "Age") +
    scale_y_continuous(breaks = scales::breaks_width(TICK)) +
    coord_cartesian(ylim = ylim) +                    # window only — no data dropped
    labs(x = "Genotype", y = if (show_ylab) "Mismatches / Mb" else NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text  = element_text(size = 11))
}


figS4 <- panel(dat, "Nuclear coding positions, no filtering", lims)

figS4

ggsave(file.path(OUT_DIR, "FigS4_mismatch_rate_genic.pdf"), figS4, width = 4.2, height = 4.4)
ggsave(file.path(OUT_DIR, "FigS4_mismatch_rate_genic.png"), figS4, width = 4.2, height = 4.4, dpi = 300)
ggsave(file.path(OUT_DIR, "FigS4_mismatch_rate_genic.svg"), figS4, width = 5.5, height = 5.5, dpi = 300)

## ---- accompanying tables -------------------------------------------------- #
group_summary <- dat %>%
  mutate(panel = paste0("S4_", FILTER, "_genic")) %>%
  group_by(panel, Age, Mutation_Status) %>%
  summarise(mean_per_Mb = mean(rate_per_Mb),
            sd_per_Mb   = sd(rate_per_Mb),
            sem_per_Mb  = sd(rate_per_Mb) / sqrt(n()),
            n = n(), .groups = "drop")
write_csv(group_summary, file.path(OUT_DIR, "FigS4_group_summary.csv"))

persample <- dat %>%
  transmute(panel = paste0("S4_", FILTER, "_genic"),
            sample_id, Age, Mutation_Status, Deletion, rate_per_Mb)
write_csv(persample, file.path(OUT_DIR, "FigS4_persample_rates.csv"))

message("Fig S4 written to: ", OUT_DIR)


