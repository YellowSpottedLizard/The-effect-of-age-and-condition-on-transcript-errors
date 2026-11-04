# ============================================================================
# Figures S5 / S6 — pN/pS by Age x Mutation status (two panels)
#   A : weighted dN/dS   (neutral line at 1 shown; 1 is a valid null here)
#   B : unweighted pN/pS (NO line at 1 — spectrum bias shifts the neutral point)
#
# One parameterized script for both supplements:
#   FILTER <- "f1_none"          -> Fig S2  (folder "Fig S2_pNpS_all")
#   FILTER <- "f8_a1pct_c001"    -> Fig S3  (folder "Fig S3_pNpS_filter")
# Set FILTER below and run; it writes into the matching folder automatically.
#
# Self-contained: loads + prepares its own data, builds the figure, and saves
# the figure + a group summary table. Requires an .Rproj/.here marker in the
# main project folder so here() anchors to the project root.
# ============================================================================
rm(list = ls())


library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)

## ---- CHOOSE THE FILTER ---------------------------------------------------- #
FILTER <- "unfiltered"          # unfiltered <- switch to "across_1pct" for Fig S3

# map filter -> output folder (matches the existing Figs/ subfolders)
FOLDER <- dplyr::case_when(
  FILTER == "unfiltered"       ~ "Fig S5_pNpS_all",
  FILTER == "across_1pct" ~ "Fig S6_pNpS_filter",
  TRUE                      ~ paste0("Fig S_pNpS_", FILTER)   # fallback
)

RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here( FOLDER)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- load + prepare data -------------------------------------------------- #
meta <- read_tsv(META_FILE, skip = 1,
                 col_names = c("row_index","sample_id","Age","Mutation_Status","Deletion"),
                 show_col_types = FALSE) %>%
  mutate(Age             = factor(Age, c("Y","O"), c("Young","Old")),
         Mutation_Status = factor(Mutation_Status, c("wt","mu"), c("Wild type","Deletion")),
         Deletion        = factor(Deletion)) %>%
  select(-row_index)

f <- file.path(RESULTS_DIR, paste0(FILTER, ".summary.tsv"))
if (!file.exists(f)) {
  avail <- sub("\\.summary\\.tsv$", "", list.files(RESULTS_DIR, "\\.summary\\.tsv$"))
  stop("No summary.tsv for '", FILTER, "'.\nAvailable:\n  ", paste(avail, collapse = ", "),
       call. = FALSE)
}

DATAFRAME <- read_tsv(f, show_col_types = FALSE) %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(Age))

## ---- panel drawing: boxplots, zoomed to +/- 0.3 around the data ---------- #
pal <- c("Wild type" = "#2c7fb8", "Deletion" = "#d95f0e")

YMAX <- DATAFRAME %>%
  select(dNdS, dNdS_weighted) %>%
  unlist() %>%
  .[is.finite(.)] %>%
  max()
YMAX <- max(YMAX, 1.03) * 1.05

YMAX <- 1.1


panel <- function(col, ylab, title) {
  df <- DATAFRAME %>% mutate(val = .data[[col]]) %>% filter(is.finite(val), val > 0)
  
  ggplot(df, aes(Age, val, fill = Mutation_Status)) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey40") +   # drawn behind the boxes
    geom_boxplot(width = 0.55, position = position_dodge(0.62),
                 outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.5,
                 colour = "grey25", linewidth = 0.5) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(x = "Age", y = ylab, title = title) +
    coord_cartesian(ylim = c(0, YMAX), expand = FALSE) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11, face = "bold"))
}

pA <- panel("dNdS",          "Unweighted mN/mS", "Unweighted")
pB <- panel("dNdS_weighted", "Weighted mN/mS",   "Weighted")



## ---- assemble (shared legend, free y per panel, A/B tags) ----------------- #
figS <- (pA | pB) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

figS

base <- paste0("Fig_pNpS_", FILTER)
ggsave(file.path(OUT_DIR, paste0(base, ".pdf")), figS, width = 8.5, height = 4.6)
ggsave(file.path(OUT_DIR, paste0(base, ".png")), figS, width = 8.5, height = 4.6, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(base, ".svg")), figS, width = 8.5, height = 4.6, dpi = 300)



## ---- accompanying table: geometric mean + t-CI on log scale --------------- #
# geometric mean matches the ratio nature of dN/dS; p_vs1 only meaningful for weighted.
dnds_ci <- function(x, test_vs1 = FALSE) {
  x  <- x[is.finite(x) & x > 0]
  tt <- t.test(log(x))
  tibble(geom_mean = exp(mean(log(x))),
         lcl = exp(tt$conf.int[1]), ucl = exp(tt$conf.int[2]),
         p_vs1 = if (test_vs1) t.test(log(x))$p.value else NA_real_,
         n = length(x))
}

group_summary <- bind_rows(
  DATAFRAME %>% group_by(Age, Mutation_Status) %>%
    summarise(dnds_ci(dNdS_weighted, TRUE),  .groups = "drop") %>% mutate(measure = "weighted_dNdS"),
  DATAFRAME %>% group_by(Age, Mutation_Status) %>%
    summarise(dnds_ci(dNdS, FALSE), .groups = "drop") %>% mutate(measure = "unweighted_pNpS")
) %>% relocate(measure)
write_csv(group_summary, file.path(OUT_DIR, paste0(base, "_group_summary.csv")))

message("Written to: ", OUT_DIR)



