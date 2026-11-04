# ============================================================================
# Figure 1 — Mismatch rate by Age x genotype (three panels)
#   A : unfiltered          total mismatches / Mb          (total depth)
#   B : across_1pct    syn + missense mismatches / Mb (genic depth)
#   C : unfiltered        synonymous mismatches / Mb     (synonymous depth)
#
# ============================================================================
rm(list = ls())

library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(glmmTMB)
library(emmeans)

## ---- paths ---------------------------------------------------------------- #
RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here("Fig 1_abc_all_freq_syn")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FILTER_A <- "unfiltered"
FILTER_B <- "across_1pct"
FILTER_C <- "unfiltered"

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
# row-wise. This is what lets panel B use synonymous + missense as one numerator.
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
      numerator   = rowSums(across(all_of(num)),  na.rm = TRUE),
      denominator = rowSums(across(all_of(den)), na.rm = TRUE),
      rate_per_Mb = numerator / denominator * 1e6
    )
}

datA <- load_rate(FILTER_A, "total_mutations",      "depth_sum_total")
datB <- load_rate(FILTER_B, "cds_mutations",        "depth_sum_cds")
datC <- load_rate(FILTER_C, "synonymous_mutations", "depth_sum_synonymous")

## ---- model-based group means + 95% CI ------------------------------------- #
# Same model family as the text: NB GLMM with a Deletion random effect and a
# log-depth offset. emmeans(offset = log(1e6)) + regrid() puts estimates and
# CIs on the additive per-Mb scale, matching the contrasts reported in the text.
# If the text model has extra terms (dispersion formula, random slopes, ...),
# copy that formula here so the bars match the reported CIs exactly.
group_ci <- function(df) {
  m <- glmmTMB(
    numerator ~ Age * Mutation_Status + (1 | Deletion) + offset(log(denominator)),
    family = nbinom2, data = df
  )
  emmeans(m, ~ Age * Mutation_Status, offset = log(1e6)) %>%
    regrid() %>%
    confint(level = 0.95) %>%
    as_tibble() %>%
    rename(mean = any_of(c("response", "rate", "emmean")),
           lo   = any_of(c("lower.CL", "asymp.LCL")),
           hi   = any_of(c("upper.CL", "asymp.UCL"))) %>%
    select(Age, Mutation_Status, mean, lo, hi)
}

ciA <- group_ci(datA)
ciB <- group_ci(datB)
ciC <- group_ci(datC)

## ---- equal-span, panel-specific y windows --------------------------------- #
# Each panel gets its own window centred on its own data, but every window is
# the SAME width in data units — so vertical distances are comparable across
# panels even though the absolute values differ.

panel_range <- function(df, ci) {
  range(c(df$rate_per_Mb, ci$lo, ci$hi), na.rm = TRUE)
}

rngs <- list(A = panel_range(datA, ciA),
             B = panel_range(datB, ciB),
             C = panel_range(datC, ciC))

# A and B share one window covering both ranges; C keeps its own position
rng_AB <- range(c(rngs$A, rngs$B))
SPAN   <- max(diff(rng_AB), diff(rngs$C))          # same span in all three panels

make_lims <- function(r, span = SPAN) {
  mid <- mean(r)
  lo  <- mid - span / 2
  hi  <- mid + span / 2
  if (lo < 0) { hi <- hi - lo; lo <- 0 }         # slide up, never clip below zero
  c(lo, hi)
}

lims <- list(A = make_lims(rng_AB), B = make_lims(rng_AB), C = make_lims(rngs$C))

# identical tick spacing in every panel (same visual grid, different offsets)
TICK <- signif(SPAN / 4, 1)

## ---- panel drawing -------------------------------------------------------- #
# x = genotype, colour/dodge/lines = age: one line per age group across genotypes
pal <- c("Young" = "#33BB36", "Old" = "#9B77BB")
W   <- 0.6

panel <- function(df, ci, title, ylim, show_ylab = FALSE) {
  ggplot(df, aes(Mutation_Status, rate_per_Mb, fill = Age, group = Age)) +
    geom_beeswarm(dodge.width = W, cex = 1.3, size = 1.5, alpha = 1,
                  priority = "density", shape = 21, colour = "white", stroke = 0.3) +
    geom_line(data = ci, inherit.aes = FALSE,
              aes(Mutation_Status, mean, group = Age),
              position = position_dodge(W), linewidth = 0.3, colour = "black") +
    geom_errorbar(data = ci, inherit.aes = FALSE,
                  aes(Mutation_Status, ymin = lo, ymax = hi, group = Age),
                  position = position_dodge(W), width = 0.14,
                  linewidth = 0.3, colour = "black") +
    geom_point(data = ci, inherit.aes = FALSE,
               aes(Mutation_Status, mean, group = Age),
               position = position_dodge(W), size = 0.6, colour = "black") +
    scale_fill_manual(values = pal, name = "Age") +
    scale_y_continuous(breaks = scales::breaks_width(TICK)) +
    coord_cartesian(ylim = ylim) +                    # window only — no data dropped
    labs(x = "Genotype", y = if (show_ylab) "Mismatches / Mb" else NULL, title = title) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_text(size = 8),
          axis.text  = element_text(size = 7))
}

pA <- panel(datA, ciA, "All callable positions",       lims$A, show_ylab = TRUE)
pB <- panel(datB, ciB, "1% filter (coding positions)", lims$B)
pC <- panel(datC, ciC, "Synonymous mismatches",        lims$C)

## ---- assemble ------------------------------------------------------------- #
pA <- pA + labs(tag = "A") + theme(plot.tag.position = c(0, 1),
                                   plot.tag = element_text(hjust = -2.6, vjust = 0))
pB <- pB + labs(tag = "B") + theme(plot.tag.position = c(0, 1),
                                   plot.tag = element_text(hjust = -1.4, vjust = 0))
pC <- pC + labs(tag = "C") + theme(plot.tag.position = c(0, 1),
                                   plot.tag = element_text(hjust = -1.7, vjust = 0))
fig1 <- (pA | pB | pC) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box.spacing = unit(1.5, "pt"),
        legend.text = element_text(size = 8))

fig1

ggsave(file.path(OUT_DIR, "Fig1_mismatch_rate_panels_u.pdf"), fig1, width = 6.27, height = 3)
ggsave(file.path(OUT_DIR, "Fig1_mismatch_rate_panels_u.png"), fig1, width = 6.27, height = 3, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig1_mismatch_rate_panels_u.svg"), fig1, width = 6.27, height = 3, dpi = 300)

## ---- accompanying tables -------------------------------------------------- #
# Raw group means/SD/SEM alongside the plotted model-based means and 95% CIs
raw_summary <- bind_rows(
  datA %>% mutate(panel = paste0("A_", FILTER_A)),
  datB %>% mutate(panel = paste0("B_", FILTER_B)),
  datC %>% mutate(panel = paste0("C_", FILTER_C))
) %>%
  group_by(panel, Age, Mutation_Status) %>%
  summarise(raw_mean_per_Mb = mean(rate_per_Mb),
            sd_per_Mb       = sd(rate_per_Mb),
            sem_per_Mb      = sd(rate_per_Mb) / sqrt(n()),
            n = n(), .groups = "drop")

model_summary <- bind_rows(
  ciA %>% mutate(panel = paste0("A_", FILTER_A)),
  ciB %>% mutate(panel = paste0("B_", FILTER_B)),
  ciC %>% mutate(panel = paste0("C_", FILTER_C))
) %>%
  rename(emm_per_Mb = mean, ci95_lo = lo, ci95_hi = hi)

group_summary <- left_join(raw_summary, model_summary,
                           by = c("panel", "Age", "Mutation_Status"))
write_csv(group_summary, file.path(OUT_DIR, "Fig1_group_summary.csv"))

persample <- bind_rows(
  datA %>% transmute(panel = paste0("A_", FILTER_A), sample_id, Age, Mutation_Status, Deletion, rate_per_Mb),
  datB %>% transmute(panel = paste0("B_", FILTER_B), sample_id, Age, Mutation_Status, Deletion, rate_per_Mb),
  datC %>% transmute(panel = paste0("C_", FILTER_C), sample_id, Age, Mutation_Status, Deletion, rate_per_Mb)
)
write_csv(persample, file.path(OUT_DIR, "Fig1_persample_rates.csv"))

message("Fig 1 written to: ", OUT_DIR)

