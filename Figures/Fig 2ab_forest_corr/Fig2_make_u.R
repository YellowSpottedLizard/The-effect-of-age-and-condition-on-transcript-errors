# ============================================================================
# Figure 2 — Mutation-spectrum changes: Age vs Mutation status
#   A : forest plot of per-channel ADDITIVE effects (errors/Mb), Age & Status
#   B : correlation of the per-channel Age vs Status additive effects
#
# Self-contained: loads + prepares data from stage3_results/, fits the
# per-channel models, builds both panels, and saves the figure + tables into
# this figure's own folder under Manuscript files/Figs/.
#
# NOTE: this fits 12 per-channel negative-binomial GLMMs, so it takes a little
# while to run. Requires an .Rproj or .here marker in the main project folder.
# ============================================================================

rm(list = ls())
library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(ggbeeswarm)
library(ggrepel)
library(patchwork)

## ---- config --------------------------------------------------------------- #
FILTER      <- "across_1pct"                       # which filter this figure uses
RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here("Fig 2ab_forest_corr")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPEC_RE <- "^spec_[ACGU]>[ACGU]$"

## ---- load + prepare data -------------------------------------------------- #
meta <- read_tsv(META_FILE, skip = 1,
                 col_names = c("row_index","sample_id","Age","Mutation_Status","Deletion"),
                 show_col_types = FALSE) %>%
  mutate(Age             = factor(Age, c("Y","O"), c("Young","Old")),
         Mutation_Status = factor(Mutation_Status, c("wt","mu"), c("Control","Deletion")),
         Deletion        = factor(Deletion)) %>%
  select(-row_index)

summ <- read_tsv(file.path(RESULTS_DIR, paste0(FILTER, ".summary.tsv")),
                 show_col_types = FALSE)

DATAFRAME <- summ %>%
  left_join(meta, by = "sample_id") %>%
  filter(!is.na(Age)) %>%
  mutate(depth        = depth_sum_genespan )

SPEC_COLS     <- grep(SPEC_RE, names(DATAFRAME), value = TRUE)
channel_names <- sub("^spec_", "", SPEC_COLS)

## ---- per-channel additive rate models (Option A, depth-offset) ------------ #
# regrid() -> per-Mb scale so contrast() gives an ADDITIVE difference (errors/Mb),
# not a ratio. revpairwise => Old - Young and Deletion - Control (positive = up).
fit_channel <- function(ch) {
  d <- DATAFRAME %>% mutate(y = .data[[paste0("spec_", ch)]])
  glmmTMB(y ~ Age + Mutation_Status + offset(log(depth)) + (1 | Deletion),
          family = nbinom2, data = d)
}
channel_fits <- setNames(map(channel_names, fit_channel), channel_names)

# helper: pull additive contrast (estimate + CI) in a name-robust way
add_contrast <- function(m, spec) {
  d <- emmeans(m, spec, offset = log(1e6)) %>% regrid() %>%
    contrast("revpairwise") %>% summary(infer = TRUE) %>% as.data.frame()
  l <- intersect(c("lower.CL","asymp.LCL"), names(d))[1]
  u <- intersect(c("upper.CL","asymp.UCL"), names(d))[1]
  tibble(estimate = d$estimate[1], lcl = d[[l]][1], ucl = d[[u]][1], p = d$p.value[1])
}

channel_summary <- map_dfr(channel_names, function(ch) {
  m  <- channel_fits[[ch]]
  ea <- add_contrast(m, ~ Age)
  es <- add_contrast(m, ~ Mutation_Status)
  tibble(channel     = ch,
         rate_per_Mb = exp(fixef(m)$cond["(Intercept)"]) * 1e6,
         age_diff = ea$estimate, age_lcl = ea$lcl, age_ucl = ea$ucl, age_p = ea$p,
         status_diff = es$estimate, status_lcl = es$lcl, status_ucl = es$ucl, status_p = es$p)
}) %>%
  mutate(age_p_fdr    = p.adjust(age_p,    "fdr"),
         status_p_fdr = p.adjust(status_p, "fdr"))

# channel ordering by baseline abundance (shared by both panels)
channel_order <- channel_summary %>% arrange(desc(rate_per_Mb)) %>% pull(channel)

## ===========================================================================
## PANEL A — forest plot of additive per-channel effects
## ===========================================================================
## ===========================================================================
## PANEL A — forest plot of additive per-channel effects
## ===========================================================================
forest_dat <- bind_rows(
  channel_summary %>% transmute(channel, diff = age_diff, lcl = age_lcl, ucl = age_ucl,
                                sig = age_p_fdr < 0.05, factor = "Age (Old - Young)"),
  channel_summary %>% transmute(channel, diff = status_diff, lcl = status_lcl, ucl = status_ucl,
                                sig = status_p_fdr < 0.05, factor = "Condition (Deletion - Control)")
) %>%
  mutate(channel = factor(channel, levels = rev(channel_order)),
         # keep both legend keys even if one category is absent
         sig     = factor(sig, levels = c(TRUE, FALSE)))

# common x range across BOTH facets, centred on zero
xr <- max(abs(c(forest_dat$lcl, forest_dat$ucl)), na.rm = TRUE) * c(-0.5, 1) * 1.05

# inside-legend placement (API changed in ggplot2 3.5.0)
legend_inside <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
  theme(legend.position = "inside",
        legend.position.inside = c(0.99, 0.02),
        legend.justification = c(1, 0))
} else {
  theme(legend.position = c(0.99, 0.02),
        legend.justification = c(1, 0))
}

pA <- ggplot(forest_dat, aes(diff, channel)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), orientation = "y", width = 0.28, linewidth = 0.5, colour = "grey40") +  geom_point(aes(colour = sig), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "#33BB36", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns"),
                      drop = FALSE, name = NULL) +
  facet_wrap(~ factor) +                                    # fixed = shared x axis
  coord_cartesian(xlim = xr) +
  scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +  # room at the bottom for the legend
  labs(x = "Mismatch differences / Mb", y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        panel.spacing    = unit(1.1, "lines"),
        legend.direction = "vertical",                       # stacked
        legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
        legend.key.height = unit(0.9, "lines"),
        legend.margin     = margin(2, 4, 2, 4)) +
  legend_inside
pA

## ===========================================================================
## PANEL B — correlation of Age vs Status additive effects (one point/channel)
## ===========================================================================
mech <- channel_summary %>% transmute(channel, age_d = age_diff, status_d = status_diff)

ct_p   <- cor.test(mech$age_d, mech$status_d, method = "pearson")
ct_s   <- cor.test(mech$age_d, mech$status_d, method = "spearman")
set.seed(1)
obs    <- cor(mech$age_d, mech$status_d)
p_perm <- (1 + sum(abs(replicate(1e4, cor(mech$age_d, sample(mech$status_d)))) >= abs(obs))) / (1e4 + 1)
r_lab  <- sprintf("italic(r) == %.2f", unname(ct_p$estimate))

pB <- ggplot(mech, aes(age_d, status_d)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_smooth(method = "lm", se = TRUE, fullrange = TRUE, colour = "#33BB36", 
              fill = "#fdd0a2", linewidth = 0.7, alpha = 0.4) +
  geom_point(size = 2.8, colour = "#9B77BB") +
  geom_text_repel(aes(label = channel), size = 3, seed = 1,
                  min.segment.length = 0, box.padding = 0.4, colour = "grey25") +
  annotate("text", x = 0, y = Inf, label = r_lab, parse = TRUE,
           hjust = -0.3, vjust = 1.6, size = 4, colour = "grey20") +
  labs(x = "Age effect (Old - Young, mismatch / Mb)",
       y = "Condition effect (Deletion - Control, mismatch / Mb)") +
  theme_classic(base_size = 12)

## ---- assemble ------------------------------------------------------------- #
fig2 <- (pA | pB) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")

fig2

ggsave(file.path(OUT_DIR, "Fig2_forest_correlation.pdf"), fig2, width = 11, height = 5.2)
ggsave(file.path(OUT_DIR, "Fig2_forest_correlation.png"), fig2, width = 11, height = 5.2, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_forest_correlation.svg"), fig2, width = 11, height = 5.2,)




## ---- accompanying tables -------------------------------------------------- #
write_csv(channel_summary, file.path(OUT_DIR, "Fig2_channel_effects.csv"))

corr_tbl <- tibble(
  statistic = c("pearson_r", "spearman_rho", "permutation_p", "pearson_p_parametric"),
  value     = c(unname(ct_p$estimate), unname(ct_s$estimate), p_perm, ct_p$p.value)
)
write_csv(corr_tbl, file.path(OUT_DIR, "Fig2_correlation_stats.csv"))

message("Fig 2 written to: ", OUT_DIR)
message(sprintf("Pearson r = %.2f | Spearman rho = %.2f | perm p = %.4f",
                ct_p$estimate, ct_s$estimate, p_perm))


