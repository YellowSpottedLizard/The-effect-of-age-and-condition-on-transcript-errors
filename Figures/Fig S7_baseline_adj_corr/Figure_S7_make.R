# ============================================================================
# Figure S7 — Baseline-adjusted correlation of Age vs Status spectra
#   x : Age effect    (Old - Young),        residual from baseline
#   y : Status effect (Deletion - Wild type), residual from baseline
#
# Per-channel additive effects (errors/Mb) come from the same NB-GLMM as
# Figs S2/S3 (depth offset, Deletion random effect, emmeans regrid contrasts).
# Because both effects scale with each channel's baseline rate, each is
# regressed on baseline through the origin and the residuals are correlated.
# Pearson r with two-sided permutation p, plus Spearman rho.
#
# ============================================================================
rm(list = ls())

library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(ggrepel)

## ---- CHOOSE THE FILTER ---------------------------------------------------- #
FILTER <- "across_1pct"          # or "across_1pct"

RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here("Fig S7_baseline_adj_corr")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SPEC_RE <- "^spec_[ACGU]>[ACGU]$"
N_PERM  <- 1e4

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
  filter(!is.na(Age)) %>%
  mutate(depth = depth_sum_genespan)

SPEC_COLS     <- grep(SPEC_RE, names(DATAFRAME), value = TRUE)
channel_names <- sub("^spec_", "", SPEC_COLS)
stopifnot(length(channel_names) == 12)

## ---- per-channel additive models ------------------------------------------ #
add_diff <- function(m, spec) {
  emmeans(m, spec, offset = log(1e6)) %>% regrid() %>%
    contrast("revpairwise") %>% summary() %>% as.data.frame() %>% pull(estimate) %>% .[1]
}

channel_summary <- map_dfr(channel_names, function(ch) {
  d <- DATAFRAME %>% mutate(y = .data[[paste0("spec_", ch)]])
  m <- glmmTMB(y ~ Age + Mutation_Status + offset(log(depth)) + (1 | Deletion),
               family = nbinom2, data = d)
  tibble(channel     = ch,
         rate_per_Mb = mean(d$y / d$depth * 1e6),   # baseline: mean per-sample rate
         age_diff    = add_diff(m, ~ Age),
         status_diff = add_diff(m, ~ Mutation_Status))
})

## ---- baseline adjustment -------------------------------------------------- #
diag <- channel_summary %>%
  transmute(channel, baseline = rate_per_Mb, age_d = age_diff, status_d = status_diff) %>%
  mutate(age_res    = resid(lm(age_d    ~ 0 + baseline)),
         status_res = resid(lm(status_d ~ 0 + baseline)))

cat(sprintf("\nr(age, baseline) = %.2f | r(status, baseline) = %.2f | r(age, status) = %.2f\n",
            cor(diag$age_d, diag$baseline), cor(diag$status_d, diag$baseline),
            cor(diag$age_d, diag$status_d)))

## ---- correlation statistics ----------------------------------------------- #
# two-sided permutation p for a correlation
perm_p <- function(x, y, method = "pearson", n = N_PERM) {
  obs <- cor(x, y, method = method)
  (1 + sum(abs(replicate(n, cor(x, sample(y), method = method))) >= abs(obs))) / (n + 1)
}

set.seed(1)
r   <- cor(diag$age_res, diag$status_res)
rho <- cor(diag$age_res, diag$status_res, method = "spearman")
pp  <- perm_p(diag$age_res, diag$status_res)
p_txt <- if (pp < 0.001) "perm. p < 0.001" else sprintf("perm. p = %.3f", pp)
lab   <- sprintf("r = %.2f (%s)\n  \u03C1 = %.2f", r, p_txt, rho)
cat(sprintf("Baseline-adjusted additive effects: %s\n", gsub("\n", " | ", lab)))

## ---- figure --------------------------------------------------------------- #
figS9 <- ggplot(diag, aes(age_res, status_res)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#33BB36",
              fill = "#fdd0a2", linewidth = 0.7, alpha = 0.4) +
  geom_point(size = 2.8, colour = "#9B77BB") +
  geom_text_repel(aes(label = channel), size = 3, seed = 1,
                  min.segment.length = 0, box.padding = 0.4, colour = "grey25") +
  annotate("text", x = -Inf, y = Inf, label = lab,
           hjust = -0.08, vjust = 1.3, size = 3.6, colour = "grey20", lineheight = 0.95) +
  labs(x = "Age effect, residual from baseline (mismatch / Mb)",
       y = "Status effect, residual from baseline (mismatch / Mb)",
       title = "Baseline-adjusted additive effects") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

figS9

base <- paste0("FigS9_baseline_adj_corr_", FILTER)
ggsave(file.path(OUT_DIR, paste0(base, ".pdf")), figS9, width = 6.27, height = 5.2)
ggsave(file.path(OUT_DIR, paste0(base, ".png")), figS9, width = 6.27, height = 5.2, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(base, ".svg")), figS9, width = 6.27, height = 5.2, dpi = 900)

## ---- accompanying tables -------------------------------------------------- #
channel_order <- c("C>A","C>G","C>U","U>A","U>C","U>G",
                   "A>C","A>G","A>U","G>A","G>C","G>U")
write_csv(diag %>% mutate(channel = factor(channel, channel_order)) %>% arrange(channel),
          file.path(OUT_DIR, paste0(base, "_channel_table.csv")))
write_csv(tibble(filter = FILTER, r = r, rho = rho, perm_p = pp,
                 n_perm = N_PERM, n_channels = nrow(diag)),
          file.path(OUT_DIR, paste0(base, "_stats.csv")))

message("Written to: ", OUT_DIR)


