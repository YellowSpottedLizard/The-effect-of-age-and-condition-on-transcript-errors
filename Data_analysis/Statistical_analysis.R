

rm(list = ls())

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(broom.mixed)
library(patchwork)
library(scales)
library(jsonlite)
library(conflicted)
library(performance)
library(conflicted)

library(here)
library(ggbeeswarm)


conflicts_prefer(dplyr::filter)


# -------- CHANGE THE FILTER HERE -------------------------------------------- #
FILTER       <- "unfiltered"   

# """ Filters available in ../stage3 :
#  across_0p01pct, across_0p1pct, across_1pct, across_5pct, unfiltered  """



RESULTS_DIR  <- "../stage3"                              # holds *.summary.tsv etc.
META_FILE    <- "../stage3/metadata_minimal.txt"


available_filters <- sort(sub("\\.summary\\.tsv$", "",
                              list.files(RESULTS_DIR, pattern = "\\.summary\\.tsv$")))

cat("Filters available in", RESULTS_DIR, ":\n  ",
    paste(available_filters, collapse = ", "), "\n")
cat("Analysing:", FILTER, "\n")


summary_path  <- file.path(RESULTS_DIR, paste0(FILTER, ".summary.tsv"))
overview_path <- file.path(RESULTS_DIR, paste0(FILTER, ".overview.tsv"))
params_path   <- file.path(RESULTS_DIR, paste0(FILTER, ".params.json"))
stopifnot("Chosen FILTER has no summary.tsv in RESULTS_DIR" = file.exists(summary_path))



summ     <- read_tsv(summary_path, show_col_types = FALSE)
overview <- if (file.exists(overview_path)) read_tsv(overview_path, show_col_types = FALSE) else NULL
filt_par <- if (file.exists(params_path)) read_json(params_path, show_col_types = FALSE) else NULL

head(overview)





# Metadata carries a leading row-index column whose header is absent, so we read
# positionally and name columns ourselves.
meta <- read_tsv(
  META_FILE, skip = 1,
  col_names = c("row_index", "sample_id", "Age", "Mutation_Status", "Deletion"),
  show_col_types = FALSE
) %>%
  mutate(
    Age             = factor(Age, levels = c("Y", "O"),  labels = c("Young", "Old")),
    Mutation_Status = factor(Mutation_Status, levels = c("wt", "mu"),
                             labels = c("Control", "Deletion")),
    Deletion        = factor(Deletion)
  ) %>%
  select(-row_index)

dat <- summ %>% left_join(meta, by = "sample_id")




unmatched <- dat %>% filter(is.na(Age)) %>% pull(sample_id)
if (length(unmatched))
  warning("Samples in summary with no metadata match: ", paste(unmatched, collapse = ", "))
SPEC_RE <- "^spec_[ACGU]>[ACGU]$"


# tiny helper to print a titled table block to the console
show <- function(df, title = NULL, digits = 3) {
  if (!is.null(title)) cat("\n===", title, "===\n")
  df <- as.data.frame(df)
  num <- sapply(df, is.numeric)
  df[num] <- lapply(df[num], function(x) round(x, digits))
  print(df, row.names = FALSE)
  invisible(df)
}


show(dat %>% count(Mutation_Status, Age) %>%
       pivot_wider(names_from = Age, values_from = n, values_fill = 0),
     "Sample counts per Age x Deletion-status cell")


DATAFRAME <- dat %>%
  mutate(
    depth_reconciles = abs((depth_sum_cds + depth_sum_noncds) - depth_sum_total) <= 1,
    spec_sum         = rowSums(across(matches(SPEC_RE))),
    spec_reconciles  = abs(total_mutations - spec_sum - spectrum_skipped) <= 1
  )






###Choose which one, coding, mito, or all positions
#TOTAL
m_full <- glmmTMB(
  total_mutations ~ Age * Mutation_Status + offset(log(depth_sum_total)) + (1 | Deletion),
  family = nbinom2, data = DATAFRAME
)


###CDS
m_full <- glmmTMB(
  cds_mutations ~ Age * Mutation_Status + offset(log(depth_sum_cds)) + (1 | Deletion),
  family = nbinom2, data = DATAFRAME
)


###MITO
m_full <- glmmTMB(
  total_mutations_mito ~ Age * Mutation_Status + offset(log(depth_sum_total_mito)) + (1 | Deletion),
  family = nbinom2, data = DATAFRAME
)


###SYNONUMOUS
m_full <- glmmTMB(
  synonymous_mutations ~ Age * Mutation_Status + offset(log(depth_sum_synonymous)) + (1 | Deletion),
  family = nbinom2, data = DATAFRAME
)



###SYNONUMOUS MITO
m_full <- glmmTMB(
  synonymous_mutations_mito ~ Age * Mutation_Status + offset(log(depth_sum_synonymous_mito)) + (1 | Deletion),
  family = nbinom2, data = DATAFRAME
)










cat("\n========== INTERACTION TEST (full model) ==========\n")
print(joint_tests(m_full))                 # type-III-style term tests

## ---- 2d. Omnibus term tests, additive model ----------------------------- ##
cat("\n========== MAIN-EFFECT TESTS (additive model) ==========\n")
m_add <- update(m_full, . ~ . - Age:Mutation_Status)
print(joint_tests(m_add))





# --- ADDITIVE differences (per Mb) with 95% CI --------------------------------
# regrid() -> response scale as identity grid, so contrast() subtracts (diff),
# NOT divides (ratio). revpairwise => Old - Young and Deletion - Control (positive = higher).
diff_age  <- emmeans(m_add, ~ Age,             offset = log(1e6)) %>% regrid() %>%
  contrast("revpairwise") %>% confint()
diff_stat <- emmeans(m_add, ~ Mutation_Status, offset = log(1e6)) %>% regrid() %>%
  contrast("revpairwise") %>% confint()

cat("\n========== ADDITIVE DIFFERENCE per Mb: Age (Old - Young), 95% CI ==========\n"); print(diff_age)
cat("\n========== ADDITIVE DIFFERENCE per Mb: Status (Del - Ctrl), 95% CI ==========\n"); print(diff_stat)






## ---- 2e-int. Interaction effect with 95% CI (ADDITIVE only) ------------ ##
int_add <- emmeans(m_full, ~ Age * Mutation_Status, offset = log(1e6)) %>%
  regrid() %>%
  contrast(interaction = c(Age = "revpairwise", Mutation_Status = "revpairwise")) %>%
  summary(infer = TRUE) %>% as.data.frame()

chi2 <- int_add$z.ratio^2                         # 1-df Wald chi-square = z^2

cat("\n===== ADDITIVE interaction per Mb: (Old-Young | Del) - (Old-Young | Ctrl), 95% CI =====\n")
print(int_add)
print(paste("X^2 value" , round(chi2, 2) ) ) 




#####MUTATION RATES -----------------------------------------------------------------------------------------------------

###TOTAL
NUM <- "total_mutations";      DEN <- "depth_sum_total";      TITLE <- "All callable positions"

###CDS
NUM <- "cds_mutations";        DEN <- "depth_sum_cds";        TITLE <- "Coding positions"

###MITO
NUM <- "total_mutations_mito"; DEN <- "depth_sum_total_mito"; TITLE <- "Mitochondrial"


plot_dat <- DATAFRAME %>%
  filter(!is.na(Age)) %>%
  mutate(rate_per_Mb = .data[[NUM]] / .data[[DEN]] * 1e6)

means <- plot_dat %>%
  group_by(Age, Mutation_Status) %>%
  summarise(mean = mean(rate_per_Mb),
            sem  = if (n() > 1) sd(rate_per_Mb) / sqrt(n()) else 0,
            .groups = "drop")


lims    <- range(c(plot_dat$rate_per_Mb, means$mean - means$sem, means$mean + means$sem), na.rm = TRUE)
lims[1] <- max(0, lims[1])
TICK    <- signif(diff(lims) / 4, 1)

pal <- setNames(c("#33BB36", "#9B77BB"), levels(plot_dat$Mutation_Status))
W   <- 0.6


ggplot(plot_dat, aes(Age, rate_per_Mb, fill = Mutation_Status, group = Mutation_Status)) +
  geom_beeswarm(dodge.width = W, cex = 1.3, size = 1.5, alpha = 1,
                priority = "density", shape = 21, colour = "white", stroke = 0.3) +
  geom_line(data = means, inherit.aes = FALSE,
            aes(Age, mean, group = Mutation_Status),
            position = position_dodge(W), linewidth = 0.3, colour = "black") +
  geom_errorbar(data = means, inherit.aes = FALSE,
                aes(Age, ymin = mean - sem, ymax = mean + sem, group = Mutation_Status),
                position = position_dodge(W), width = 0.14,
                linewidth = 0.3, colour = "black") +
  geom_point(data = means, inherit.aes = FALSE,
             aes(Age, mean, group = Mutation_Status),
             position = position_dodge(W), size = 0.6, colour = "black") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(breaks = scales::breaks_width(TICK)) +
  coord_cartesian(ylim = lims) +
  labs(x = "Age", y = "Mismatches / Mb", title = paste0(TITLE, " (", FILTER, ")")) +
  theme_classic(base_size = 8) +
  theme(legend.position    = "bottom",
        legend.box.spacing = unit(1.5, "pt"),
        legend.text        = element_text(family = "Arial", size = 8),
        plot.title         = element_text(size = 8, hjust = 0.5),
        axis.title         = element_text(size = 8),
        axis.text          = element_text(size = 7))







#####MUTATION SPECTRA -----------------------------------------------------------------------------------------------------

library(ggrepel)

spec_dat      <- DATAFRAME %>% filter(!is.na(Age))
SPEC_COLS     <- grep(SPEC_RE, names(spec_dat), value = TRUE)
channel_names <- sub("^spec_", "", SPEC_COLS)



# regrid() -> per-Mb scale, so contrast() gives an ADDITIVE difference (errors/Mb).
# revpairwise => Old - Young and Deletion - Control (positive = up).
channel_fits    <- list()
channel_summary <- tibble()

for (ch in channel_names) {
  d <- spec_dat %>% mutate(y = .data[[paste0("spec_", ch)]])
  m <- glmmTMB(y ~ Age + Mutation_Status + offset(log(depth_sum_genespan)) + (1 | Deletion),
               family = nbinom2, data = d)
  channel_fits[[ch]] <- m
  
  ea <- emmeans(m, ~ Age, offset = log(1e6)) %>% regrid() %>%
    contrast("revpairwise") %>% summary(infer = TRUE) %>% as.data.frame()
  es <- emmeans(m, ~ Mutation_Status, offset = log(1e6)) %>% regrid() %>%
    contrast("revpairwise") %>% summary(infer = TRUE) %>% as.data.frame()
  
  # CI column names differ between asymptotic and finite-df output
  lcl <- intersect(c("lower.CL", "asymp.LCL"), names(ea))[1]
  ucl <- intersect(c("upper.CL", "asymp.UCL"), names(ea))[1]
  
  channel_summary <- bind_rows(channel_summary, tibble(
    channel     = ch,
    rate_per_Mb = exp(fixef(m)$cond["(Intercept)"]) * 1e6,
    age_diff    = ea$estimate[1], age_lcl    = ea[[lcl]][1], age_ucl    = ea[[ucl]][1], age_p    = ea$p.value[1],
    status_diff = es$estimate[1], status_lcl = es[[lcl]][1], status_ucl = es[[ucl]][1], status_p = es$p.value[1]
  ))
}

channel_summary <- channel_summary %>%
  mutate(age_p_fdr    = p.adjust(age_p,    "fdr"),
         status_p_fdr = p.adjust(status_p, "fdr"))

# channel order by baseline abundance (shared by both panels)
channel_order <- channel_summary %>% arrange(desc(rate_per_Mb)) %>% pull(channel)

show(channel_summary, "Per-channel additive effects (per Mb)")


forest_dat <- bind_rows(
  channel_summary %>% transmute(channel, diff = age_diff, lcl = age_lcl, ucl = age_ucl,
                                sig = age_p_fdr < 0.05, factor = "Age (Old - Young)"),
  channel_summary %>% transmute(channel, diff = status_diff, lcl = status_lcl, ucl = status_ucl,
                                sig = status_p_fdr < 0.05, factor = "Condition (Deletion - Control)")
) %>%
  mutate(channel = factor(channel, levels = rev(channel_order)),
         sig     = factor(sig, levels = c(TRUE, FALSE)))   # keep both legend keys

# common x range across both facets, centred on zero
xr <- max(abs(c(forest_dat$lcl, forest_dat$ucl)), na.rm = TRUE) * c(-1, 1) * 1.05

# inside-legend placement (API changed in ggplot2 3.5.0)
if (utils::packageVersion("ggplot2") >= "3.5.0") {
  legend_inside <- theme(legend.position = "inside",
                         legend.position.inside = c(0.99, 0.02),
                         legend.justification = c(1, 0))
} else {
  legend_inside <- theme(legend.position = c(0.99, 0.02),
                         legend.justification = c(1, 0))
}

pA <- ggplot(forest_dat, aes(diff, channel)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), orientation = "y",
                width = 0.28, linewidth = 0.5, colour = "grey40") +
  geom_point(aes(colour = sig), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "#33BB36", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "ns"),
                      drop = FALSE, name = NULL) +
  facet_wrap(~ factor) +
  coord_cartesian(xlim = xr) +
  scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +
  labs(x = "Mismatch differences / Mb", y = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background  = element_blank(),
        strip.text        = element_text(face = "bold"),
        panel.spacing     = unit(1.1, "lines"),
        legend.direction  = "vertical",
        legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
        legend.key.height = unit(0.9, "lines"),
        legend.margin     = margin(2, 4, 2, 4)) +
  legend_inside


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
  geom_smooth(method = "lm", se = TRUE, colour = "#33BB36", fill = "#fdd0a2",
              linewidth = 0.7, alpha = 0.4) +
  geom_point(size = 2.8, colour = "#9B77BB") +
  geom_text_repel(aes(label = channel), size = 3, seed = 1,
                  min.segment.length = 0, box.padding = 0.4, colour = "grey25") +
  annotate("text", x = -Inf, y = Inf, label = r_lab, parse = TRUE,
           hjust = -0.3, vjust = 1.6, size = 4, colour = "grey20") +
  labs(x = "Age effect (Old - Young, missmatch / Mb)",
       y = "Condition effect (Deletion - Control, missmatch / Mb)") +
  theme_classic(base_size = 12)



(pA | pB) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A", title = FILTER)

cat(sprintf("\nPearson r = %.2f (p = %.4f) | Spearman rho = %.2f | permutation p = %.4f\n",
            ct_p$estimate, ct_p$p.value, ct_s$estimate, p_perm))






diag <- channel_summary %>%
  transmute(channel, baseline = rate_per_Mb, age_d = age_diff, status_d = status_diff)

cat(sprintf("\nr(age, baseline) = %.2f | r(status, baseline) = %.2f | r(age, status) = %.2f\n",
            cor(diag$age_d, diag$baseline), cor(diag$status_d, diag$baseline),
            cor(diag$age_d, diag$status_d)))

# remove the proportional-to-baseline part (regression through the origin),
# then correlate what is left
diag <- diag %>%
  mutate(age_res    = resid(lm(age_d    ~ 0 + baseline)),
         status_res = resid(lm(status_d ~ 0 + baseline)))

ct_adj <- cor.test(diag$age_res, diag$status_res)
cat(sprintf("Baseline-adjusted: r = %.2f (p = %.3f) | Spearman rho = %.2f\n",
            ct_adj$estimate, ct_adj$p.value,
            cor(diag$age_res, diag$status_res, method = "spearman")))









#### ---- Spectrum bar plots (Fig S4/S5 style) ---------------------------- ####
# fixed channel order for this figure; separate name so the forest plot's
# abundance-based channel_order is not overwritten
spec_order <- c("C>A","C>G","C>U","U>A","U>C","U>G",
                "A>C","A>G","A>U","G>A","G>C","G>U")

spec_long <- spec_dat %>%
  select(sample_id, Age, Mutation_Status, depth_sum_genespan, all_of(SPEC_COLS)) %>%
  pivot_longer(all_of(SPEC_COLS), names_to = "channel", values_to = "count") %>%
  mutate(channel = factor(sub("^spec_", "", channel), levels = spec_order),
         rate    = count / depth_sum_genespan * 1e6)

DW  <- 0.8      # dodge width; 2 groups -> bar centres at +/- DW/4
OFF <- DW / 4


panel_specs <- list(
  list(group = "Age", p_col = "age_p_fdr",
       title = "Young vs Old",
       pal   = c("Young" = "#9ecae1", "Old" = "#3182bd")),
  list(group = "Mutation_Status", p_col = "status_p_fdr",
       title = paste(levels(spec_dat$Mutation_Status), collapse = " vs "),
       pal   = setNames(c("#2c7fb8", "#d95f0e"), levels(spec_dat$Mutation_Status)))
)
spec_panels <- list()

for (ps in panel_specs) {
  d    <- spec_long %>% mutate(grp = .data[[ps$group]])
  yrng <- diff(range(d$rate))
  
  # brackets for FDR-significant channels
  brk <- channel_summary %>%
    transmute(channel = factor(channel, levels = spec_order),
              p_fdr   = .data[[ps$p_col]]) %>%
    filter(p_fdr < 0.05) %>%
    left_join(d %>% group_by(channel) %>% summarise(ytop = max(rate), .groups = "drop"),
              by = "channel") %>%
    mutate(idx  = as.integer(channel),
           x1   = idx - OFF, x2 = idx + OFF,
           y    = ytop + 0.04 * yrng, tick = 0.015 * yrng,
           star = as.character(cut(p_fdr, c(-Inf, .001, .01, .05, Inf),
                                   c("***", "**", "*", ""))))
  
  p <- ggplot(d, aes(channel, rate, fill = grp)) +
    stat_summary(fun = mean, geom = "col", position = position_dodge(DW),
                 width = 0.7, colour = "grey30", linewidth = 0.3) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = DW),
               shape = 21, size = 1.1, alpha = 0.5, stroke = 0.3,
               colour = "grey20", show.legend = FALSE) +
    scale_fill_manual(values = ps$pal, name = NULL) +
    labs(x = NULL, y = "Mismatches / Mb", title = ps$title) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          plot.title  = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(face = "bold"))
  
  if (nrow(brk)) {
    p <- p +
      geom_segment(data = brk, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = brk, aes(x = x1, xend = x1, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = brk, aes(x = x2, xend = x2, y = y, yend = y - tick),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = brk, aes(x = idx, y = y, label = star),
                inherit.aes = FALSE, vjust = -0.1, size = 3.2)
  }
  
  spec_panels[[ps$group]] <- p
}


(spec_panels$Age / spec_panels$Mutation_Status) +
  plot_annotation(tag_levels = "A", title = FILTER)





diag <- channel_summary %>%
  transmute(channel, baseline = rate_per_Mb, age_d = age_diff, status_d = status_diff)

cat(sprintf("\nr(age, baseline) = %.2f | r(status, baseline) = %.2f | r(age, status) = %.2f\n",
            cor(diag$age_d, diag$baseline), cor(diag$status_d, diag$baseline),
            cor(diag$age_d, diag$status_d)))

# remove the proportional-to-baseline part (regression through the origin),
# then correlate what is left
diag <- diag %>%
  mutate(age_res    = resid(lm(age_d    ~ 0 + baseline)),
         status_res = resid(lm(status_d ~ 0 + baseline)))

ct_adj <- cor.test(diag$age_res, diag$status_res)
cat(sprintf("Baseline-adjusted: r = %.2f (p = %.3f) | Spearman rho = %.2f\n",
            ct_adj$estimate, ct_adj$p.value,
            cor(diag$age_res, diag$status_res, method = "spearman")))






#### ---- pN/pS: weighted and unweighted ------------------------------ ####
PNPS_COLS <- c("Unweighted" = "dNdS", "Weighted" = "dNdS_weighted")

pnps_dat <- DATAFRAME %>%
  filter(!is.na(Age)) %>%
  select(sample_id, Age, Mutation_Status, Deletion, all_of(PNPS_COLS)) %>%   # renames to Unweighted / Weighted
  pivot_longer(all_of(names(PNPS_COLS)), names_to = "measure", values_to = "pNpS")

n_bad <- sum(!is.finite(pnps_dat$pNpS))
if (n_bad) warning(n_bad, " non-finite pN/pS values dropped (e.g. pS = 0)")
pnps_dat <- pnps_dat %>% filter(is.finite(pNpS))

pnps_desc <- bind_rows(
  pnps_dat %>%
    group_by(measure, Age, Mutation_Status) %>%
    summarise(n = n(), mean = mean(pNpS), sd = sd(pNpS), .groups = "drop") %>%
    mutate(across(c(Age, Mutation_Status), as.character)),
  pnps_dat %>%
    group_by(measure) %>%
    summarise(n = n(), mean = mean(pNpS), sd = sd(pNpS), .groups = "drop") %>%
    mutate(Age = "All", Mutation_Status = "All")
) %>%
  mutate(sem = sd / sqrt(n),
         lcl = mean - qt(0.975, n - 1) * sem,
         ucl = mean + qt(0.975, n - 1) * sem)

pnps_table <- pnps_desc %>%
  mutate(value = sprintf("%.3f [%.3f, %.3f]", mean, lcl, ucl)) %>%
  select(Age, Mutation_Status, measure, value, n) %>%
  pivot_wider(names_from = measure, values_from = c(value, n),
              names_glue = "{measure}_{.value}")

show(pnps_table, "pN/pS by group: mean [95% CI]")


# Model per pN/pS version: raw (additive) scale, Deletion line as random effect.
# Age:    Old - Young,        pooled over Control + Deletion
# Status: Deletion - Control, pooled over Young + Old
# Weighted only: overall pN/pS (all samples) tested against 1.
pnps_results <- tibble()

for (lab in names(PNPS_COLS)) {
  d <- pnps_dat %>% filter(measure == lab)
  m <- glmmTMB(pNpS ~ Age + Mutation_Status + (1 | Deletion),
               family = gaussian, data = d)
  
  for (fac in c("Age", "Mutation_Status")) {
    r <- emmeans(m, fac) %>% contrast("revpairwise") %>%
      summary(infer = TRUE) %>% as.data.frame()
    lcl <- intersect(c("lower.CL", "asymp.LCL"), names(r))[1]
    ucl <- intersect(c("upper.CL", "asymp.UCL"), names(r))[1]
    
    r_out <- r %>%
      transmute(measure = lab, test = "Group difference",
                comparison = as.character(contrast),
                estimate, lcl = .data[[lcl]], ucl = .data[[ucl]], p = p.value)
    
    show(r_out %>% mutate(p = format.pval(p, digits = 3)),
         paste0(lab, " pN/pS: ", r_out$comparison), digits = 4)
    pnps_results <- bind_rows(pnps_results, r_out)
  }
  
  if (lab == "Weighted") {
    v1 <- emmeans(m, ~ Age * Mutation_Status) %>%
      summary(infer = TRUE, null = 1, adjust = "holm") %>% as.data.frame()
    lcl <- intersect(c("lower.CL", "asymp.LCL"), names(v1))[1]
    ucl <- intersect(c("upper.CL", "asymp.UCL"), names(v1))[1]
    
    v1_out <- v1 %>%
      transmute(measure = lab, test = "Group pN/pS vs 1 (Holm)",
                comparison = paste(Age, Mutation_Status, "vs 1"),
                estimate = emmean, lcl = .data[[lcl]], ucl = .data[[ucl]], p = p.value)
    
    show(v1_out %>% mutate(p = format.pval(p, digits = 3)),
         "Weighted pN/pS: each group vs 1", digits = 4)
    pnps_results <- bind_rows(pnps_results, v1_out)
  }
}




show(pnps_results %>% mutate(p = format.pval(p, digits = 3)),
     "All pN/pS tests", digits = 4)



#### ---- pN/pS boxplots (Fig S2/S3 style) ------------------------------- ####
pal_pnps <- setNames(c("#2c7fb8", "#d95f0e"), levels(pnps_dat$Mutation_Status))

pnps_specs <- list(
  list(measure = "Weighted",   ylab = "Weighted dN/dS",   show_null1 = TRUE),
  list(measure = "Unweighted", ylab = "Unweighted pN/pS", show_null1 = FALSE)  # no line at 1: spectrum bias moves the neutral point
)
pnps_panels <- list()

for (ps in pnps_specs) {
  df <- pnps_dat %>% filter(measure == ps$measure)
  
  # zoom window: median +/- 0.2 (view only; stats use all data)
  ctr <- median(df$pNpS)
  ylo <- ctr - 0.2
  yhi <- ctr + 0.2
  if (ps$show_null1) { ylo <- min(ylo, 0.97); yhi <- max(yhi, 1.03) }   # keep the neutral line visible
  
  p <- ggplot(df, aes(Age, pNpS, fill = Mutation_Status)) +
    geom_boxplot(width = 0.55, position = position_dodge(0.62),
                 outlier.shape = 21, outlier.size = 1.2, outlier.alpha = 0.5,
                 colour = "grey25", linewidth = 0.5) +
    scale_fill_manual(values = pal_pnps, name = NULL) +
    labs(x = "Age", y = ps$ylab, title = ps$measure) +
    coord_cartesian(ylim = c(ylo, yhi)) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 11, face = "bold"))
  
  if (ps$show_null1) {
    p <- p +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0.5, y = 1, label = "neutral (= 1)",
               hjust = 0, vjust = -0.5, size = 3, colour = "grey40")
  }
  
  pnps_panels[[ps$measure]] <- p
}


(pnps_panels$Weighted | pnps_panels$Unweighted) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", title = FILTER) &
  theme(legend.position = "bottom")



filt_check <- tibble()

for (flt in available_filters) {
  s <- read_tsv(file.path(RESULTS_DIR, paste0(flt, ".summary.tsv")), show_col_types = FALSE) %>%
    filter(sample_id %in% meta$sample_id, is.finite(dNdS_weighted), is.finite(dNdS))
  
  filt_check <- bind_rows(filt_check, tibble(
    filter      = flt,
    mean_w      = mean(s$dNdS_weighted), median_w = median(s$dNdS_weighted),
    mean_u      = mean(s$dNdS),          median_u = median(s$dNdS),
    median_syn  = if ("synonymous_mutations" %in% names(s)) median(s$synonymous_mutations) else NA,
    median_mis  = if ("missense_mutations"   %in% names(s)) median(s$missense_mutations)   else NA
  ))
}

show(filt_check, "pN/pS across filters: mean vs median, with typical counts per sample")

# columns needed for a pooled (sum N / sum S) pN/pS
grep("syn|miss|dn|ds|sites|N_|S_", names(DATAFRAME), value = TRUE, ignore.case = TRUE)




show(DATAFRAME %>%
       summarise(depth_per_pos_nuc   = median(depth_sum_total / positions_covered),
                 depth_per_pos_mito  = median(depth_sum_total_mito / positions_covered_mito),
                 reads_per_allele_nuc  = median(reads_per_allele),
                 reads_per_allele_mito = median(reads_per_allele_mito)),
     "Depth and reads per allele: nuclear vs mito")




#### ---- Table S1: pN/pS by age and genotype, all filters -> Excel ------- ####
library(openxlsx)

TABLE_FILTERS <- c("Unfiltered" = "unfiltered", "5%" = "across_5pct",
                   "1%" = "across_1pct", "0.1%" = "across_0p1pct")
GROUP_ORDER   <- c("Young", "Old", "Wild type", "Deletion")

# same model as the main pN/pS analysis; groups are pooled marginal means
tab_long <- tibble()

for (flab in names(TABLE_FILTERS)) {
  d_f <- read_tsv(file.path(RESULTS_DIR, paste0(TABLE_FILTERS[[flab]], ".summary.tsv")),
                  show_col_types = FALSE) %>%
    inner_join(meta, by = "sample_id")
  
  for (lab in names(PNPS_COLS)) {
    d <- d_f %>% mutate(pNpS = .data[[PNPS_COLS[[lab]]]]) %>% filter(is.finite(pNpS))
    m <- glmmTMB(pNpS ~ Age + Mutation_Status + (1 | Deletion), family = gaussian, data = d)
    
    em <- bind_rows(
      emmeans(m, ~ Age) %>% summary(null = 1, infer = c(FALSE, TRUE)) %>%
        as.data.frame() %>% transmute(group = as.character(Age), emmean, SE, p.value),
      emmeans(m, ~ Mutation_Status) %>% summary(null = 1, infer = c(FALSE, TRUE)) %>%
        as.data.frame() %>% transmute(group = as.character(Mutation_Status), emmean, SE, p.value)
    ) %>%
      mutate(filter = flab, measure = lab,
             p_holm = p.adjust(p.value, "holm"))     # vs 1, Holm over the 4 groups
    
    tab_long <- bind_rows(tab_long, em)
  }
}

tab_long <- tab_long %>%
  mutate(group = recode(group, "Control" = "Wild type"),
         cell  = sprintf("%.3f \u00b1 %.3f", emmean, SE),
         bold  = measure == "Weighted" & p_holm < 0.05)

show(tab_long %>% select(filter, measure, group, cell, p_holm, bold), "Table S1 values")

## ---- wide layout ---------------------------------------------------------- #
tab_wide <- tab_long %>%
  mutate(col = paste(measure, group)) %>%
  select(filter, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  mutate(filter = factor(filter, names(TABLE_FILTERS))) %>%
  arrange(filter) %>%
  transmute(Filter = as.character(filter),
            !!!syms(paste("Unweighted", GROUP_ORDER)),
            spacer = "",
            !!!syms(paste("Weighted", GROUP_ORDER)))

## ---- write Excel ---------------------------------------------------------- #
wb <- createWorkbook()
addWorksheet(wb, "Table S1")

caption <- paste("Table S1. Unweighted and weighted ratios of nonsynonymous to synonymous",
                 "mismatches (pN/pS) by age and genotype (condition). Mean estimates and",
                 "standard errors are reported for unfiltered data, as well as for data filtered",
                 "for positions with mismatch frequencies above a threshold of 5, 1 and 0.1%.",
                 "Estimates in bold are significantly different from 1.")
writeData(wb, 1, caption, startRow = 1, startCol = 1)

writeData(wb, 1, t(c("", "pN/pS unweighted", "", "", "", "", "pN/pS weighted", "", "", "")),
          startRow = 3, colNames = FALSE)
mergeCells(wb, 1, cols = 2:5,  rows = 3)
mergeCells(wb, 1, cols = 7:10, rows = 3)
writeData(wb, 1, t(c("Filter", GROUP_ORDER, "", GROUP_ORDER)), startRow = 4, colNames = FALSE)
writeData(wb, 1, tab_wide, startRow = 5, colNames = FALSE)

addStyle(wb, 1, createStyle(textDecoration = "bold", halign = "center"),
         rows = 3:4, cols = 1:10, gridExpand = TRUE)
addStyle(wb, 1, createStyle(halign = "center"),
         rows = 5:(4 + nrow(tab_wide)), cols = 2:10, gridExpand = TRUE)

# bold the weighted estimates that differ from 1
b <- tab_long %>%
  filter(bold) %>%
  mutate(r = 4 + match(filter, names(TABLE_FILTERS)),
         c = ifelse(measure == "Weighted", 6, 1) + match(group, GROUP_ORDER))
if (nrow(b))
  addStyle(wb, 1, createStyle(textDecoration = "bold", halign = "center"),
           rows = b$r, cols = b$c, gridExpand = FALSE, stack = TRUE)

setColWidths(wb, 1, cols = 1:10, widths = c(12, rep(15, 4), 3, rep(15, 4)))
saveWorkbook(wb, "Table_S1_pNpS.xlsx", overwrite = TRUE)




