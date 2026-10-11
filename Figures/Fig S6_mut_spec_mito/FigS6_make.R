# ============================================================================
# Figure S6 — Mutation spectra (12 channels), MITOCHONDRIAL transcripts
#   TOP    : Young vs Old          (bracket = FDR-sig per-channel Age effect)
#   BOTTOM : Wild type vs Deletion (bracket = FDR-sig per-channel Status effect)
#
# Bars = mean per-channel rate (errors/Mb); transparent dots = per-sample rates;
# brackets + stars = channels whose rate differs significantly between the two
# groups (per-channel additive NB-GLMM, depth-offset, FDR-corrected).
#
# ============================================================================

library(here)
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)

## ---- CHOOSE THE FILTER ---------------------------------------------------- #
FILTER <- "unfiltered"

FOLDER <- dplyr::case_when(
  FILTER == "unfiltered" ~ "Fig S6_mut_spec_mito",
  TRUE                        ~ paste0("Fig S_mut_spec_mito_", FILTER)
)

RESULTS_DIR <- here("Data_Results")
META_FILE   <- here("Data_Results", "metadata_minimal.txt")
OUT_DIR     <- here( FOLDER)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPEC_RE    <- "^spec_[ACGU]>[ACGU]_mito$"
DEPTH_COL  <- "depth_sum_total_mito"          # <- mito depth, not nuclear

# strips BOTH the prefix and the _mito suffix -> "A>C"
strip_spec <- function(x) sub("_mito$", "", sub("^spec_", "", x))

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
  mutate(depth = .data[[DEPTH_COL]]) %>%
  filter(!is.na(depth), depth > 0)

SPEC_COLS     <- grep(SPEC_RE, names(DATAFRAME), value = TRUE)
channel_names <- strip_spec(SPEC_COLS)                  # "A>C", "A>G", ...

if (length(SPEC_COLS) != 12) {
  warning("Expected 12 mito spectrum columns, found ", length(SPEC_COLS), ".")
}

channel_order <- c("C>A","C>G","C>U","U>A","U>C","U>G",
                   "A>C","A>G","A>U","G>A","G>C","G>U")

## ---- per-channel additive models: Age & Status significance --------------- #
add_p <- function(m, spec) {
  emmeans(m, spec, offset = log(1e6)) %>% regrid() %>%
    contrast("revpairwise") %>% summary() %>% as.data.frame() %>% pull(p.value) %>% .[1]
}

fit_channel <- function(ch) {
  d <- DATAFRAME %>% mutate(y = .data[[paste0("spec_", ch, "_mito")]])
  # channels that are all-zero (or near) will not fit; return NA rather than stop
  if (sum(d$y, na.rm = TRUE) == 0) {
    return(tibble(channel = ch, age_p = NA_real_, status_p = NA_real_))
  }
  m <- try(glmmTMB(y ~ Age + Mutation_Status + offset(log(depth)) + (1 | Deletion),
                   family = nbinom2, data = d), silent = TRUE)
  if (inherits(m, "try-error")) {
    warning("Model failed for channel ", ch)
    return(tibble(channel = ch, age_p = NA_real_, status_p = NA_real_))
  }
  tibble(channel = ch, age_p = add_p(m, ~ Age), status_p = add_p(m, ~ Mutation_Status))
}

sig <- map_dfr(channel_names, fit_channel) %>%
  mutate(age_p_fdr    = p.adjust(age_p, "fdr"),
         status_p_fdr = p.adjust(status_p, "fdr"))

## ---- long per-sample rates ------------------------------------------------ #
spec_long <- DATAFRAME %>%
  select(sample_id, Age, Mutation_Status, depth, all_of(SPEC_COLS)) %>%
  pivot_longer(all_of(SPEC_COLS), names_to = "channel", values_to = "count") %>%
  mutate(channel = strip_spec(channel),
         rate    = count / depth * 1e6)

stopifnot(all(unique(spec_long$channel) %in% channel_order))   # catches the NA bug early
spec_long$channel <- factor(spec_long$channel, levels = channel_order)

## ---- panel builder -------------------------------------------------------- #
DW  <- 0.8            # dodge width; 2 groups -> bar centres at +/- DW/4
OFF <- DW / 4
stars <- function(p) cut(p, c(-Inf,.001,.01,.05,Inf), c("***","**","*",""))

spectrum_panel <- function(group_col, pal, p_col, title) {
  d <- spec_long %>% mutate(grp = .data[[group_col]])
  
  # bracket positions for FDR-significant channels
  sig_ch <- sig %>%
    transmute(channel = factor(channel, channel_order),
              p_fdr   = .data[[p_col]]) %>%
    filter(!is.na(p_fdr), p_fdr < 0.05)
  
  ytop <- d %>% group_by(channel) %>% summarise(ytop = max(rate), .groups = "drop")
  yrng <- diff(range(d$rate))
  brk  <- sig_ch %>% left_join(ytop, by = "channel") %>%
    mutate(idx  = as.integer(channel),
           x1   = idx - OFF, x2 = idx + OFF,
           y    = ytop + 0.04 * yrng, tick = 0.015 * yrng,
           star = as.character(stars(p_fdr)))
  
  p <- ggplot(d, aes(channel, rate, fill = grp)) +
    stat_summary(fun = mean, geom = "col", position = position_dodge(DW),
                 width = 0.7, colour = "grey30", linewidth = 0.3) +
    geom_point(aes(fill = grp),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = DW),
               shape = 21, size = 1.1, alpha = 0.5, stroke = 0.3,
               colour = "grey20", show.legend = FALSE) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_colour_manual(values = pal, guide = "none") +
    # headroom so the brackets are not clipped by the panel edge
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = NULL, y = "Rate (errors / Mb)", title = title) +
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
  p
}

## ---- build the two stacked panels ----------------------------------------- #
pal_age  <- c("Young" = "#9ecae1", "Old" = "#3182bd")
pal_stat <- c("Wild type" = "#2c7fb8", "Deletion" = "#d95f0e")

pTop <- spectrum_panel("Age",             pal_age,  "age_p_fdr",    "Young vs Old")
pBot <- spectrum_panel("Mutation_Status", pal_stat, "status_p_fdr", "Wild type vs Deletion")

figS <- (pTop / pBot) + plot_annotation(tag_levels = "A")

figS

base <- paste0("Fig_mut_spectra_mito_", FILTER)
ggsave(file.path(OUT_DIR, paste0(base, ".pdf")), figS, width = 9, height = 7.5)
ggsave(file.path(OUT_DIR, paste0(base, ".png")), figS, width = 9, height = 7.5, dpi = 300)
ggsave(file.path(OUT_DIR, paste0(base, ".svg")), figS, width = 9, height = 7.5, dpi = 300)


## ---- accompanying tables -------------------------------------------------- #
group_means <- spec_long %>%
  group_by(channel, Age, Mutation_Status) %>%
  summarise(mean_rate = mean(rate), .groups = "drop")
write_csv(group_means, file.path(OUT_DIR, paste0(base, "_group_means.csv")))
write_csv(sig %>% mutate(channel = factor(channel, channel_order)) %>% arrange(channel),
          file.path(OUT_DIR, paste0(base, "_significance.csv")))

message("Written to: ", OUT_DIR)


