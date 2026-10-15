##############################Michael Doerfler##################################
##############################2026-10-14########################################
#Figures for the K = 2 asymptotic-behaviour study.
#Reads sims/<study>/summary.rds (output of 3_summarize.R) and writes PDFs to report/figures/.

library(here)
library(dplyr)
library(ggplot2)

study <- "k2_asymptotics_pilot"          # switch to "k2_asymptotics" for the full run
fig_dir <- here::here("report", "figures", study)
perf  <- readRDS(here::here("sims", study, "summary.rds"))
design <- readRDS(here::here("sims", study, "design.rds"))
perf   <- left_join(perf, design %>% select(scenario_id, N1, N2, psi), by = "scenario_id")

# Methods to compare and how to label them ---------------------------------------
method_labels <- c(mle_wald     = "MLE, Wald",
                   mle_logit    = "MLE, logit",
                   bayes_median = "Bayes")

# Coverage plot -------------------------------------------------------------------
cov_df <- perf %>%
  filter(method %in% names(method_labels)) %>%
  mutate(method = factor(method, levels = names(method_labels), labels = method_labels),
         target = factor(target, levels = c("N1", "N2"),
                         labels = c("Cohort 1 (with anchor)", "Cohort 2")),
         lo = coverage - 1.96 * coverage_mcse,      # Monte Carlo error bars
         hi = coverage + 1.96 * coverage_mcse)

dodge <- position_dodge(width = 0.04)   # in log10 units: spreads the methods slightly apart

p_cov <- ggplot(cov_df, aes(x = N1, y = coverage, colour = method, shape = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
  geom_line(position = dodge) +
  geom_point(position = dodge, size = 1.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = dodge) +
  geom_vline(xintercept = 159, linetype = "dotted", colour = "grey60") +   # CRISP sample size
  scale_x_log10(breaks = c(50, 200, 1000, 5000, 20000), labels = scales::comma) +
  facet_wrap(~ target) +
  labs(x = expression(N^{(1)} ~ "(log scale)"),
       y = "Coverage of 95% interval",
       colour = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p_cov

ggsave(file.path(fig_dir, "coverage.pdf"), p_cov, width = 7, height = 3.8)

# Interval width bar plot ------------------------------------------------------------
# Width grows like sqrt(N), so raw widths span orders of magnitude across the grid.
# Plotting width relative to the true N (width / N) puts all scenarios on one scale;
# set use_relative <- FALSE for raw widths.
use_relative <- TRUE

width_df <- perf %>%
  filter(method %in% names(method_labels)) %>%
  mutate(true_N = ifelse(target == "N1", N1, N2),
         y      = if (use_relative) width / true_N else width,
         y_mcse = if (use_relative) width_mcse / true_N else width_mcse,
         method = factor(method, levels = names(method_labels), labels = method_labels),
         target = factor(target, levels = c("N1", "N2"),
                         labels = c("Cohort 1 (with anchor)", "Cohort 2")),
         N1_f   = factor(N1, levels = sort(unique(N1))))     # bars need a discrete x

p_width <- ggplot(width_df, aes(x = N1_f, y = y, fill = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_errorbar(aes(ymin = y - 1.96 * y_mcse, ymax = y + 1.96 * y_mcse),
                position = position_dodge(width = 0.85), width = 0.3, linewidth = 0.3) +
  facet_wrap(~ target, ncol = 1) +
  labs(x = expression(N^{(1)}),
       y = if (use_relative) "Mean interval width / N" else "Mean interval width",
       fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_width

ggsave(file.path(fig_dir, "width.pdf"),    p_width, width = 7, height = 3.8)
