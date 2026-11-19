# Produce replication tables and figures.
# Usage: Rscript R/report.R [batch.rds] [out_prefix]
source("R/analysis.R")
source("R/compare.R")
suppressPackageStartupMessages({ library(ggplot2); library(mgcv) })

a <- commandArgs(TRUE)
batch_file <- if (length(a) >= 1) a[1] else "outputs/sim_main.rds"
pre <- if (length(a) >= 2) a[2] else "outputs/"
b <- readRDS(batch_file)
o <- readRDS("data/ggs_original.rds")
w_b <- enpp_weights(b$runs$enpp); w_o <- enpp_weights(o$runs$enpp)
col_ggs <- "#2a78d6"; col_rep <- "#eb6834"; col_real <- "#555555"
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "top"))

## ---- Tables -----------------------------------------------------------------
cmp <- compare_batch(b, o)
saveRDS(cmp, paste0(pre, "comparison.rds"))

se_tab <- function(ct) sprintf("%.3f (%.3f)", ct[, 1], ct[, 2])
t2b <- table2(b$runs, b$gov, w_b); t3b <- table3(b$runs, b$gov, w_b); t4b <- table4(b$runs, w_b)

md <- c("# Replication of Golder, Golder & Siegel (2012): results", "",
        sprintf("Replication batch: %s (%d party systems, %d governing-party observations).",
                batch_file, nrow(b$runs), nrow(b$gov)), "",
        "Columns:", "",
        "- **Published**: the values reported in the article.",
        "- **GGS data**: our analysis code run on the authors' deposited simulation output.",
        "- **Replication**: our R implementation of the model.", "",
        "## Headline statistics", "",
        "| Statistic | Published | GGS data | Replication | Diff (Rep − GGS) |",
        "|---|---:|---:|---:|---:|",
        with(cmp$headline, sprintf("| %s | %.3f | %.3f | %.3f | %+.3f |", stat, published,
                                   ggs_data, replication, replication - ggs_data)),
        "", "## Replication regressions with clustered SEs", "",
        "Table 2, unweighted:", "", "```", capture.output(print(t2b$unweighted$ct)), "```",
        "Table 2, weighted:", "", "```", capture.output(print(t2b$weighted$ct)), "```",
        "Table 3, formateur:", "", "```", capture.output(print(t3b$formateur$ct)),
        sprintf("R2 = %.3f", t3b$formateur$r2), "```",
        "Table 3, parameters:", "", "```", capture.output(print(t3b$params$ct)),
        sprintf("R2 = %.3f", t3b$params$r2), "```",
        "Table 4, OLS:", "", "```", capture.output(print(coeftest(t4b$ols))),
        sprintf("R2 = %.3f", summary(t4b$ols)$r.squared), "```",
        "Table 4, Cox (Efron ties):", "", "```", capture.output(print(coef(summary(t4b$cox)))), "```",
        "", "## Conditional on the number of parties m (unweighted)", "",
        "GGS data:", "", "```", capture.output(print(cmp$by_m$ggs)), "```",
        "Replication:", "", "```", capture.output(print(cmp$by_m$rep)), "```",
        "", "## Conditional on beta (unweighted)", "",
        "GGS data:", "", "```", capture.output(print(cmp$by_beta$ggs)), "```",
        "Replication:", "", "```", capture.output(print(cmp$by_beta$rep)), "```")
writeLines(md, paste0(pre, "results.md"))

## ---- Figure 1: single run --------------------------------------------------
if (file.exists("outputs/fig1_run.rds")) {
  f1 <- read.delim("data/original/GGS_JOP_GovForm_Data_Fig1.tab", check.names = FALSE)
  r1 <- readRDS("outputs/fig1_run.rds")
  d1 <- rbind(
    data.frame(src = "GGS data", tick = f1$tick, coef = "Slope", v = f1[[29]]),
    data.frame(src = "GGS data", tick = f1$tick, coef = "Intercept", v = f1[[30]]),
    data.frame(src = "Replication", tick = seq_along(r1$run_slope), coef = "Slope", v = r1$run_slope),
    data.frame(src = "Replication", tick = seq_along(r1$run_int), coef = "Intercept", v = r1$run_int))
  p1 <- ggplot(subset(d1, tick >= 10), aes(tick, v, colour = src)) +
    geom_line(linewidth = 0.6) + geom_vline(xintercept = 2500, linetype = 2, colour = "grey50") +
    facet_wrap(~coef, scales = "free_y") +
    scale_colour_manual(values = c("GGS data" = col_ggs, "Replication" = col_rep), name = NULL) +
    labs(x = "Government formations", y = "Cumulative coefficient",
         title = "Figure 1: cumulative Gamson regression for one party system")
  ggsave(paste0(pre, "fig1.png"), p1, width = 8, height = 3.6, dpi = 150)
}

## ---- Figure 2: government type vs beta (weighted GAMs) ----------------------
gam_curve <- function(runs, w, type_set, src) {
  y <- as.numeric(runs$type %in% type_set); ok <- !is.na(runs$type)
  dd <- data.frame(y = y[ok], beta = runs$beta[ok], w = w[ok] / mean(w[ok]))
  g <- gam(y ~ s(beta), data = dd, weights = w)
  nd <- data.frame(beta = seq(0, 1, by = 0.01))
  p <- predict(g, nd, type = "terms", se.fit = TRUE)
  data.frame(src = src, beta = nd$beta, fit = p$fit[, 1], se = p$se.fit[, 1])
}
panels <- list("(A) Single-party minority" = "spmin", "(B) Minority coalition" = "mincoal",
               "(C) Minimal winning coalition" = "mwc", "(D) Surplus majority" = "surplus")
d2 <- do.call(rbind, lapply(names(panels), function(pn) {
  rbind(cbind(panel = pn, gam_curve(o$runs, w_o, panels[[pn]], "GGS data")),
        cbind(panel = pn, gam_curve(b$runs, w_b, panels[[pn]], "Replication")))
}))
p2 <- ggplot(d2, aes(beta, fit, colour = src, fill = src)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.7) + facet_wrap(~panel) + coord_cartesian(ylim = c(-0.1, 0.1)) +
  scale_colour_manual(values = c("GGS data" = col_ggs, "Replication" = col_rep), name = NULL) +
  scale_fill_manual(values = c("GGS data" = col_ggs, "Replication" = col_rep), name = NULL) +
  labs(x = expression(beta), y = "Change in proportion (GAM smooth term)",
       title = "Figure 2: policy weight and government type (weighted)")
ggsave(paste0(pre, "fig2.png"), p2, width = 8, height = 6, dpi = 150)

## ---- Figure 3: delay distribution ------------------------------------------
real <- read.delim("data/original/GGS_JOP_GovForm_Data_real1.tab", check.names = FALSE,
                   fileEncoding = "latin1")
dens <- function(x, w, src) { dd <- density(x, weights = w / sum(w), bw = 5, from = 0, to = 230)
  data.frame(src = src, x = dd$x, y = dd$y) }
rd <- real$bargaining_duration[!is.na(real$bargaining_duration)]
d3 <- rbind(dens(rd, rep(1, length(rd)), "Real world (West Europe)"),
            dens(2.08 * o$runs$delay, w_o, "GGS data (weighted)"),
            dens(2.08 * b$runs$delay, w_b, "Replication (weighted)"))
p3 <- ggplot(d3, aes(x, y, colour = src)) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Real world (West Europe)" = col_real,
                                 "GGS data (weighted)" = col_ggs,
                                 "Replication (weighted)" = col_rep), name = NULL) +
  labs(x = "Bargaining delay (days; simulated periods x 2.08)", y = "Density",
       title = "Figure 3: distribution of bargaining delays")
ggsave(paste0(pre, "fig3.png"), p3, width = 7, height = 4, dpi = 150)

## ---- Distribution of delay (simulated periods), side by side ----------------
qs <- c(.1, .25, .5, .75, .9)
cat("Delay quantiles GGS:", quantile(o$runs$delay, qs), "\n")
cat("Delay quantiles Rep:", quantile(b$runs$delay, qs), "\n")
cat("KS test of delay distributions: D =",
    round(suppressWarnings(ks.test(o$runs$delay, b$runs$delay))$statistic, 4), "\n")
message("Wrote ", pre, "results.md and figures")
