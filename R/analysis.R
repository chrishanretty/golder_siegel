# Analysis functions: reproduce Tables 1-4 and delay statistics of
# Golder, Golder & Siegel (2012) from a set of simulated party systems.
#
# Input format (common to the original GGS data and our replication):
#   runs:  data.frame, one row per party system, with columns
#          run, m, beta, enpp, diversity, delay, formateur (1-based party id),
#          type (one of "spmaj","spmin","mincoal","mwc","surplus")
#   gov:   data.frame, one row per *governing* party in the final government,
#          with columns run, party, seat_contrib (seat share of government),
#          cab_share, formateur (0/1)

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(survival)
})

gov_types <- c("spmaj", "mwc", "mincoal", "spmin", "surplus")
gov_type_labels <- c(spmaj = "Single Party Majority",
                     mwc = "Minimal Winning Coalition",
                     mincoal = "Minority Coalition",
                     spmin = "Single Party Minority",
                     surplus = "Surplus Majority")

## ---------------------------------------------------------------------------
## Real-world ENPP weights
## GGS use np::npudens with a likelihood-cross-validated bandwidth
## (npudensbw default: Gaussian kernel, cv.ml). We reimplement that here.
## ---------------------------------------------------------------------------
real_world_enpp <- function(path = "data/original/GGS_JOP_GovForm_Data_real1.tab") {
  r <- read.delim(path, check.names = FALSE, fileEncoding = "latin1")
  r$enpp[!is.na(r$enpp)]
}

bw_cv_ml <- function(x) {
  n <- length(x)
  D2 <- outer(x, x, "-")^2
  nll <- function(logh) {
    h <- exp(logh)
    K <- dnorm(sqrt(D2) / h) / h
    diag(K) <- 0
    -sum(log(rowSums(K) / (n - 1)))
  }
  exp(optimize(nll, log(c(0.01, 5)))$minimum)
}

enpp_weights <- function(enpp_sim, enpp_real = real_world_enpp()) {
  h <- bw_cv_ml(enpp_real)
  vapply(enpp_sim, function(e) mean(dnorm((e - enpp_real) / h)) / h, numeric(1))
}

## ---------------------------------------------------------------------------
## Read original GGS simulated data into the common format
## ---------------------------------------------------------------------------
read_ggs_original <- function(path = "data/original/GGS_JOP_GovForm_Data.tab") {
  d <- read.delim(path, check.names = FALSE)
  n <- nrow(d)
  type_cols <- c(spmaj = 68, spmin = 69, mincoal = 70, mwc = 71, surplus = 72, unity = 73)
  tm <- as.matrix(d[, type_cols])
  colnames(tm) <- names(type_cols)
  type <- names(type_cols)[max.col(tm, ties.method = "first")]
  type[rowSums(tm) == 0] <- NA
  type[type == "unity"] <- "surplus"     # GGS merge unity into surplus
  runs <- data.frame(run = seq_len(n), m = d$numParties, beta = d[[2]],
                     enpp = d[[24]], diversity = d[[25]], delay = d[[66]],
                     formateur = d[[3]] + 1, type = type,
                     n_gov = d[[67]])
  contrib <- as.matrix(d[, 46:55]); cab <- as.matrix(d[, 56:65])
  idx <- which(!is.na(cab) & cab > 0, arr.ind = TRUE)
  gov <- data.frame(run = idx[, 1], party = idx[, 2],
                    seat_contrib = contrib[idx], cab_share = cab[idx])
  gov$formateur <- as.integer(gov$party == runs$formateur[gov$run])
  gov <- gov[order(gov$run, gov$party), ]
  list(runs = runs, gov = gov)
}

## ---------------------------------------------------------------------------
## Tables
## ---------------------------------------------------------------------------
table1 <- function(runs, w) {
  ok <- !is.na(runs$type)
  sapply(gov_types, function(t) c(
    unweighted = mean(runs$type[ok] == t),
    weighted   = weighted.mean(runs$type[ok] == t, w[ok])))
}

clustered_lm <- function(formula, data, weights = NULL, cluster) {
  data$.w <- if (is.null(weights)) 1 else weights
  fit <- lm(formula, data = data, weights = .w)
  V <- vcovCL(fit, cluster = cluster, type = "HC0", cadjust = FALSE)
  list(fit = fit, ct = coeftest(fit, vcov. = V), r2 = summary(fit)$r.squared,
       n = nobs(fit))
}

table2 <- function(runs, gov, w) {
  g <- gov; g$w <- w[g$run]
  list(unweighted = clustered_lm(cab_share ~ seat_contrib, g, NULL, g$run),
       weighted   = clustered_lm(cab_share ~ seat_contrib, g, g$w, g$run))
}

table3 <- function(runs, gov, w) {
  g <- gov; g$w <- w[g$run]
  g$beta <- runs$beta[g$run]; g$m <- runs$m[g$run]
  g$enpp <- runs$enpp[g$run]; g$diversity <- runs$diversity[g$run]
  list(formateur = clustered_lm(cab_share ~ seat_contrib + formateur, g, g$w, g$run),
       params = clustered_lm(cab_share ~ seat_contrib * beta + seat_contrib * m +
                               seat_contrib * enpp + seat_contrib * diversity,
                             g, g$w, g$run))
}

table4 <- function(runs, w) {
  r <- runs; r$w <- w
  ols <- lm(delay ~ beta + m + enpp + diversity, data = r, weights = w)
  r$event <- as.integer(r$delay < 100)
  cox <- coxph(Surv(delay, event) ~ beta + m + enpp + diversity, data = r,
               weights = w, ties = "efron")
  list(ols = ols, cox = cox)
}

delay_stats <- function(runs, w) {
  c(mean_unweighted = mean(runs$delay), mean_weighted = weighted.mean(runs$delay, w),
    median = median(runs$delay), p_censored = mean(runs$delay == 100),
    p_one = mean(runs$delay == 1))
}

## Compact named vector of headline statistics (for comparing two datasets)
headline <- function(runs, gov, w = enpp_weights(runs$enpp)) {
  t1 <- table1(runs, w)
  t2 <- table2(runs, gov, w)
  t3 <- table3(runs, gov, w)
  t4 <- table4(runs, w)
  ds <- delay_stats(runs, w)
  c(setNames(t1["unweighted", ], paste0("T1u_", colnames(t1))),
    setNames(t1["weighted", ], paste0("T1w_", colnames(t1))),
    T2u_slope = unname(coef(t2$unweighted$fit)[2]), T2u_int = unname(coef(t2$unweighted$fit)[1]),
    T2w_slope = unname(coef(t2$weighted$fit)[2]), T2w_int = unname(coef(t2$weighted$fit)[1]),
    T2_N = t2$weighted$n,
    T3a_slope = unname(coef(t3$formateur$fit)["seat_contrib"]),
    T3a_int = unname(coef(t3$formateur$fit)["(Intercept)"]),
    T3a_form = unname(coef(t3$formateur$fit)["formateur"]),
    setNames(coef(t3$params$fit), paste0("T3b_", names(coef(t3$params$fit)))),
    setNames(coef(t4$ols), paste0("T4ols_", names(coef(t4$ols)))),
    T4ols_R2 = summary(t4$ols)$r.squared,
    setNames(coef(t4$cox), paste0("T4cox_", names(coef(t4$cox)))),
    delay_mean_u = unname(ds["mean_unweighted"]), delay_mean_w = unname(ds["mean_weighted"]),
    delay_median = unname(ds["median"]), delay_p100 = unname(ds["p_censored"]),
    delay_p1 = unname(ds["p_one"]),
    mean_n_gov = nrow(gov) / nrow(runs))
}

## Published values from the article (for the comparison table)
published <- c(
  T1u_spmaj = .086, T1u_mwc = .196, T1u_mincoal = .461, T1u_spmin = .150, T1u_surplus = .107,
  T1w_spmaj = .157, T1w_mwc = .257, T1w_mincoal = .254, T1w_spmin = .230, T1w_surplus = .102,
  T2u_slope = .744, T2u_int = .102, T2w_slope = .771, T2w_int = .117, T2_N = 125996,
  T3a_slope = .714, T3a_int = .107, T3a_form = .076,
  `T3b_(Intercept)` = .084, T3b_seat_contrib = 1.099, T3b_beta = -.010, T3b_m = .027,
  T3b_enpp = -.030, T3b_diversity = .058, `T3b_seat_contrib:beta` = .040,
  `T3b_seat_contrib:m` = -.098, `T3b_seat_contrib:enpp` = .046,
  `T3b_seat_contrib:diversity` = -.130,
  `T4ols_(Intercept)` = 84.8, T4ols_beta = 19.5, T4ols_m = 8.82, T4ols_enpp = -20.3,
  T4ols_diversity = -40.7, T4ols_R2 = .21,
  T4cox_beta = -.82, T4cox_m = -.27, T4cox_enpp = .53, T4cox_diversity = 1.42,
  delay_mean_u = 29.13, delay_mean_w = 29.84)
