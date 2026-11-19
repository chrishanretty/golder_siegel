# Compare a replication batch with the original GGS simulated data,
# overall (headline statistics) and conditionally on m and on beta.
source("R/analysis.R")

compare_batch <- function(b, o = readRDS("data/ggs_original.rds"), digits = 3) {
  w_b <- enpp_weights(b$runs$enpp); w_o <- enpp_weights(o$runs$enpp)
  hb <- headline(b$runs, b$gov, w_b); ho <- headline(o$runs, o$gov, w_o)
  nm <- names(published)
  head_tab <- data.frame(stat = nm, published = published[nm],
                         ggs_data = ho[nm], replication = hb[nm], row.names = NULL)
  by_m <- function(r, g) {
    g$m <- r$m[g$run]
    gl <- split(g, g$m)
    sl <- sapply(gl, function(x) coef(lm(cab_share ~ seat_contrib, x))[2])
    rbind(delay = tapply(r$delay, r$m, mean),
          p100 = tapply(r$delay == 100, r$m, mean),
          n_gov = tapply(r$n_gov, r$m, mean),
          spmaj = tapply(r$type == "spmaj", r$m, mean, na.rm = TRUE),
          spmin = tapply(r$type == "spmin", r$m, mean, na.rm = TRUE),
          mincoal = tapply(r$type == "mincoal", r$m, mean, na.rm = TRUE),
          mwc = tapply(r$type == "mwc", r$m, mean, na.rm = TRUE),
          surplus = tapply(r$type == "surplus", r$m, mean, na.rm = TRUE),
          slope = sl)
  }
  bb <- cut(b$runs$beta, c(0, .05, .1, .2, .4, .7, 1), include.lowest = TRUE)
  bo <- cut(o$runs$beta, c(0, .05, .1, .2, .4, .7, 1), include.lowest = TRUE)
  by_beta <- function(r, f) rbind(delay = tapply(r$delay, f, mean),
                                  minority = tapply(r$type %in% c("spmin", "mincoal"), f, mean),
                                  n_gov = tapply(r$n_gov, f, mean))
  list(headline = head_tab,
       by_m = list(ggs = round(by_m(o$runs, o$gov), digits),
                   rep = round(by_m(b$runs, b$gov), digits)),
       by_beta = list(ggs = round(by_beta(o$runs, bo), digits),
                      rep = round(by_beta(b$runs, bb), digits)))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  b <- readRDS(a[1])
  cmp <- compare_batch(b)
  options(width = 150)
  print(cmp$headline, digits = 3)
  print(cmp$by_m); print(cmp$by_beta)
  if (length(a) >= 2) saveRDS(cmp, a[2])
}
