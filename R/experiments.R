# Quick experiments on implementation variants (cabinet-share draw, policy range,
# voting metric). Each variant is run on the same set of party systems; the
# stationary distribution is estimated from formations 201..n_form of each run.
source("R/ggs_model.R")
suppressPackageStartupMessages(library(parallel))

el_metric <- function(z, V, metric = "L2") {
  dx <- abs(outer(V[, 1], z[, 1], "-")); dy <- abs(outer(V[, 2], z[, 2], "-"))
  D <- switch(metric, L2 = dx^2 + dy^2, Linf = pmax(dx, dy), L1 = dx + dy)
  tabulate(max.col(-D, "random"), nrow(z)) / nrow(V)
}

run_variant <- function(systems, cab_method = 0, policy_hw = 3, metric = "L2",
                        n_form = 700, burn = 200, cores = 8, seed = 1) {
  V <- voter_grid()
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  res <- mclapply(seq_along(systems), function(k) {
    ps <- systems[[k]]
    s <- el_metric(ps$z, V, metric)
    r <- ggs_run_cpp(ps$z, s, ps$beta, n_form = n_form, history = TRUE,
                     cab_method = cab_method, policy_half_width = policy_hw)
    keep <- (burn + 1):n_form
    samp <- seq(burn + 50, n_form, by = 50)     # government snapshots
    gl <- do.call(rbind, lapply(samp, function(f) {
      cs <- r$cab_hist[f, ]; g <- which(cs > 0)
      data.frame(sys = k, f = f, sc = s[g] / sum(s[g]), cs = cs[g],
                 form = as.integer(g == r$formateur[f]))
    }))
    list(sys = data.frame(sys = k, m = ps$m, beta = ps$beta,
                          delay = mean(r$delay[keep]), p100 = mean(r$delay[keep] == 100),
                          spmaj = mean(r$type[keep] == 1), spmin = mean(r$type[keep] == 2),
                          mincoal = mean(r$type[keep] == 3), mwc = mean(r$type[keep] == 4),
                          surplus = mean(r$type[keep] >= 5), n_gov = mean(r$n_gov[keep])),
         gov = gl)
  }, mc.cores = cores, mc.set.seed = TRUE)
  sy <- do.call(rbind, lapply(res, `[[`, "sys")); gv <- do.call(rbind, lapply(res, `[[`, "gov"))
  gv$m <- sy$m[gv$sys]
  bym <- rbind(sapply(split(sy[, -(1:3)], sy$m), colMeans),
               slope = sapply(split(gv, gv$m), function(x) coef(lm(cs ~ sc, x))[2]),
               form = sapply(split(gv, gv$m), function(x) coef(lm(cs ~ sc + form, x))[3]))
  list(overall = c(colMeans(sy[, -(1:3)]), slope = unname(coef(lm(cs ~ sc, gv))[2]),
                   int = unname(coef(lm(cs ~ sc, gv))[1]),
                   form = unname(coef(lm(cs ~ sc + form, gv))[3])),
       by_m = round(bym, 3))
}

make_systems <- function(n_per_m = 150, seed = 7) {
  set.seed(seed)
  unlist(lapply(3:10, function(m) replicate(n_per_m, {
    z <- matrix(runif(2 * m, -3, 3), m, 2); list(m = m, beta = runif(1), z = z)
  }, simplify = FALSE)), recursive = FALSE)
}
