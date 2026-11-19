# Golder, Golder & Siegel (2012) zero-intelligence government formation model.
# Single-run code: election stage, government-formation stage (C++ core and a
# pure-R reference implementation) and run-level summaries.

suppressPackageStartupMessages(library(Rcpp))

ggs_defaults <- list(
  half_width  = 3,        # policy space X = [-3,3]^2
  voters_nx   = 40,       # voter grid (nx * ny = 1000 voters, cf. GGS data:
  voters_ny   = 25,       #   all vote shares are multiples of 0.001)
  n_form      = 2500,     # government formations per run
  max_periods = 100,      # periods before the caretaker is reinstalled
  quota       = 0.5,      # majority: strictly more than half the seats
  metric      = "L2",     # voter distance: "L2" (paper) or "Linf" (sensitivity)
  m_range     = 3:10,
  beta_range  = c(0, 1)
)

.ggs_dir <- local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f)) "R" else dirname(f)
})

load_ggs_cpp <- function(dir = .ggs_dir) {
  if (!exists("ggs_run_cpp", mode = "function"))
    sourceCpp(file.path(dir, "model.cpp"), env = globalenv(),
              cacheDir = file.path(dir, ".rcpp_cache"))
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## Election stage
## ---------------------------------------------------------------------------
voter_grid <- function(nx = ggs_defaults$voters_nx, ny = ggs_defaults$voters_ny,
                       hw = ggs_defaults$half_width) {
  gx <- -hw + (2 * hw) * (seq_len(nx) - 0.5) / nx   # cell-centred, evenly spaced
  gy <- -hw + (2 * hw) * (seq_len(ny) - 0.5) / ny
  as.matrix(expand.grid(x = gx, y = gy))
}

# Proximity voting, ties broken by a fair die; seat share = vote share.
election <- function(z, voters = voter_grid(), metric = "L2") {
  dx <- abs(outer(voters[, 1], z[, 1], "-")); dy <- abs(outer(voters[, 2], z[, 2], "-"))
  d2 <- switch(metric, L2 = dx^2 + dy^2, Linf = pmax(dx, dy), L1 = dx + dy)
  dmin <- apply(d2, 1, min)
  choice <- vapply(seq_len(nrow(d2)), function(v) {
    w <- which(d2[v, ] <= dmin[v] + 1e-12)
    if (length(w) == 1) w else w[sample.int(length(w), 1)]
  }, integer(1))
  tabulate(choice, nbins = nrow(z)) / nrow(voters)
}

enp <- function(s) 1 / sum(s^2)

# Inverse of Alvarez & Nagler (2004) weighted compactness:
#   sum_j s_j |z_j - zbar| / sigma_c, zbar = seat-weighted mean party position,
#   |.| Euclidean distance, sigma_c = sd of voter ideal points (see summary doc).
diversity <- function(z, s, voters = voter_grid(), sigma_c = NULL) {
  # sigma_c = 2.5 is calibrated: the ratio of our raw (unscaled) measure to
  # GGS's published diversity variable is 2.50 for every m (cf. sqrt(6) = 2.45,
  # the 2-D standard deviation of a uniform voter density on [-3,3]^2).
  if (is.null(sigma_c)) sigma_c <- 2.5
  zbar <- colSums(z * s)
  sum(s * sqrt(rowSums(sweep(z, 2, zbar)^2))) / sigma_c
}

draw_party_system <- function(m = NULL, beta = NULL, P = ggs_defaults,
                              voters = voter_grid(P$voters_nx, P$voters_ny, P$half_width)) {
  if (is.null(m)) m <- P$m_range[sample.int(length(P$m_range), 1)]
  if (is.null(beta)) beta <- runif(1, P$beta_range[1], P$beta_range[2])
  z <- matrix(runif(2 * m, -P$half_width, P$half_width), m, 2)
  s <- election(z, voters, if (is.null(P$metric)) "L2" else P$metric)
  list(m = m, beta = beta, z = z, seats = s, enpp = enp(s),
       diversity = diversity(z, s, voters))
}

## ---------------------------------------------------------------------------
## Government-formation stage: pure-R reference implementation.
## Mirrors model.cpp draw-for-draw (used for testing; ~100x slower).
## ---------------------------------------------------------------------------
classify_gov <- function(in_gov, seats, quota = 0.5) {
  k <- sum(in_gov); tot <- sum(seats[in_gov == 1])
  if (k == 0) return(0L)
  maj <- tot > quota
  if (k == 1) return(if (maj) 1L else 2L)
  if (!maj) return(3L)
  if (any(tot - seats[in_gov == 1] > quota)) return(if (k == length(seats)) 6L else 5L)
  4L
}

ggs_run_R <- function(z, seats, beta, n_form = 2500, max_periods = 100,
                      half_width = 3, quota = 0.5, sq_policy = c(0, 0)) {
  m <- length(seats)
  sq_cab <- numeric(m); sq_in <- integer(m); sq_x <- sq_policy; sq_form <- 0L
  delay <- type <- formateur <- n_gov <- integer(n_form)
  for (f in seq_len(n_form)) {
    u_sq <- sq_cab - beta * ((sq_x[1] - z[, 1])^2 + (sq_x[2] - z[, 2])^2)
    chosen <- 0L
    for (t in seq_len(max_periods)) {
      mem <- matrix(0L, m, m); cab <- matrix(0, m, m); pol <- matrix(0, m, 2)
      for (i in seq_len(m)) {
        for (j in seq_len(m)) mem[i, j] <- if (j == i) 1L else as.integer(runif(1) < 0.5)
        for (j in which(mem[i, ] == 1L)) cab[i, j] <- -log(runif(1))
        cab[i, ] <- cab[i, ] / sum(cab[i, ])
        pol[i, ] <- -half_width + 2 * half_width * runif(2)
      }
      qual <- integer(0)
      for (i in seq_len(m)) {
        u <- cab[i, ] - beta * ((pol[i, 1] - z[, 1])^2 + (pol[i, 2] - z[, 2])^2)
        pref <- u > u_sq
        if (all(pref[mem[i, ] == 1L]) && sum(seats[pref]) > quota) qual <- c(qual, i)
      }
      if (length(qual)) {
        chosen <- if (length(qual) == 1) qual else qual[as.integer(runif(1) * length(qual)) + 1]
        break
      }
    }
    if (chosen > 0) {
      sq_cab <- cab[chosen, ]; sq_in <- mem[chosen, ]; sq_x <- pol[chosen, ]; sq_form <- chosen
    }
    delay[f] <- t; formateur[f] <- sq_form
    type[f] <- classify_gov(sq_in, seats, quota); n_gov[f] <- sum(sq_in)
  }
  list(delay = delay, type = type, formateur = formateur, n_gov = n_gov,
       final_cab = sq_cab, final_in = sq_in, final_policy = sq_x)
}

## ---------------------------------------------------------------------------
## One complete run: draw a party system, then 2,500 formations.
## Returns the outcomes at the final (2,500th) formation, as in GGS.
## ---------------------------------------------------------------------------
type_codes <- c("none", "spmaj", "spmin", "mincoal", "mwc", "surplus", "unity")

ggs_one_run <- function(run_id = 1L, P = ggs_defaults, ps = NULL, engine = c("cpp", "R"),
                        voters = voter_grid(P$voters_nx, P$voters_ny, P$half_width)) {
  engine <- match.arg(engine)
  if (is.null(ps)) ps <- draw_party_system(P = P, voters = voters)
  f <- if (engine == "cpp") ggs_run_cpp else ggs_run_R
  r <- f(ps$z, ps$seats, ps$beta, n_form = P$n_form, max_periods = P$max_periods,
         half_width = P$half_width, quota = P$quota)
  nf <- P$n_form
  gi <- which(r$final_in == 1L)
  tc <- tabulate(r$type + 1L, nbins = 7) / nf          # all-formation shares
  run <- data.frame(run = run_id, m = ps$m, beta = ps$beta, enpp = ps$enpp,
                    diversity = ps$diversity, delay = r$delay[nf],
                    formateur = r$formateur[nf],
                    type = type_codes[r$type[nf] + 1L], n_gov = length(gi),
                    mean_delay_all = mean(r$delay),
                    p100_all = mean(r$delay == max(P$max_periods)),
                    cum_spmaj = tc[2], cum_spmin = tc[3], cum_mincoal = tc[4],
                    cum_mwc = tc[5], cum_surplus = tc[6] + tc[7],
                    policy_x = r$final_policy[1], policy_y = r$final_policy[2])
  gov <- data.frame(run = rep(run_id, length(gi)), party = gi,
                    seat_contrib = ps$seats[gi] / sum(ps$seats[gi]),
                    cab_share = r$final_cab[gi],
                    formateur = as.integer(gi == r$formateur[nf]))
  parties <- data.frame(run = run_id, party = seq_len(ps$m), x = ps$z[, 1], y = ps$z[, 2],
                        seat = ps$seats)
  list(run = run, gov = gov, parties = parties)
}
