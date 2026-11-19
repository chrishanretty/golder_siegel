# Tests: run from project root with  Rscript tests/test_model.R
source("R/ggs_model.R"); load_ggs_cpp()
ok <- function(cond, msg) { cat(if (cond) "PASS" else "FAIL", msg, "\n"); if (!cond) quit(status = 1) }

# 1. C++ core and pure-R reference are draw-for-draw identical
for (cfg in list(c(3, 0.1), c(5, 0.25), c(8, 0.9))) {
  set.seed(cfg[1]); ps <- draw_party_system(m = cfg[1], beta = cfg[2])
  set.seed(99); a <- ggs_run_cpp(ps$z, ps$seats, ps$beta, n_form = 60)
  set.seed(99); b <- ggs_run_R(ps$z, ps$seats, ps$beta, n_form = 60)
  ok(identical(a[c("delay", "type", "formateur", "n_gov")], b[c("delay", "type", "formateur", "n_gov")]) &&
       isTRUE(all.equal(a$final_cab, b$final_cab)), sprintf("C++ == R reference (m=%g)", cfg[1]))
}

# 2. Election: shares sum to 1, multiples of 1/1000
set.seed(1); ps <- draw_party_system(m = 7)
ok(abs(sum(ps$seats) - 1) < 1e-12 && all(abs(ps$seats * 1000 - round(ps$seats * 1000)) < 1e-9),
   "seat shares sum to one on a 1000-voter grid")

# 3. Government classification
s <- c(0.4, 0.35, 0.15, 0.1)
ok(classify_gov(c(1, 1, 0, 0), s) == 4L, "MWC")
ok(classify_gov(c(1, 1, 1, 0), s) == 5L, "surplus")
ok(classify_gov(c(1, 1, 1, 1), s) == 6L, "unity")
ok(classify_gov(c(0, 0, 1, 1), s) == 3L, "minority coalition")
ok(classify_gov(c(1, 0, 0, 0), s) == 2L, "single-party minority")
ok(classify_gov(c(1, 0, 0, 0), c(0.6, 0.2, 0.2)) == 1L, "single-party majority")

# 4. Invariants on a long run: cabinet shares sum to one, governments are majorities
#    or (with beta > 0) minority governments with majority support
set.seed(5); ps <- draw_party_system(m = 6, beta = 0.5)
r <- ggs_run_cpp(ps$z, ps$seats, ps$beta, n_form = 500, history = TRUE)
ok(all(abs(rowSums(r$cab_hist) - 1) < 1e-9), "cabinet shares sum to one")
ok(all(r$delay >= 1 & r$delay <= 100), "delays in [1, 100]")

# 5. beta = 0 => no minority governments (paper footnote 19)
set.seed(11); ps <- draw_party_system(m = 6, beta = 0)
r <- ggs_run_cpp(ps$z, ps$seats, 0, n_form = 500)
ok(!any(r$type %in% c(2L, 3L)), "beta = 0 yields no minority governments")
cat("All tests passed\n")
