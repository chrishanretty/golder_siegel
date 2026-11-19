# Parallel harness: simulate many party systems (runs) of the GGS model.
#
# * Runs are split into chunks; chunk k uses the k-th L'Ecuyer-CMRG stream
#   derived from `seed`, so results are reproducible and independent of the
#   number of cores.
# * Each finished chunk is written to `out_dir`, so an interrupted batch
#   resumes where it stopped.
# * Forking (mclapply) on Unix; a PSOCK cluster elsewhere.

source(file.path(.ggs_dir_h <- {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(f)) "R" else dirname(f)
}, "ggs_model.R"))

chunk_seeds <- function(seed, n_chunks) {
  old <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  s <- vector("list", n_chunks); s[[1]] <- .Random.seed
  if (n_chunks > 1) for (k in 2:n_chunks) s[[k]] <- parallel::nextRNGStream(s[[k - 1]])
  RNGkind(old)
  s
}

run_chunk <- function(k, ids, seed_k, P, out_dir) {
  f <- file.path(out_dir, sprintf("chunk_%05d.rds", k))
  if (file.exists(f)) return(f)
  load_ggs_cpp()
  RNGkind("L'Ecuyer-CMRG"); assign(".Random.seed", seed_k, envir = globalenv())
  voters <- voter_grid(P$voters_nx, P$voters_ny, P$half_width)
  res <- lapply(ids, function(i) ggs_one_run(i, P = P, voters = voters))
  out <- list(runs = do.call(rbind, lapply(res, `[[`, "run")),
              gov = do.call(rbind, lapply(res, `[[`, "gov")),
              parties = do.call(rbind, lapply(res, `[[`, "parties")))
  saveRDS(out, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  f
}

run_batch <- function(n_runs = 50000, out_dir = "outputs/sim_main", seed = 20120401,
                      chunk_size = 250, n_cores = parallel::detectCores(),
                      P = ggs_defaults) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(n_runs = n_runs, seed = seed, chunk_size = chunk_size, P = P),
          file.path(out_dir, "config.rds"))
  load_ggs_cpp()                                   # compile once in the parent
  chunks <- split(seq_len(n_runs), ceiling(seq_len(n_runs) / chunk_size))
  seeds <- chunk_seeds(seed, length(chunks))
  todo <- which(!file.exists(file.path(out_dir, sprintf("chunk_%05d.rds", seq_along(chunks)))))
  message(sprintf("%d runs in %d chunks (%d to do) on %d cores",
                  n_runs, length(chunks), length(todo), n_cores))
  t0 <- Sys.time()
  if (length(todo)) {
    job <- function(k) run_chunk(k, chunks[[k]], seeds[[k]], P, out_dir)
    if (.Platform$OS.type == "unix" && n_cores > 1) {
      r <- parallel::mclapply(todo, job, mc.cores = n_cores, mc.preschedule = FALSE)
    } else {
      cl <- parallel::makeCluster(n_cores); on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c(".ggs_dir"), envir = globalenv())
      parallel::clusterEvalQ(cl, { source(file.path(.ggs_dir, "ggs_model.R")); load_ggs_cpp() })
      r <- parallel::parLapplyLB(cl, todo, job)
    }
    bad <- vapply(r, inherits, logical(1), "try-error")
    if (any(bad)) stop("Chunks failed: ", paste(todo[bad], collapse = ", "))
  }
  message("Elapsed: ", format(round(Sys.time() - t0, 1)))
  collect_batch(out_dir)
}

collect_batch <- function(out_dir) {
  fs <- sort(list.files(out_dir, "^chunk_\\d+\\.rds$", full.names = TRUE))
  l <- lapply(fs, readRDS)
  out <- list(runs = do.call(rbind, lapply(l, `[[`, "runs")),
              gov = do.call(rbind, lapply(l, `[[`, "gov")),
              parties = do.call(rbind, lapply(l, `[[`, "parties")))
  out$runs <- out$runs[order(out$runs$run), ]
  out$runs$type[out$runs$type == "unity"] <- "surplus"   # as GGS
  out$runs$type[out$runs$type == "none"] <- NA
  out
}

if (sys.nframe() == 0L) {                 # Rscript R/harness.R [n_runs] [out_dir] [cores]
  a <- commandArgs(TRUE)
  n <- if (length(a) >= 1) as.integer(a[1]) else 50000
  od <- if (length(a) >= 2) a[2] else "outputs/sim_main"
  nc <- if (length(a) >= 3) as.integer(a[3]) else parallel::detectCores()
  P <- ggs_defaults
  if (length(a) >= 4) P$metric <- a[4]            # optional voting metric
  b <- run_batch(n, od, n_cores = nc, P = P)
  saveRDS(b, paste0(od, ".rds"))
}
