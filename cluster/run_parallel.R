# Parallel, reproducible replacement for generate_data() in gen_data_proj_cp.R.
#
# The original generate_data() runs `num_runs` independent replications in a
# serial for-loop, writing to data_path/<n_i>/<timestamp>.csv after each run.
# This reimplements the same per-run computation (same helper functions, same
# super_eigs cache reset per run) but runs replications in parallel chunks of
# `mc.cores` via parallel::mclapply, and only writes to disk once per chunk,
# in the parent process, to avoid concurrent writes from forked workers.
#
# Reproducibility: uses RNGkind("L'Ecuyer-CMRG") + a single set.seed() call
# before any mclapply, with mc.set.seed = TRUE. Each replication gets its own
# RNG substream keyed to its run index, independent of mc.cores or chunking,
# so re-running this script with the same parameters and seed reproduces the
# same numbers.
#
# Env vars OPENBLAS_NUM_THREADS=1 / OMP_NUM_THREADS=1 must be set in the shell
# that launches this script, so each forked worker's BLAS calls stay single
# threaded (otherwise mc.cores processes x multi-threaded BLAS oversubscribes
# the machine). data.table is also pinned to 1 thread for the same reason.

suppressMessages(library(parallel))
suppressMessages(library(data.table))
setDTthreads(1)

generate_data_parallel <- function(num_runs, mc.cores, seed, chunk_size = mc.cores){
    n_max <- max(ns)
    dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
    for (n_i in ns) dir.create(file.path(data_path, n_i), recursive = TRUE, showWarnings = FALSE)

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    run_one <- function(run){
        super_eigs <<- list()   # fresh cache per run, same as the serial version

        if (type == "bm"){
            sample_max <- generate_brownian_motion_sample(n = n_max, S = S)
        } else if (type == "ns"){
            sample_max <- generate_nonseparable_sample(n = n_max, S = S)
        } else if (type == "ma"){
            sample_max <- generate_ma_sample(n = n_max, S = S)
        }

        out <- list()
        for (n_i in ns){
            samp <- get_lambda_subsample(n_i / n_max, samp = sample_max)
            out[[paste(n_i)]] <- get_export_matrix(run_num = run, samp = samp)
        }
        out
    }

    n_chunks <- ceiling(num_runs / chunk_size)
    all_results <- vector("list", num_runs)
    t_start <- Sys.time()

    for (chunk in 1:n_chunks){
        idx <- ((chunk - 1) * chunk_size + 1):min(chunk * chunk_size, num_runs)
        t0 <- Sys.time()

        chunk_results <- mclapply(idx, run_one, mc.cores = mc.cores,
                                   mc.preschedule = FALSE, mc.set.seed = TRUE)

        errs <- which(vapply(chunk_results, function(x) inherits(x, "try-error"), logical(1)))
        if (length(errs) > 0){
            stop(paste("run(s) failed:", paste(idx[errs], collapse = ", "),
                       "\nfirst error:", attr(chunk_results[[errs[1]]], "condition")$message))
        }

        all_results[idx] <- chunk_results

        elapsed <- as.numeric(Sys.time() - t0, units = "secs")
        total_elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
        eta <- (total_elapsed / chunk) * (n_chunks - chunk)
        cat(sprintf("[%s] chunk %d/%d (runs %d-%d) done in %.1fs | total %.1fs | ETA %.1fs\n",
                    format(Sys.time(), "%H:%M:%S"), chunk, n_chunks, min(idx), max(idx),
                    elapsed, total_elapsed, eta))
    }

    for (n_i in ns){
        dt_list <- lapply(all_results, function(r) {
            dt <- data.table(r[[paste(n_i)]], keep.rownames = TRUE)
            colnames(dt) <- c("", colnames(dt)[2:length(colnames(dt))])
            dt
        })
        combined <- rbindlist(dt_list)
        f_name <- file.path(data_path, n_i,
                             paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"))
        fwrite(combined, f_name, row.names = FALSE, col.names = TRUE)
        cat("wrote", f_name, "(", nrow(combined), "rows )\n")
    }
}
