# Command-line runner for scripts/allInsects_fit.R
#
# Usage (from the BITS repository root):
#   Rscript scripts/run_allInsects_server.R <master_dir> [n_chains] [n_iter] [n_burnin] [n_thin] [parallel]
#
#   master_dir  folder holding 00_masterScript.R's project root
#               (the `kwaku_work_here` folder, which `here::here()` resolves)
#   defaults    n_chains = 4, n_iter = 6000, n_burnin = 3000, n_thin = 3,
#               parallel = TRUE
#
# Example:
#   Rscript scripts/run_allInsects_server.R /p/826001.../kwaku_work_here 4 6000 3000 3 TRUE
#
# Output (fit, PPC, residual and map figures) goes to scripts/allInsects_output/.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript run_allInsects_server.R <master_dir> [n_chains n_iter n_burnin n_thin parallel]")

repo <- normalizePath(getwd())
stopifnot(file.exists(file.path(repo, "DESCRIPTION")))

master_dir <- args[1]
opt <- function(i, default, f = as.numeric) if (length(args) >= i) f(args[i]) else default
n_chains_fit <- opt(2, 4)
n_iter_fit   <- opt(3, 6000)
n_burnin_fit <- opt(4, 3000)
n_thin_fit   <- opt(5, 3)
parallel_fit <- opt(6, TRUE, as.logical)

# Keep compiled models and BLAS single-threaded on a shared server
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# Data preparation: run the master script up to (not including) the model fit.
# Its server branch expects exactly two arguments (version, flag).
setwd(master_dir)
master <- readLines(file.path("bits", "modelFit", "00_masterScript.R"))
last <- grep("allInsectData.R", master)[1]
# Keep everything before the `source(here::here(file.path("bits", ... allInsectData.R` call
prep <- master[seq_len(last - 5)]
prep <- prep[!grepl("^\\s*setwd\\(", prep)]
prep <- sub("^USER <- .*", 'USER <- "server"', prep)
prep <- c(prep, "}")  # close the `if (taxa_in_focus ...)` block

commandArgs <- function(trailingOnly = TRUE) c("1", "TRUE")
eval(parse(text = prep), envir = globalenv())
rm(commandArgs)
stopifnot(exists("all_insect_data"))
message("Rows in all_insect_data: ", nrow(all_insect_data))

# Load BINTS from source so that local edits are used
setwd(repo)
devtools::load_all(repo, quiet = TRUE)

t_run <- system.time(source(file.path("scripts", "allInsects_fit.R")))
message("Finished in ", round(t_run[["elapsed"]] / 60, 1), " minutes. Output: ",
        file.path(repo, "scripts", "allInsects_output"))
