## Refresh the small inputs that Pixi cannot infer from pixi.lock alone.
## Always run this task before cache lookup; unchanged values keep identical
## bytes. No participant data or analysis results are read here.
cache_dir <- "cache/analysis"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

write_if_changed <- function(value, name) {
    path <- file.path(cache_dir, name)
    text <- capture.output(dput(value))
    if (!file.exists(path) || !identical(readLines(path, warn = FALSE), text)) {
        tmp <- tempfile(pattern = paste0(name, "-"), tmpdir = cache_dir)
        writeLines(text, tmp)
        if (!file.rename(tmp, path)) stop("Cannot update cache context: ", path)
    }
}

## Include packages installed by BiocManager, not only Pixi-managed packages.
## Use the active library order so changing which installation is selected
## also invalidates the cache. Scanning DESCRIPTION metadata does not load R
## packages (or open graphics devices).
packages <- installed.packages(noCache = TRUE)[, c("Package", "Version", "LibPath", "Built"), drop = FALSE]
packages <- packages[order(packages[, "Package"], packages[, "LibPath"]), , drop = FALSE]
write_if_changed(list(
    R = R.version.string,
    platform = R.version$platform,
    libraries = .libPaths(),
    packages = packages,
    locale = Sys.getlocale(),
    environment = Sys.getenv(c("TZ", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"))
), "runtime.R")
write_if_changed(Sys.getenv(c("BETA_DIV_TEST_N", "BETA_DIV_N_CORES")), "beta.R")
write_if_changed(Sys.getenv("DIFF_AB_TEST_N"), "diffab.R")

## Persist the refresh token: unsetting PIPELINE_FORCE after a forced run
## must not trigger another rebuild. One token is shared by the entire DAG.
refresh_path <- file.path(cache_dir, "refresh.R")
if (identical(Sys.getenv("PIPELINE_FORCE"), "1")) {
    write_if_changed(paste(Sys.time(), Sys.getpid(), tempfile()), "refresh.R")
} else if (!file.exists(refresh_path)) {
    write_if_changed("initial", "refresh.R")
}
