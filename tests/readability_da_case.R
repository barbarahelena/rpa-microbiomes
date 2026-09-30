## Exercise the original and split reporting boundaries with fixed model outputs.
## Model fitting is stubbed deliberately: these checks test wiring and RNG order,
## not scientific inference. Every output stays inside the test sandbox.
args <- commandArgs(TRUE)
version <- args[1]
level <- args[2]
scenario <- args[3]
load_before_sites <- function(path) {
    for (e in parse(path)) {
        if (is.call(e) && identical(e[[1]], as.name("<-")) &&
            identical(e[[2]], as.name("sites"))) break
        eval(e, envir = .GlobalEnv)
    }
}
if (version == "before") {
    path <- if (level == "asv") "scripts/9_differential_abundance_16s.R" else
                                "scripts/9b_differential_abundance_genus_16s.R"
} else {
    path <- if (level == "asv") "scripts/12_differential_abundance_16s_asv_ethnicity_compute.R" else
                                "scripts/13_differential_abundance_16s_genus_ethnicity_compute.R"
}
load_before_sites(path)
Maaslin2 <- function(input_data, input_metadata, output, fixed_effects, reference, ...) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(counts = input_data, metadata = input_metadata,
                 fixed_effects = fixed_effects, reference = reference, settings = list(...)),
            file.path(output, "test_model_inputs.rds"))
    n <- ncol(input_data)
    q <- rep(0.8, n)
    if (scenario == "single") q[1] <- 0.01
    if (scenario == "multiple") q[1:3] <- c(0.001, 0.01, 0.1)
    results <- data.frame(feature = colnames(input_data), metadata = "EthnicityTotal",
                          value = "Moroccan", coef = seq(-1, 1, length.out = n),
                          stderr = 0.1, pval = q / 2, qval = q)
    readr::write_tsv(results, file.path(output, "all_results.tsv"))
    invisible(NULL)
}
ps <- readRDS("data/processed/ps_throat_rarefied.RDS")
ps <- prune_samples(sample_data(ps)$EthnicityTotal %in% c("Turkish", "Moroccan"), ps)
meta <- as(sample_data(ps), "data.frame")
meta$SeqBatch <- factor(rep(c("run1", "run2"), length.out = nrow(meta)))
sample_data(ps) <- sample_data(meta)
set.seed(20260928)
result <- if (level == "asv") {
    run_da_pair(ps, "throat", "Turkish", "Moroccan", character(), "SeqBatch", outdir)
} else {
    run_da_pair(ps, "throat", "Turkish", "Moroccan", character(), outdir)
}
saveRDS(.Random.seed, "compute_rng.rds")
if (version == "before") {
    saveRDS(result, "summary.rds")
} else {
    saveRDS(result$summary, "summary.rds")
    cache <- list(pairs = list(result$report), confounders = list(), summary = result$summary)
    if (level == "genus") cache$agglomeration_summary <- tibble()
    dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(cache, file.path(outdir, "cache", paste0("da_", level, ".rds")))
    report <- sub("_compute.R", "_report.R", path, fixed = TRUE)
    source(report, print.eval = TRUE)
}
