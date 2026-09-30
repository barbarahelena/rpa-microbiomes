## Check that the refactor preserves R expressions against the Git baseline.
## Parse baseline/current scripts, ignoring comments and whitespace;
## compare cleaning code and shared/analysis functions, allowing the theme extraction.

## ---- Define parsers and extract functions for comparison ----
original <- function(name) parse(file.path("baseline", name))
current <- function(name) parse(file.path("scripts", name))
functions <- function(expressions) {
    result <- list()
    for (e in expressions) {
        if (is.call(e) && identical(e[[1]], as.name("<-")) &&
            is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
            result[[as.character(e[[2]])]] <- e[[3]]
    }
    result
}
## ---- Compare cleaning code, the extracted theme, and beta-diversity code ----
stopifnot(identical(original("1a_datacleaning_helius.R"), current("01_clean_metadata.R")))
a <- original("1b_datacleaning_biome.R")
b <- current("02_clean_microbiome.R")
is_theme <- function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
    identical(e[[2]], as.name("theme_Publication"))
is_source <- function(e) is.call(e) && identical(e[[1]], as.name("source"))
stopifnot(identical(a[!vapply(a, is_theme, logical(1))], b[!vapply(b, is_source, logical(1))]))
stopifnot(identical(functions(a)$theme_Publication, functions(current("lib/plot_style.R"))$theme_Publication))
stopifnot(identical(original("7a_beta_diversity_16s_compute.R"), current("10_beta_diversity_16s_ethnicity_compute.R")))
## ---- Compare analysis functions across the compute/report split ----
pairs <- list(
    c("6_alpha_diversity_16s_ethnicity.R", "08_alpha_diversity_16s_ethnicity_compute.R"),
    c("11_alpha_diversity_shotgun.R", "09_alpha_diversity_shotgun_ethnicity_compute.R"),
    c("8_beta_diversity_16s_migration.R", "11_beta_diversity_16s_migration_compute.R"),
    c("9_differential_abundance_16s.R", "12_differential_abundance_16s_asv_ethnicity_compute.R"),
    c("9b_differential_abundance_genus_16s.R", "13_differential_abundance_16s_genus_ethnicity_compute.R")
)
for (pair in pairs) {
    a <- functions(original(pair[1])); b <- functions(current(pair[2]))
    for (name in setdiff(intersect(names(a), names(b)), "run_da_pair"))
        stopifnot(identical(a[[name]], b[[name]]))
}
cat("Cleaning expressions, theme, ethnicity-beta computation and statistical helpers match.\n")
