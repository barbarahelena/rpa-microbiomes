## Build small synthetic inputs for the structural-refactor checks (no cohort data).
## Generate metadata, ASV counts, taxonomy, and trees with a fixed random seed;
## save nose/throat phyloseq fixtures under data/processed in the test workspace.

library(phyloseq)
library(ape)
set.seed(101)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
## ---- Generate synthetic metadata, counts, taxonomy, and trees ----
for (site in c("throat", "nose")) {
    groups <- rep(c("Dutch", "Turkish", "Moroccan", "Ghanaian"), c(53, 51, 51, 12))
    n <- length(groups)
    first <- rep(c(TRUE, FALSE), length.out = n)
    duration <- runif(n, 1, 40)
    meta <- data.frame(
        EthnicityTotal = factor(groups),
        MigrationGen = factor(ifelse(groups == "Dutch", NA,
                                    ifelse(first, "1st generation", "2nd generation"))),
        ResidenceDuration_BA = ifelse(first, duration, NA),
        DifficultyDutch_BA = rnorm(n),
        CultFeelBerrys_BA = factor(rep(c("a", "b"), length.out = n)),
        CultOrientBerrys_BA = rnorm(n),
        CultNetworkBerrys_BA = rnorm(n),
        CultDistMeanScore0_BA = rnorm(n),
        DiscrMean_BA = rnorm(n),
        Antibiotics_FU = factor(rep("No", n)),
        row.names = paste0("sample", seq_len(n))
    )
    counts <- matrix(rpois(12 * n, 30), nrow = 12,
                     dimnames = list(paste0("ASV_", 1:12), rownames(meta)))
    counts[1, ] <- round(duration * 12 + 1)
    tree <- rtree(12, tip.label = rownames(counts))
    tax <- matrix(c(rep("Bacteria", 12), rep("Firmicutes", 12),
                    paste0("Genus", 1:12), paste0("Taxon", 1:12), rownames(counts)),
                  nrow = 12, dimnames = list(rownames(counts),
                      c("Kingdom", "Phylum", "Genus", "Tax", "ASV")))
    ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(meta),
                   tax_table(tax), phy_tree(tree))
    ## ---- Save the synthetic phyloseq fixture for this site ----
    saveRDS(ps, paste0("data/processed/ps_", site, "_rarefied.RDS"))
}
