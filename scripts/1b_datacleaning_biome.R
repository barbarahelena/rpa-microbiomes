## Data cleaning of microbiome data

## Libraries
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(vegan)
library(decontam)
library(Biostrings)

## Functions
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
      }

# 16S dataset
meta <- rio::import("data/raw/sample_sheet_withmeta.csv")
names(meta)
meta$sample <- str_c("S", meta$sample)
ps <- readRDS("data/raw/phyloseq/complete/phyloseq.RDS")

# Clean sample names
sample_names(ps) <- str_remove(sample_names(ps), "_T1")
sample_names(ps) <- str_replace(sample_names(ps), "X", "S") # IDs are now S[0-9]* instead of X[0-9]*
ps # 2741 taxa and 4501 samples (throat and nose)
sample_sums(ps) # unrarefied
negctrlsvec <- sample_names(ps)[str_detect(sample_names(ps), "NEG")]
negctrls <- prune_samples(negctrlsvec, ps)

df <- data.frame(sample = sample_names(ps))
df$Ctrl <- case_when(df$sample %in% negctrlsvec ~ TRUE, .default = FALSE)
df$Count <- sample_sums(ps)
df <- left_join(df, meta)
rownames(df) <- df$sample
head(df)
sample_data(ps) <- sample_data(df)
ps

# Remove too long amplicons
seq_lengths <- width(refseq(ps))
summary(seq_lengths)
summary(seq_lengths > 265) # 14 taxa with amplicons > 265
seq_lengths[which(seq_lengths > 265)] # few from 274-286, rest > 350
pstoolong <- prune_taxa(seq_lengths > 265, ps)
taxa_sums(pstoolong) # most abundant ASV > 265 has ~800 counts
ps <- prune_taxa(seq_lengths <= 265, ps)

# Remove mitochondria
tax_table(ps) |> as.data.frame() |> filter(str_detect(Family, "Mitochondria|mitochondria"))
ps <- subset_taxa(ps, !str_detect(Family, "Mitochondria") | is.na(Family))

# Frequency filter decontam
ps_freq <- prune_samples(!is.na(sample_data(ps)$Nucl_Acid_Conc) & sample_data(ps)$Nucl_Acid_Conc > 0, ps)
contam_freq <- isContaminant(ps_freq, method = "frequency", conc = "Nucl_Acid_Conc")
cont1 <- rownames(contam_freq)[which(contam_freq$contaminant == TRUE)]
table(contam_freq$contaminant)
plot_frequency(ps_freq, taxa_names(ps_freq)[sample(which(contam_freq$contaminant),3)], conc="Nucl_Acid_Conc") +
    xlab("DNA Concentration")

# Prevalence filter compared to neg ctrl
contamdf_prev <- isContaminant(ps, method="prevalence", neg="Ctrl")
table(contamdf_prev$contaminant)
cont2 <- rownames(contamdf_prev)[which(contamdf_prev$contaminant == TRUE)]

ps_pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
psneg <- prune_samples(sample_data(ps_pa)$Ctrl == TRUE, ps_pa)
pssample <- prune_samples(sample_data(ps_pa)$Ctrl == FALSE, ps_pa)
dfpa <- data.frame(pssample=taxa_sums(pssample), psneg=taxa_sums(psneg),
                      contaminant=contamdf_prev$contaminant)
ggplot(data=dfpa, aes(x=psneg, y=pssample, color=contaminant)) + geom_point() +
  scale_color_manual(values = c("royalblue4", "firebrick2")) + 
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + theme_Publication()

# Remove contaminants
contams <- vctrs::vec_c(cont1, cont2)
length(contams) # remove 35 ASVs
ps_noncontam <- prune_taxa(!taxa_names(ps) %in% contams, ps)
ps_noncontam

# Remove neg controls
psnoneg <- prune_samples(!sample_names(ps_noncontam) %in% negctrlsvec, ps_noncontam)
psnoneg

# Remove duplicate samples, keeping the sample with the highest read count
dup_df <- tibble(sample = sample_names(psnoneg), reads = sample_sums(psnoneg)) |>
  mutate(base_id = str_remove(sample, "_2$"),
         is_dup = str_detect(sample, "_2$")) |>
  group_by(base_id) |>
  filter(n() > 1) |>
  mutate(keep = reads == max(reads),
         msg = paste0(sample, ": ", reads, " reads")) |>
  summarise(message = paste(msg, collapse = " vs "),
            to_remove = sample[!keep])

walk(dup_df$message, \(x) cat(x, "\n"))
psnoneg <- prune_samples(!sample_names(psnoneg) %in% dup_df$to_remove, psnoneg)
sample_names(psnoneg) <- str_remove(sample_names(psnoneg), "_2")

# Identify throat and nose samples
throat_samples <- sample_names(psnoneg)[str_detect(sample_names(psnoneg), "Throat")]
nose_samples <- sample_names(psnoneg)[str_detect(sample_names(psnoneg), "Nose")]

# Subset phyloseq objects
ps_throat <- prune_samples(throat_samples, psnoneg)
sample_names(ps_throat) <- str_remove(sample_names(ps_throat), "_Throat")
ps_throat

ps_nose <- prune_samples(nose_samples, psnoneg)
sample_names(ps_nose) <- str_remove(sample_names(ps_nose), "_Nose")
ps_nose

# Save the split phyloseq objects
saveRDS(ps_throat, "data/processed/ps_throat.RDS")
saveRDS(ps_nose, "data/processed/ps_nose.RDS")

# Rarefaction nose samples
otu_nose <- t(as(otu_table(ps_nose), "matrix"))
gghistogram(sample_sums(ps_nose), fill = "royalblue4") + scale_x_log10()
rarecurve(otu_nose, step = 100, sample = min(rowSums(otu_nose)), label = FALSE)

grDevices::pdf(NULL) # this is to prevent the next line from printing all plots
rare_data <- lapply(seq_len(nrow(otu_nose)), function(i) {
  vegan::rarecurve(
    otu_nose[i, , drop = FALSE],
    step   = 1000,
    sample = min(rowSums(otu_nose)),
    label  = FALSE,
    plot   = FALSE
  )
})
dev.off()
df <- bind_rows(lapply(1:length(rare_data), function(i){
  tibble(
    sample = rownames(otu)[i],
    reads = as.numeric(str_remove(names(rare_data[[i]][[1]]), "N")),
    richness = as.numeric(rare_data[[i]][[1]])
  )
}))

df_summary <- df %>%
  filter(reads %in% c(0, seq(1, 40001, by = 1000))) %>%
  group_by(reads) %>%
  summarise(
    mean_richness = mean(richness),
    lowerse = mean(richness) - sd(richness) / sqrt(n()),
    higher = mean(richness) + sd(richness) / sqrt(n())
  )

ggplot(df, aes(x = reads, y = richness)) +
  geom_line(aes(group = sample), alpha = 0.3, color = "steelblue") +
  geom_line(data = df_summary, aes(x = reads, y = mean_richness), color = "darkblue", size = 1) +
  theme_Publication() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 40000)) + 
  labs(x = "Sequencing depth", y = "ASV richness") + 
  geom_vline(xintercept = 5000, linetype = "dashed", color = "black")
ggsave("results/rarefaction/arefactioncurve_nose_long.pdf", width = 6, height = 4)

ggplot(df, aes(x = reads, y = richness)) +
  geom_line(aes(group = sample), alpha = 0.3, color = "steelblue") +
  geom_line(data = df_summary, aes(x = reads, y = mean_richness), color = "darkblue", size = 1) +
  theme_Publication() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 10000)) + 
  labs(x = "Sequencing depth", y = "ASV richness") + 
  geom_vline(xintercept = 2500, linetype = "dashed", color = "black")
ggsave("results/rarefaction/rarefactioncurve_nose_short.pdf", width = 6, height = 4)


# Rarefaction throat samples
otu_throat <- t(as(otu_table(ps_throat), "matrix"))
gghistogram(sample_sums(ps_throat), fill = "royalblue4") + scale_x_log10()
rarecurve(otu_throat, step = 1000, sample = min(rowSums(otu_throat)), label = FALSE)

grDevices::pdf(NULL) # to prevent it from printing all plots
rare_data <- lapply(seq_len(nrow(otu_throat)), function(i) {
  vegan::rarecurve(
    otu_throat[i, , drop = FALSE],
    step   = 1000,
    sample = min(rowSums(otu_throat)),
    label  = FALSE,
    plot   = FALSE # this doesn't seem to help
  )
})
dev.off()
# Convert rarecurve output into a tidy data frame
df <- bind_rows(lapply(1:length(rare_data), function(i){
  tibble(
    sample = rownames(otu)[i],
    reads = as.numeric(str_remove(names(rare_data[[i]][[1]]), "N")),
    richness = as.numeric(rare_data[[i]][[1]])
  )
}))

df_summary <- df %>%
  filter(reads %in% c(0, seq(1, 50001, by = 1000))) %>%
  group_by(reads) %>%
  summarise(
    mean_richness = mean(richness),
    lowerse = mean(richness) - sd(richness) / sqrt(n()),
    higher = mean(richness) + sd(richness) / sqrt(n())
  )

ggplot(df, aes(x = reads, y = richness)) +
  geom_line(aes(group = sample), alpha = 0.3, color = "steelblue") +
  geom_line(data = df_summary, aes(x = reads, y = mean_richness), color = "darkblue", size = 1) +
  theme_Publication() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 50000)) + 
  labs(x = "Sequencing depth", y = "ASV richness") + 
  geom_vline(xintercept = 7500, linetype = "dashed", color = "black")
ggsave("results/rarefaction/rarefactioncurve_throat_long.pdf", width = 6, height = 4)

ggplot(df, aes(x = reads, y = richness)) +
  geom_line(aes(group = sample), alpha = 0.3, color = "steelblue") +
  geom_line(data = df_summary, aes(x = reads, y = mean_richness), color = "darkblue", size = 1) +
  theme_Publication() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 20000)) + 
  labs(x = "Sequencing depth", y = "ASV richness") + 
  geom_vline(xintercept = 7500, linetype = "dashed", color = "black")
ggsave("results/rarefaction/rarefactioncurve_throat_short.pdf", width = 6, height = 4)

# Rarefy phyloseq objects
set.seed(123)
ps_nose_rarefied <- rarefy_even_depth(ps_nose, sample.size = 2000, rngseed = 123) # 612 samples removed
ps_throat_rarefied <- rarefy_even_depth(ps_throat, sample.size = 7500, rngseed = 123) # 70 samples removed

ps_nose_rarefied # 2243 taxa and 1242 samples
ps_throat_rarefied # 1626 taxa and 2390 samples
sample_names(ps_nose_rarefied)
sample_names(ps_throat_rarefied)

# Save rarefied phyloseq objects
saveRDS(ps_nose_rarefied, "data/processed/ps_nose_rarefied.RDS")
saveRDS(ps_throat_rarefied, "data/processed/ps_throat_rarefied.RDS")

# Metagenomics: merge batches
batch1 <- rio::import("data/raw/combined_table.txt")
dim(batch1) # 3 samples in the last batch (broken fastqs)
batch2 <- rio::import("data/raw/combined_table2.txt")
dim(batch2) # 317 samples in the first batch
tot <- full_join(batch1, batch2)
tot <- tot %>% filter(str_detect(clade_name, "s__") & !str_detect(clade_name, "t__")) # select species

# Clean clade names
clade <- tot$clade_name
cladesplit <- str_split(clade, "\\|", n = 8, simplify = TRUE)
cladesplit <- as.data.frame(cladesplit[,-8])
colnames(cladesplit) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cladesplit <- cladesplit %>% mutate(across(everything(.), ~str_remove_all(.x, "[a-z]__")))
cladesplit$rowname <- clade
saveRDS(cladesplit, "data/processed/shotgun_taxtable.RDS")
rownames(tot) <- cladesplit$Species[match(cladesplit$rowname, tot$clade_name)]
tot$clade_name <- NULL
tot <- t(as.matrix(tot))
head(tot)[1:5,1:5]

# Replace IDs with helius IDs
meta <- readRDS("data/processed/HELIUSmetadata_clean.RDS")
rownames(tot) 
meta$TongueSampleID
sampleID <- case_when(rownames(tot) %in% meta$TongueSampleID ~ str_c(str_extract(rownames(tot), "[0-9]*(_)"), "Tongue"),
                      rownames(tot) %in% meta$ThroatSampleID ~ str_c(str_extract(rownames(tot), "[0-9]*(_)"), "Throat"))
all(str_extract(rownames(tot), "[0-9]*(_)") == str_extract(sampleID, "[0-9]*(_)")) # check if matches: TRUE
sampleID
rownames(tot) <- sampleID
tot <- as.data.frame(tot)
tot <- tot |> mutate(across(everything(), function(x) case_when(is.na(x) ~ 0, .default = x)))
rowSums(tot)
saveRDS(tot, "data/processed/metaphlantable_tot.RDS")

# Extract tongue data
tongue <- tot[ str_detect(rownames(tot), "Tongue"), ]
rownames(tongue) <- str_remove(rownames(tongue), "_Tongue")
head(tongue)[1:5,1:5]
saveRDS(tongue, "data/processed/metaphlan_tongue.RDS")

# Extract throat data
throat <- tot[ str_detect(rownames(tot), "Throat"), ]
rownames(throat) <- str_remove(rownames(throat), "_Throat")
head(throat)[1:5,1:5]
saveRDS(throat, "data/processed/metaphlan_throat.RDS")
