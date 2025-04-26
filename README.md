# HiMEx    High-Throughput Microvolume Extraction Method


**From Microscale to Microbial Insights: Validating High-Throughput Microvolume Extraction Methods (HiMEx) for Marine Microbial Ecology**

## Abstract

Extracting and amplifying DNA from small-volume, low-biomass samples enables rapid, ultra high-throughput analyses, and facilitates studying microbial communities where large-volume sample collection is challenging. This method can aid where ‘conventional’ filtration-based methods miss capturing smaller microbes, or where microscale variability matters, such as the ocean. Here, we develop and validate physical and chemical-based DNA extractions from microvolumes with universal rRNA gene amplicons and metagenomic sequencing for all domains and viruses, on natural seawater and experimentally manipulated marine waters. Compared to 500-mL filter-based extraction, direct PCR of 3 μL of lysed microvolume samples consistently captured comparable microbial community composition and diversity, with reliable amplification and little to no contamination. Metagenomic results of 10 μL of lysed microvolume samples captured high-quality assemblies, bacterial genomes, and diversity on par or better than our conventional extraction (48 vs 44 representative MAGs), and substantially more putative complete circular viral genomes (38 vs 1). Our approach enables processing hundreds of samples in per day for a fraction of the cost of conventional methods and complements existing microvolume approaches by increasing throughput and removing unnecessary expenses, like excess plasticware and expensive bead clean-up, expanding opportunities for more comprehensive microbial community monitoring, controlled laboratory experiments and field-based investigations.


## Lab Protocol

High-throughput Microvolume Extraction (HiMEx) methods were tested using:
- Physical extraction through thermal shock (freeze-thaw, FT),
- Combination of physical and chemical extraction (FT + proteinase K, FTP),
- Further enhanced disruption with surfactant (FTP + IGEPAL, FTPIG).


![Protocol](https://github.com/user-attachments/assets/283d3e5d-f8d2-402f-ba16-e34de1b96675)


## Citation
Ghotbi M, et al. (2024). From Microscale to Microbial Insights: Validating High-Throughput Microvolume Extraction Methods (HiMEx) for Marine Microbial Ecology. bioRxiv.
🔗 [https://doi.org/10.1101/2024.12.02.626238](https://doi.org/10.1101/2024.12.02.626238)



# Gene Marker Analysis Pipeline

**Author**: Marjan Ghotbi  
**Date**: April 23, 2025  

---



# Prokaryotes Amplicon Analysis 
Performs:
Metadata cleanup and matching
Rarefaction
Beta diversity (ordination, betadisper)
Differential abundance via DESeq2
Alpha diversity (Shannon, modeling with lme4)


# Load Packages

```{r}

# Load libraries
# install.packages("devtools")

devtools::install_github("mghotbi/DspikeIn", build_vignettes = TRUE, dependencies = TRUE)
browseVignettes("DspikeIn")



library(DspikeIn)
library(ggplot2)
library(phyloseq)
library(decontam)
library(ggrepel)
library(vegan)
library(DESeq2)
library(lme4)
library(lmerTest)
library(emmeans)
library(microbiome)
library(dplyr)
library(tidyverse)
library(extrafont)

```


# Load + Prep data

```{r}

# ───────────────────────────────────────────────────────────────
# Step 1: Load phyloseq objs
# ───────────────────────────────────────────────────────────────
library(here)

physeq_16S <- readRDS(here::here("data", "physeq16S.rds"))
# chloroplast <- readRDS("data/Chloroplast.rds")
# physeq_18S <- readRDS("data/physeq18S.rds")



```

# QC & Contaminant Detection

```{r}
# Pull metadata from phyloseq object
df <- microbiome::meta(physeq_16S)
df$LibrarySize <- sample_sums(physeq_16S)

# Order by LibrarySize and assign an index
df <- df[order(df$LibrarySize), ]
df$Index <- seq_len(nrow(df))
```

# Plot: Library Size per Sample (Colored by Method)
step 1: decontamination

```{r}
# Custom color palette for Methods
method_palette <- c(
  "conventional" = "#3D52A4",
  "FTPIG" = "#A6D38D",
  "FTP" = "#52B848",
  "FT" = "#8CCCF0",
  "Spike" = "#5C9DD5",
  "blank" = "steelblue2",
  "neg" = "red"
)

ggplot(df, aes(x = Index, y = LibrarySize, color = Method)) +
  geom_point(size = 4) +
  scale_color_manual(values = method_palette) +
  theme_minimal() +
  labs(
    title = "Library Size per Sample",
    x = "Sample Index",
    y = "Library Size",
    color = "Method"
  )
```

# Identify Blanks and Check Their Library Size

```{r}
# Subset to blanks
blank <- subset_samples(physeq_16S, Sample_blank == "Blank")

# Recalculate and visualize
df2 <- microbiome::meta(blank)
df2$LibrarySize <- sample_sums(blank)
df2$SampleID <- rownames(df2)
df2 <- df2[order(df2$LibrarySize), ]
df2$Index <- seq_len(nrow(df2))

# Plot library sizes for blanks
ggplot(df2, aes(x = Index, y = LibrarySize, color = Sample_blank)) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleID), hjust = 1.1, vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "Library Size by Sample (Blanks)",
    x = "Index",
    y = "Library Size"
  )

# Export for manual inspection if needed
write.csv(df2, "Method_lib_size_standard.csv")
```

# Identify Contaminants with Decontam (Prevalence Method)

```{r}
# Add logical vector for blanks to sample_data
sample_data(physeq_16S)$is.neg <- sample_data(physeq_16S)$Sample_blank == "Blank"

# Run decontam with default threshold
contamdf.prev <- isContaminant(physeq_16S, method = "prevalence", neg = "is.neg")
table(contamdf.prev$contaminant)

# Stricter threshold
contamdf.prev05 <- isContaminant(physeq_16S, method = "prevalence", neg = "is.neg", threshold = 0.05)
table(contamdf.prev05$contaminant)

# Extract contaminant ASV IDs
contaminants <- rownames(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
contaminants_prevalence <- contamdf.prev05[contamdf.prev05$contaminant == TRUE, ]

# Save list of contaminants
write.csv(contaminants_prevalence, "contaminants_prevalence_0.05_rarefaction4000.csv")
```
# Curated ASVs to Remove (Manual + Decontam + Spiked bactertia removal)
# Spiked Bacteria Codes:  
"1ae1bd31e0995f907062cc851d803bd7","39e27cbe03e2b71dbc8a5a1c8b8e7171","ac3a3e27aa0a24e16938ecd9de1c1060","1b468de6094a372d7fbc7d58b16f81c6"


```{r}

ASVs_to_remove <- c(
  # decontam + manual additions
  "ac164649e33237f799c1ffee4f344d5d",
  "92fb114641e27921122496428e2ef0dd",
  "3a8a7409bdc33da14a01ec3e3313b9a7",
  "3fcfe36ad281b523c2a5b9793512f5ad",
  "863728e1cec6befd5ba02d15baef4c36",
  "7c081f19c943c8d11819cdce94984719",
  "1c16af05dfb4ec260717d24d0b9b8274",
  "9076f0dc1a57c84a55fd88e445511903",
  "1af47a3186d12f284691296661bc7310",
  "f2d099ef5556f02539ded15181a3d994",
  "fd55e4491e6b3cf2259827d9f78367d6",
  "01cc2516697cf1b59d3bdcdc881c233e",
  "637ddfd0682928e2e0c4e0001883425e",
  "c728db0a39ffcd40e58dd84831880fde",
  "e1aee885aa820fc1bbc8eea6c27cdc3d",
  "ac70ecd361a0c67c771380f7aac46cc7",
  "79f37fee0660e917bd1debe546718bad",
  "4db5ac9fbbb89c5fc979475084c8c596",
  "645dd8b575c7e7c0952933bbe90b4bb9",
  "dd4ac7657b2665699cffd5a9b939a260",
  "ac9b52e609fb05e7ede04ae1e6a6ff2f",
  "e157cf6c27f73286aa20b0107e08e7d1",
  "1ae1bd31e0995f907062cc851d803bd7",
  "39e27cbe03e2b71dbc8a5a1c8b8e7171",
  "ac3a3e27aa0a24e16938ecd9de1c1060",
  "1b468de6094a372d7fbc7d58b16f81c6"
)





physeq_16S_clean <- prune_taxa(!taxa_names(physeq_16S) %in% ASVs_to_remove, physeq_16S)

physeq_16S_clean

```
# Step 2: Filter Non-informative Samples

```{r}
# Start with samples that have at least one read
physeq2 <- prune_samples(sample_sums(physeq_16S_clean) > 0, physeq_16S_clean)

# Remove blanks
physeq3 <- subset_samples(physeq2, Type != "Blank")

# Remove synthetic spike-ins
physeq4 <-subset_samples(physeq3, Type != "Spike")
physeq4


# Final object for rarefaction
physeq4
meta_raw<-sample_data(physeq4)
meta_raw
```

# Compute the Spearman Correlation

```{r}
# Extract MV and bulk samples
meta_raw <- sample_data(physeq4) %>% 
  as("data.frame") %>% 
  tibble::rownames_to_column("SampleID")

# Add library size
meta_raw$LibrarySize <- sample_sums(physeq4)
```

# Normalization

rle 
Tightly clustered points near y = 0: That shows consistent scaling across samples.
Symmetry around the zero line: No bias toward high or low values in specific samples.


DESeq
DESeq2 computes per-sample scaling factors that adjust for sequencing depth and compositional biases.
Put samples on the same "scale" without rarefying, retaining more power.
The distribution of values across samples is reasonably uniform.
No major skews or samples standing out with inflated variance.

clr
Values near 0 mean the taxon abundance is close to the geometric mean of that sample.
Most of the values cluster close to 0: expected in CLR — the transformation centers data per sample.
Some deep negative outliers -> taxa nearly absent in those samples (low relative abundance).


```{r}

library(DspikeIn)

# Choose grouping variable — if none, can leave it NULL or use 'Size'
group_var <- "Size"

# CLR normalization (best for compositional comparison)
result_clr <- normalization_set(physeq4, method = "clr", groups = group_var)

# DESeq normalization
result_deseq <- normalization_set(physeq4, method = "DESeq", groups = group_var)

# RLE normalization
result_rle <- normalization_set(physeq4, method = "rle", groups = group_var)

# Visualize for sanity check
library(ggplot2)
boxplot(otu_table(result_clr$dat.normed), main = "CLR")
boxplot(otu_table(result_deseq$dat.normed), main = "DESeq2 size factors")
boxplot(otu_table(result_rle$dat.normed), main = "RLE")

```
# Visualizing library Size Differences.

```{r paired-correlation-DESeq2, message=FALSE, warning=FALSE}

library(dplyr)
library(tidyr)
library(ggplot2)
library(microbiome)

# Step 1: Use DESeq2-normalized phyloseq object
physeq7 <- result_clr$dat.normed

# Step 2: Extract metadata + add normalized library sizes
meta_raw <- microbiome::meta(physeq7)
meta_raw$LibrarySize <- sample_sums(physeq7)

# Step 3: Create paired dataset (mean per mv/bulk per KFTno)
paired_df <- meta_raw %>%
  filter(Size %in% c("mv", "bulk")) %>%
  select(KFTno, Size, LibrarySize) %>%
  group_by(KFTno, Size) %>%
  summarise(LibrarySize = mean(LibrarySize), .groups = "drop") %>%  
  pivot_wider(names_from = Size, values_from = LibrarySize) %>%
  filter(!is.na(mv), !is.na(bulk))

# Step 4: Correlation test
cor_result <- cor.test(paired_df$mv, paired_df$bulk, method = "spearman")
print(cor_result)

# Step 5: Visualization of MV vs Bulk read counts
paired_df_long <- paired_df %>%
  pivot_longer(cols = c("mv", "bulk"), names_to = "Extraction", values_to = "LibrarySize")

ggplot(paired_df_long, aes(x = as.factor(KFTno), y = LibrarySize, fill = Extraction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Normalized Library Size Comparison (MV vs Bulk)",
       x = "KFT Number", y = "Normalized Read Count", fill = "Extraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +scale_fill_manual(values = DspikeIn::color_palette$mix_MG)


```


# Rarefaction 
Step 3: Rarefy dataset for fair comparisons
visually assess sample richness and sequencing depth sufficiency

```{r}
otu.rare <- otu_table(physeq4)
otu.rare <- as.data.frame(t(otu.rare))
sample_names <- rownames(otu.rare)

library(vegan)

# Fix: Set cex to a small positive number (like 0.5), or omit it
vegan::rarecurve(otu.rare, step = 100, cex = 0.5, label = TRUE)


```


# Choosing Rarefaction Depth Based on Sequencing Depth Distribution

```{r}

# Example 1: Rarefy to 4000 reads 
physeq_rare_4000 <- rarefy_even_depth(
  physeq4,
  rngseed = 02032020,
  sample.size = 4000,
  verbose = FALSE
)


saveRDS(physeq_rare_4000, "physeq_rare_4000rared.rds")

```

# Ordination

```{r}
# Create ordination object if not already done
ord <- ordinate(physeq_rare_4000, method = "PCoA", distance = "bray")

# Ensure consistent factor levels for plotting
physeq_rare_4000@sam_data$KFTno <- factor(
  physeq_rare_4000@sam_data$KFTno,
  levels = c(83, 85, 92, 100, 101, 102, 114)
)

# Plot
pcoa_plot <- plot_ordination(physeq_rare_4000, ord, type = "samples", color = "KFTno", shape = "Method") +
  scale_colour_manual(values = c(
    "83" = "#CC79A7",
    "85" = "pink",
    "92" = "blue",
    "100" = "red",
    "101" = "#56B4E9",
    "102" = "#009E73",
    "114" = "#0072B2"
  )) +
  geom_point(size = 5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCoA of Bray-Curtis Distances",
    color = "KFT Number",
    shape = "Method"
  ) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

pcoa_plot


```

# Evaluating Whether Extraction Volume Affects Microbial Community Structure

```{r}

# Convert variables to factors for plotting
physeq_rare_4000@sam_data$Extracted_volume <- factor(
  physeq_rare_4000@sam_data$Extracted_volume, 
  levels = c(100, 200, 400, 1000, 50000)
)

physeq_rare_4000@sam_data$KFTno <- factor(
  physeq_rare_4000@sam_data$KFTno, 
  levels = c(83, 85, 92, 100, 101, 102, 114)
)

# Run PCoA ordination
ord <- ordinate(physeq_rare_4000, method = "PCoA", distance = "bray")

# Create the plot
pcoa_plot <- plot_ordination(physeq_rare_4000, ord, type = "samples", 
                             color = "KFTno", shape = "Method") +
  geom_point(aes(size = Extracted_volume), alpha = 0.7) +
  
  # Method shapes
  scale_shape_manual(values = c(
    "conventional" = 19,  # Circle
    "FT"           = 15,  # Square
    "FTP"          = 17,  # Triangle
    "FTPIG"        = 8    # Star
  )) +
  
  # KFT colors
  scale_color_manual(values = c(
    "83" = "#CC33CC",
    "85" = "#FFA07A",
    "92" = "#4682B4",
    "100" = "red",
    "101" = "#56B4E9",
    "102" = "#009E73",
    "114" = "#0072B2"
  )) +

  # Size mapped to Extracted Volume
  scale_size_manual(
    name = "Extracted Volume (µL)",
    values = c("100" = 3, "200" = 4, "400" = 5, "1000" = 6, "50000" = 7),
    breaks = c("100", "200", "400", "1000", "50000"),
    labels = c("100 µL", "200 µL", "400 µL", "1000 µL", "50,000 µL")
  ) +

  # Confidence ellipses per group
  stat_ellipse(aes(group = KFTno, color = KFTno), 
               type = "t", linetype = "solid", size = 0.8, level = 0.95, alpha = 0.3) +

  # Theme and labels
  labs(
    title = "PCoA of Bray-Curtis Distances",
    color = "KFT Sample Number",
    shape = "Extraction Method",
    size  = "Extracted Volume"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

# Display plot
pcoa_plot

```
# Barplot

```{r taxa_barplot_cleaned, fig.width=16, fig.height=8, message=FALSE, warning=FALSE}

library(DspikeIn)
library(ggplot2)

# Ensure sample codes are readable
physeq_rare_4000@sam_data$KFTno <- factor(
  physeq_rare_4000@sam_data$KFTno,
  levels = c(83, 85, 92, 100, 101, 102, 114)
)

bp_ab <- taxa_barplot(
  physeq_rare_4000,
  target_glom        = "Order",
  treatment_variable = "Code",
  abundance_type     = "relative",
  x_angle            = 0,
  fill_variable      = "Order",
  facet_variable     = "KFTno",
  top_n_taxa         = 30,
  palette            = color_palette$mix_MG
)

# Customize theme and improve layout
bp_ab$barplot <- bp_ab$barplot +
  facet_grid(~KFTno, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.x = unit(0.3, "cm"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  labs(
    title = "Relative Abundance of Top 30 Orders per Sample (Faceted by KFTno)",
    y = "Relative Abundance (%)",
    x = "Sample Code",
    fill = "Order"
  )

# Print the adjusted plot
print(bp_ab$barplot)

```




# Function for betadispers, Tukey HSD and Adonis


```{r}
run_beta_div_tests <- function(physeq_obj, group_var = "Method", strata_var = "KFTno") {
  library(vegan)
#   # First, install devtools if you don't have it
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# 
# # Then install pairwiseAdonis correctly
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

  library(pairwiseAdonis)
  library(microbiome)

  meta_df <- microbiome::meta(physeq_obj)
  dist_bray <- phyloseq::distance(physeq_obj, method = "bray", weighted = TRUE)

  message("Running betadisper for dispersion homogeneity")
  betadispmod <- betadisper(dist_bray, meta_df[[group_var]])
  print(anova(betadispmod))

  message("Running permutest for dispersion differences")
  perm_result <- permutest(betadispmod, pairwise = TRUE, permutations = 9999, p.adjust = "fdr")
  print(perm_result)

  if (length(unique(meta_df[[group_var]])) > 2) {
    message("Running Tukey HSD post-hoc on dispersion")
    hsd <- TukeyHSD(betadispmod)
    hsd_df <- as.data.frame(hsd$group)
    hsd_df$Significance <- cut(hsd_df$`p adj`,
                               breaks = c(0, 0.001, 0.01, 0.05, 1),
                               labels = c("***", "**", "*", "ns"))
    print(hsd_df)
  }

  message("Running PERMANOVA (adonis2)")
  formula <- as.formula(paste("dist_bray ~", group_var))
  adonis_res <- adonis2(formula, data = meta_df, strata = meta_df[[strata_var]], permutations = 9999)
  print(adonis_res)

  message("Running pairwise PERMANOVA with Holm correction")
  pairwise_res <- pairwise.adonis2(formula,
                                   data = meta_df,
                                   strata = meta_df[[strata_var]],
                                   p.adjust = "holm")
  print(pairwise_res)
}

```

# Subset

```{r}
# Set rarefied object
physeq_rare <- physeq_rare_4000

# Subset 1: Early September to Mid-November
subset_2_30Sep_15Nov <- subset_samples(physeq_rare, KFTno %in% c(92, 100, 114))

# Subset 2: Mid-August
subset_5_12Aug <- subset_samples(physeq_rare, KFTno %in% c(83, 85))

# Subset 3: Early October
subset_4_7Oct <- subset_samples(physeq_rare, KFTno %in% c(101, 102))

```


```{r}
# Run on August subset by Method
run_beta_div_tests(subset_5_12Aug, group_var = "Method")

# Run on September–November subset by Volume
run_beta_div_tests(subset_2_30Sep_15Nov, group_var = "Extracted_volume")

# Run on 7 October subset by Method
run_beta_div_tests(subset_4_7Oct, group_var = "Method")

```

# Beta Dispersion with Tukey HSD seprately & Visualization

```{r}
# ──────────────────────────────────────────────────────────────
# Beta Dispersion Test (Bray-Curtis)
# ──────────────────────────────────────────────────────────────

bray_dist <- phyloseq::distance(physeq_rare_4000, method = "bray")
metadata <- sample_data(physeq_rare_4000)

# Beta dispersion model
betadisp_mod <- betadisper(bray_dist, group = metadata$Method)

# ANOVA and Tukey HSD
anova(betadisp_mod)

HSD <- TukeyHSD(betadisp_mod)
HSD_results <- as.data.frame(HSD$group)
HSD_results$Significance <- cut(HSD_results$`p adj`, 
                                breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                labels = c("***", "**", "*", "ns"))
print(HSD_results)
plot(HSD)

```

# Visualization of Dispersion

```{r}
# ──────────────────────────────────────────────
# Visualization of Dispersion
# ──────────────────────────────────────────────

# Extract scores and distances
betscores <- as.data.frame(scores(betadisp_mod, display = "sites"))
betscores$Distance <- betadisp_mod$distances
betscores$Method <- physeq_rare_4000@sam_data$Method
betscores$KFTno <- physeq_rare_4000@sam_data$KFTno
betscores$Extracted_volume <- as.factor(physeq_rare_4000@sam_data$Extracted_volume)

betscores$Method <- factor(betscores$Method, levels = c("conventional", "FT", "FTP", "FTPIG"))

# Plot

```


```{r}
dispersplot <- ggplot(betscores, aes(x = Method, y = Distance)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(aes(shape = Method, color = KFTno), size = 4) +  # fixed

  scale_shape_manual(values = c(
    "conventional" = 19,
    "FT" = 15,
    "FTP" = 17,
    "FTPIG" = 9
  )) +
  
  scale_colour_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +
  
  labs(shape = "Extraction Method", color = "Sample (KFTno)") +
  theme_bw() +
  ylab("Distance to Centroid") +
  ggtitle("Beta Dispersion") +
  theme(
    text = element_text(family = "sans", size = 12),  # fixed font
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    strip.background = element_rect(colour = "black", fill = "gray98")
  ) 

dispersplot

```

# PERMANOVA

```{r}

library(phyloseq)
library(vegan)
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Load devtools
library(devtools)

# Install the pairwiseAdonis package
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(cluster)

# Extract metadata and distance
dist_bray <- phyloseq::distance(subset_5_12Aug, method = "bray")
metadata <- microbiome::meta(subset_5_12Aug)

# PERMANOVA for Method
adonis_method <- adonis2(dist_bray ~ Method, data = metadata, strata = metadata$KFTno, permutations = 9999)
print(adonis_method)

# Pairwise PERMANOVA for Method
pairwise_method <- pairwise.adonis2(
  dist_bray ~ Method,
  data = metadata,
  strata = metadata$KFTno,
  p.adjust = "holm"
)
print(pairwise_method)

# PERMANOVA for Volume
adonis_volume <- adonis2(dist_bray ~ Extracted_volume, data = metadata, strata = metadata$KFTno, permutations = 9999)
print(adonis_volume)

# Pairwise PERMANOVA for Volume
pairwise_volume <- pairwise.adonis2(
  dist_bray ~ Extracted_volume,
  data = metadata,
  strata = metadata$KFTno,
  p.adjust = "holm"
)
print(pairwise_volume)

```

# DESeq2 Differential Abundance

```{r}
# DESeq2 Analysis Pipeline

# ───────────────────────────────────────────────
# Step 1: Subset your phyloseq object to compare only FTPIG vs conventional
# ───────────────────────────────────────────────
Subset_Dseq <- subset_samples(physeq_rare_4000, Method %in% c("conventional", "FTPIG"))
Subset_Dseq <- prune_samples(sample_sums(Subset_Dseq) > 0, Subset_Dseq)  # remove empty

# Optional sanity check
table(sample_data(Subset_Dseq)$Method)

# ───────────────────────────────────────────────
# Step 2: Convert to DESeq2 object and run DESeq
# ───────────────────────────────────────────────
dds <- phyloseq_to_deseq2(Subset_Dseq, ~ Method)
dds <- DESeq(dds)

# ───────────────────────────────────────────────
# Step 3: Variance Stabilizing Transformation (VST)
# ───────────────────────────────────────────────
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Create new phyloseq object from VST data
otu_vst <- otu_table(assay(vsd), taxa_are_rows = TRUE)
ps_vst <- phyloseq(otu_vst, tax_table(Subset_Dseq), sample_data(Subset_Dseq))

# ───────────────────────────────────────────────
# Step 4: Run differential abundance with visualization via DspikeIn
# ───────────────────────────────────────────────
results_DESeq2 <- perform_and_visualize_DA(
  obj = ps_vst,
  method = "DESeq2",
  group_var = "Method",
  contrast = c("FTPIG", "conventional"),
  output_csv_path = "results_DESeq2_all.csv",
  target_glom = "Genus",
  significance_level = 0.05
)

# Inspect significant taxa
head(results_DESeq2$results)
results_DESeq2$bar_plot
results_DESeq2$plot
results_DESeq2$obj_significant

# Save filtered results
write.csv(results_DESeq2$results, "results_DESeq2_significant.csv")

# ───────────────────────────────────────────────
# inspect taxonomy of significant taxa
# ───────────────────────────────────────────────
tax_data <- as.data.frame(tax_table(results_DESeq2$obj_significant))
head(tax_data)

rm(dds, vsd, otu_vst)

```

# Alpha Diversity (Shannon)

```{r fig.width=10, fig.height=6, message=FALSE, warning=FALSE}

# ───────────────────────────────────────────────
# Step 1: Estimate Shannon Diversity from rarefied phyloseq object
# ───────────────────────────────────────────────
diversity <- estimate_richness(physeq_rare_4000, measures = "Shannon")
diversity_data <- cbind(diversity, microbiome::meta(physeq_rare_4000))

# ───────────────────────────────────────────────
# Step 2: Test normality of Shannon diversity
# ───────────────────────────────────────────────
shapiro_result <- shapiro.test(diversity_data$Shannon)
print(shapiro_result)

# scale for statistical comparison or plot normalization
diversity_data$Shannon_scaled <- scale(diversity_data$Shannon)

# Build the plot
bp_shannon <- ggplot(diversity_data, aes(x = KFTno, y = Shannon, fill = as.factor(KFTno))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(
    aes(
      color = as.factor(KFTno),
      shape = Method,
      fill = as.factor(KFTno),
      size = Extracted_volume
    ),
    stroke = 1.5,
    alpha = 0.7,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)
  ) +
  facet_grid(~KFTno, scales = "free_x") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  labs(
    x = "KFT Sample Number",
    y = "Shannon Diversity Index",
    title = "Shannon Diversity Boxplots with Extracted Volume Mapping"
  ) +

  # Fill colors for KFTno
  scale_fill_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +

  # Outline color (black for all)
  scale_colour_manual(values = rep("black", length(unique(diversity_data$KFTno)))) +

  # Shapes for Method
  scale_shape_manual(values = c(
    "conventional" = 21,
    "FT" = 22,
    "FTP" = 24,
    "FTPIG" = 23,
    "FTPIG-X10" = 25,
    "FTPIG-ch-X10" = 21
  )) +

  # Size scale for Extracted Volume
  scale_size_manual(
    name = "Extracted Volume (µL)",
    values = c("100" = 3, "200" = 5, "400" = 7, "1000" = 9, "50000" = 12),
    breaks = c("100", "200", "400", "1000", "50000"),
    labels = c("100 µL", "200 µL", "400 µL", "1000 µL", "50000 µL")
  )

# Print the plot
print(bp_shannon)

# Optional: Save to file
# ggsave("shannon_plot_final.pdf", bp_shannon, width = 10, height = 6)


```

# LMM Models
Fixed effects: KFTno, Method
Random effect: (1 | Month:Temp) (nested environmental covariates)

```{r}
# Perform pairwise comparisons for Method
# Load necessary libraries
library(lme4)
library(lmerTest)
library(emmeans)

# Ensure categorical variables are factors
diversity_data$KFTno <- as.factor(diversity_data$KFTno)
diversity_data$Method <- as.factor(diversity_data$Method)
diversity_data$Month <- as.factor(diversity_data$Month)
diversity_data$Temp <- as.factor(diversity_data$Temp)
diversity_data$Purification <- as.factor(diversity_data$Purification)
diversity_data$Extracted_volume <- as.factor(diversity_data$Extracted_volume)

# ================================
# Model definitions
# ================================

# Null model: Random effect of KFTno only
lmm4 <- lmer(Shannon ~ 1 + (1 | KFTno), data = diversity_data)

# Model with Method, KFTno, and Month:Temp interaction as random effect
lmm6 <- lmer(
  Shannon ~ KFTno + Method + (1 | Month:Temp),
  data = diversity_data
)

# Model with extra random effect for purification
lmm5 <- lmer(
  Shannon ~ KFTno + Method + (1 | Month:Temp) + (1 | Purification),
  data = diversity_data
)

# Model with Extracted Volume + Month:Temp
lmm3 <- lmer(
  Shannon ~ KFTno + Method + (1 | Extracted_volume) + (1 | Month:Temp),
  data = diversity_data
)

# ================================
# Model comparison
# ================================
AIC(lmm3, lmm4, lmm5, lmm6)
BIC(lmm3, lmm4, lmm5, lmm6)

# ================================
# Summary & ANOVA of best model (e.g., lmm6)
# ================================
summary(lmm6)
anova(lmm6)

# ================================
# Pairwise comparisons
# ================================

# 1. Method overall
pairwise_method <- emmeans(lmm6, pairwise ~ Method, adjust = "holm")
summary(pairwise_method)

# 2. Method within each KFTno (interaction)
pairwise_KFT_method <- emmeans(lmm6, pairwise ~ KFTno:Method, adjust = "holm")
summary(pairwise_KFT_method)

# Optional: visualize contrasts
plot(emmeans(lmm6, "Method"), comparisons = TRUE)

```
# Shannon Diversity Across Methods and Timepoints with Both Boxplots and Violin Plots

```{r}

library(ggplot2)

ggplot(diversity_data, aes(x = Method, y = Shannon)) +
geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), alpha = 0.15) +
  
  # color = KFTno
  geom_jitter(aes(color = as.factor(KFTno)), width = 0.2, size = 3, alpha = 0.7) +
  
  # Mean points
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5, color = "black", fill = "white") +
  stat_summary(
    fun.data = function(x) {
      mean_val <- mean(x)
      se_val <- sd(x) / sqrt(length(x))
      data.frame(y = mean_val, ymin = mean_val - se_val, ymax = mean_val + se_val)
    },
    geom = "errorbar", width = 0.25, color = "black"
  ) +
  
  # Color palette for KFTno
  scale_colour_manual(values = c(
    "100" = "red", 
    "101" = "#56B4E9",  
    "102" = "#009E73",  
    "11"  = "#F0E442",  
    "114" = "#0072B2",  
    "31"  = "#D55E00",  
    "83"  = "#CC33CC",  
    "84"  = "#CC79A7",  
    "85"  = "#FFA07A",  
    "86"  = "#FDE725",  
    "90"  = "orange",  
    "92"  = "blue",
    "95"  = "darkgreen"
  )) +
  
  labs(
    title = "Shannon Diversity Index by Method and KFT Sample",
    x = "Extraction Method", 
    y = "Shannon Diversity Index",
    color = "KFT Sample"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12)
  )

```

# session info

```{r}

sessionInfo()
```
