library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readr)
library(ggrepel)
library(limma)
library(tibble)
library(scales)
library(cowplot)
library(ggforce)
library(ggnewscale)
library(ggh4x)

# Formatting preamble
# ----
theme_marcus<- function() {
  theme_classic(base_size = 12, base_family = "serif") %+replace%
    theme(
      # Solid black L-shaped axes and ticks
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      
      # Black text for axes
      axis.text = element_text(color = "black", size = 20),
      axis.title = element_text(color = "black", size = 20, face = "bold"),
      
      # Clean legend with no boxes
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 20, color = "black"),
      legend.title = element_text(size = 20, face = "bold", color = "black"),
      
      # Remove all grid lines and background colors
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Center the plot title
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      
      # Small margin around the plot
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}

# ----
# Colors and levels
# ----
strain_colors <- c(
  "Wild Type"         = "#3C5488", # NPG Dark Blue (Often used for controls)
  "Evolved Wild Type" = "#8491B4", # NPG Muted Blue/Grey
  "ΔGFD2"             = "#E64B35", # NPG Red
  "ΔYAP1"             = "#4DBBD5", # NPG Light Blue
  "ΔSSN2"             = "#00A087", # NPG Teal/Green
  "ΔSSN8"             = "#F39B7F", # NPG Peach/Orange
  "ΔSSN2/SSN8"        = "#7E6148", # NPG Brown
  "ΔINM1"             = "#8c0f0f"  # NPG Dark Red
)

strain_levels <- c(
  "Wild Type",
  "Evolved Wild Type",
  "ΔGFD2",
  "ΔYAP1",
  "ΔSSN2",
  "ΔSSN8",
  "ΔSSN2/SSN8",
  "ΔINM1")

comparison_to_strain <- c(
  "BvsA" = "Evolved Wild Type",
  "CvsA" = "ΔGFD2",
  "DvsA" = "ΔYAP1",
  "EvsA" = "ΔSSN2",
  "FvsA" = "ΔSSN8",
  "GvsA" = "ΔSSN2/SSN8",
  "HvsA" = "ΔINM1"
)

# ----
# Functions
# ----
impute_perseus <- function(x, width = 0.3, downshift = 2.5) {
  missing_idx <- is.na(x)
  if(sum(missing_idx) == 0) return(x) 
  
  sample_mean <- mean(x, na.rm = TRUE)
  sample_sd <- sd(x, na.rm = TRUE)
  new_mean <- sample_mean - (downshift * sample_sd)
  new_sd <- width * sample_sd
  
  # Fills NAs with tiny variations of low numbers, preventing 0 variance
  x[missing_idx] <- rnorm(sum(missing_idx), mean = new_mean, sd = new_sd)
  return(x)
}

plotvolcano <- function(
    limma_results,
    fdr_cutoff  = 0.05,
    fc_cutoff   = 1,
    contrast    = NULL
) {
  
  results <- limma_results %>%
    rownames_to_column("Protein") %>%
    mutate(neg_log10_fdr = -log10(adj.P.Val))
  
  if (!is.null(contrast)) {
    results <- results %>%
      left_join(
        mnar_lookup %>% filter(Comparison == contrast) %>% select(PG.Genes, MNAR_type),
        by = c("Protein" = "PG.Genes")
      ) %>%
      mutate(MNAR_type = replace_na(MNAR_type, "MCAR"))
  } else {
    results <- results %>% mutate(MNAR_type = "MCAR")
  }
  
  results <- results %>%
    mutate(
      Expression = case_when(
        MNAR_type == "Only in Wild Type"     & adj.P.Val <= fdr_cutoff ~ "Only in Wild Type",
        MNAR_type == "Only in mutant" & adj.P.Val <= fdr_cutoff ~ "Only in mutant",
        logFC >=  fc_cutoff           & adj.P.Val <= fdr_cutoff ~ "Upregulated",
        logFC <= -fc_cutoff           & adj.P.Val <= fdr_cutoff ~ "Downregulated",
        TRUE                                                    ~ "Not Significant"
      ),
      Expression = factor(Expression, levels = c("Downregulated", "Not Significant",
                                                 "Only in mutant", "Only in Wild Type",
                                                 "Upregulated")),
      Point_shape = ifelse(MNAR_type %in% c("Only in Wild Type", "Only in mutant"), "MNAR", "MCAR")
    )
  
  to_label <- results %>%
    filter(adj.P.Val <= fdr_cutoff) %>%
    mutate(score = abs(logFC) * neg_log10_fdr) %>%
    filter(
      (logFC >= fc_cutoff  & rank(-score) <= 10) |
        (logFC <= -fc_cutoff & rank(-score) <= 20)
    )
  
  ggplot(results, aes(x = logFC, y = neg_log10_fdr, color = Expression, shape = Point_shape)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_manual(
      values = c(
        "Downregulated" = "blue",
        "Not Significant" = "grey70",
        "Only in mutant" = "darkgreen",
        "Only in Wild Type" = "#E25202",
        "Upregulated" = "red"
      ),
      name = "",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("MCAR" = 16, "MNAR" = 17),
      name   = ""
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = c(16, 16, 17, 17, 16)
        )
      ),
      shape = "none"
    ) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = to_label, aes(label = Protein),
      color = "black", size = 5,
      box.padding = 0.5, max.overlaps = Inf
    ) +
    labs(x = "Log2 Fold Change", y = "-Log10(FDR)", color = "") +
    theme_marcus() +
    theme(legend.position = "right")
}

# ----
# Read Data
# ----
setwd("~/Library/CloudStorage/OneDrive-Chalmers/Masters thesis/Results/Proteomics analysis/Finalresultsproteomics")
data <- read_tsv('finalproteomics.tsv')
# ----
#Data annotation/ordering
# ----
data_abundance <- data[,58:105]
# Ordering not in numerical order due to reruns of specific samples
colnames(data_abundance) <- c('1','7','13','19','25','31','37',
                              '8','14','20','26','32','38','44',
                              '9','15','21','27','33','39','45',
                              '4','10','16'             ,'40','46',
                              '11','17'              ,'41','47',
                              '12','18','24','30','36','42','48', 
                              '2','22','3','23','5','28','6','29', '34','43','35')
#Ordering
abundanceorder <- order(as.numeric(colnames(data_abundance)))
data_abundance <- data_abundance[, abundanceorder]

#Combining with protein info
analysisdata <- bind_cols(data[,c(2:4, 6:9)], data_abundance)

# ----
# Data imputation
# ----
#Makes sure that 0 values are regarded as NA, as log2(0) = -inf, not N/A potentially skewing results
analysisdata[analysisdata == 0] <- NA
#A global log2 is applied on all abundances
analysisdata <- bind_cols(analysisdata[, 1:7], log2(analysisdata[, 8:55]))

#We remove proteins that are prevalent in less than 50% of cases
keep_proteins <- rowMeans(is.na(analysisdata[, 8:55])) < 0.5
imputed_analysisdata <- analysisdata[keep_proteins, ]

#Then we replace NA values using a perseus imputation, where we replace 0s with a random variation of downshifted numbers
set.seed(42) #Makes imputation reproducable between runs
imputed_analysisdata[, 8:55] <- apply(imputed_analysisdata[, 8:55], 2, impute_perseus)

# ----
# PCA
# ----
pca <- prcomp(t(imputed_analysisdata[, 8:55]), scale. = FALSE) #GTranspose so rows = samples, cols = proteins (required format for prcomp)
pca_scores <- as.data.frame(pca$x)
summary(pca)
pca_scores$Sample <- 1:48
pca_scores$Group <- factor(
  rep(c("Wild Type", "Evolved Wild Type", "ΔGFD2", "ΔYAP1", "ΔSSN2", "ΔSSN8", "ΔSSN2/SSN8", "ΔINM1"),
      times = c(6, 6, 6, 6, 6, 6, 6, 6)),
  levels = strain_levels
)

set.seed(42)
k <- 9  # 48 samples / 3 per cluster

km <- kmeans(
  pca_scores[, c("PC1", "PC2")],
  centers  = k,
  nstart   = 50,   # many random starts = more stable solution with small clusters
  iter.max = 300   # extra iterations for convergence
)

pca_scores$Cluster <- factor(km$cluster)


ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_mark_circle(aes(group = Cluster), 
                 fill = NA, colour = "grey50",
                 linetype = "dashed", linewidth = 0.5,
                 expand = unit(3, "mm"),     # padding around points
                 radius = unit(3, "mm")) +   # corner smoothness
  geom_point(aes(colour = Group), size = 3) +
  scale_color_manual(values = strain_colors) +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")) +
  labs(colour = "Strain") +
  theme_marcus()
ggsave('FinalFigures/PCA.pdf', width = 10, height = 6, device = cairo_pdf)

# Extract loadings (top proteins driving PC)
loadings <- as.data.frame(pca$rotation)
loadings$Protein <- rownames(loadings)

# Top proteins driving PC1
loadings <- as.data.frame(pca$rotation)
loadings$Protein <- imputed_analysisdata[[2]]

loadings %>%
  select(Protein, PC2) %>%
  arrange(desc(abs(PC2))) %>%
  head(20)

top_loadings <- loadings %>%
  select(Protein, PC2) %>%
  arrange(desc(abs(PC2))) %>%
  head(20)

ggplot(top_loadings, aes(x = reorder(Protein, PC2), y = PC2)) +
  geom_col() +
  coord_flip() +
  xlab(NULL) +
  theme_marcus() +
  theme(axis.text.y = element_text(size = 8))  # shrink if crowded
# ----
# Limma (compares replicates of S1 and S2 to each other)
# ----
#Create correct data form
limmadata <- imputed_analysisdata[,8:55] %>% 
  as.matrix()
rownames(limmadata) <- imputed_analysisdata$PG.Genes

# Here the experimental design is described, in this case 16 groups with 3 samples in each
condition <- factor(rep(c('A1', 'A2', 'B1','B2','C1','C2','D1','D2','E1', 'E2', 'F1', 'F2', 'G1', 'G2', 'H1', 'H2'), 
                        times = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)))
design <- model.matrix(~0 + condition) #Creates null matrix 
colnames(design) <- levels(condition) #Assigns each strain a column
stopifnot(ncol(limmadata) == length(condition))

# Fit
linear_fit <- lmFit(limmadata, design) #fits linear model for each strain

# Contrasts
contrast_matrix <- makeContrasts(
  B1vsA1 = B1 - A1,
  C1vsA1 = C1 - A1,
  D1vsA1 = D1 - A1,
  E1vsA1 = E1 - A1,
  F1vsA1 = F1 - A1,
  G1vsA1 = G1 - A1,
  H1vsA1 = H1 - A1,
  B2vsA2 = B2 - A2,
  C2vsA2 = C2 - A2,
  D2vsA2 = D2 - A2,
  E2vsA2 = E2 - A2,
  F2vsA2 = F2 - A2,
  G2vsA2 = G2 - A2,
  H2vsA2 = H2 - A2,
  levels = design
)

fit2 <- contrasts.fit(linear_fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Automatically get all the contrast names defined ("B1vsA1", "B2vsA2", etc.)
all_comparisons <- colnames(contrast_matrix)

# Run topTable for every single contrast
results_list <- lapply(all_comparisons, function(x) {
  topTable(fit2, coef = x, number = Inf, sort.by = "P")
})

# Name the list elements with their respective contrast names
names(results_list) <- all_comparisons

# ----
# Separate MNAR (missing not at random) from MCAR (missing completely at random)
# ----
# Group column indices for each condition (3 replicates each)
condition_cols <- list(
  A1 = 1:3,   A2 = 4:6,
  B1 = 7:9,   B2 = 10:12,
  C1 = 13:15, C2 = 16:18,
  D1 = 19:21, D2 = 22:24,
  E1 = 25:27, E2 = 28:30,
  F1 = 31:33, F2 = 34:36,
  G1 = 37:39, G2 = 40:42,
  H1 = 43:45, H2 = 46:48
)

# For each contrast, flag proteins where one group is completely absent
# Returns a named list of character vectors (one per contrast)
get_mnar_proteins <- function(data, contrast_name) {
  groups  <- strsplit(contrast_name, "vs")[[1]]
  g1_cols <- condition_cols[[groups[1]]] + 7 #+7 is offset for first columns in analysisdata
  g2_cols <- condition_cols[[groups[2]]] + 7
  
  data %>%
    rowwise() %>%
    mutate(
      all_missing_g1 = all(is.na(c_across(all_of(g1_cols)))),
      all_missing_g2 = all(is.na(c_across(all_of(g2_cols)))),
      MNAR_type = case_when(
        all_missing_g1 & !all_missing_g2 ~ "Only in Wild Type",
        all_missing_g2 & !all_missing_g1 ~ "Only in mutant",
        TRUE                             ~ "MCAR"
      )
    ) %>%
    ungroup() %>%
    select(PG.Genes, MNAR_type)
}

# Build a lookup for all contrasts
mnar_list <- lapply(all_comparisons, function(comp) {
  get_mnar_proteins(analysisdata[keep_proteins, ], comp) %>%
    mutate(Comparison = comp)
})
mnar_lookup <- bind_rows(mnar_list)
# ----
# Plotting
# ----
#Plots Sample 1 and sample 2 side by side for each comparison
#Strips legend from plots
plotvolcano_nolegend <- function(...) {
  plotvolcano(...) + theme(legend.position = "none")
}

# Pair up the comparisons: B1vsA1 with B2vsA2, etc.
strains <- c("B", "C", "D", "E", "F", "G", "H")

for (s in strains) {
  comp1 <- paste0(s, "1vsA1")
  comp2 <- paste0(s, "2vsA2")
  
  p1 <- plotvolcano_nolegend(results_list[[comp1]], contrast = comp1) +
    labs(tag = "A", title = "Diauxic shift") +
    theme(
      plot.tag.position = c(0, 1),
      plot.tag          = element_text(size = 24, face = "bold", vjust = 2),
      plot.title        = element_text(size = 24, face = "bold", hjust = 0.5, vjust = 2),
      plot.margin       = margin(t = 5, r = 5, b = 20, l = 5)
    )
  
  p2 <- plotvolcano_nolegend(results_list[[comp2]], contrast = comp2) +
    labs(tag = "B", title = "Post-diauxic growth") +
    theme(
      plot.tag.position = c(0, 1),
      plot.tag          = element_text(size = 24, face = "bold", vjust = 2),
      plot.title        = element_text(size = 24, face = "bold", hjust = 0.5, vjust = 2),
      plot.margin       = margin(t = 20, r = 5, b = 5, l = 5)
    )
  
  legend_data <- data.frame(
    logFC          = c(-2,  0,   1,               1,                2),
    neg_log10_fdr  = c( 2,  1,   2,               2,                2),
    Expression     = factor(
      c("Downregulated", "Not Significant", "Only in mutant", "Only in Wild Type", "Upregulated"),
      levels = c("Downregulated", "Not Significant", "Only in mutant", "Only in Wild Type", "Upregulated")
    ),
    Point_shape    = c("MCAR", "MCAR", "MNAR", "MNAR", "MCAR")
  )
  
  legend_plot <- ggplot(legend_data,
                        aes(x = logFC, y = neg_log10_fdr,
                            color = Expression, shape = Point_shape)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(
      values = c(
        "Downregulated"     = "blue",
        "Not Significant"   = "grey70",
        "Only in mutant"    = "darkgreen",
        "Only in Wild Type" = "#E25202",
        "Upregulated"       = "red"
      ),
      name = "", drop = FALSE
    ) +
    scale_shape_manual(values = c("MCAR" = 16, "MNAR" = 17)) +
    guides(
      shape = "none",
      color = guide_legend(
        override.aes = list(shape = c(16, 16, 17, 17, 16))
      )
    ) +
    theme_marcus()
  
  legend <- get_legend(legend_plot)
  
  combined <- (p1 | p2 | wrap_elements(legend)) +
  plot_layout(widths = c(1, 1, 0.25))
  
  ggsave(
    paste0('FinalFigures/volcano_', s, 'vsA.pdf'),
    plot   = combined,
    width  = 24,
    height = 8,
    device = cairo_pdf
  )
}
# ----
# Common proteins S1 and S2 compared separately
# ----
comparison_to_strain <- c(
  "B1vsA1" = "Evolved Wild Type Diauxic shift",
  "C1vsA1" = "ΔGFD2 Diauxic shift",
  "D1vsA1" = "ΔYAP1 Diauxic shift",
  "E1vsA1" = "ΔSSN2 Diauxic shift",
  "F1vsA1" = "ΔSSN8 Diauxic shift",
  "G1vsA1" = "ΔSSN2/SSN8 Diauxic shift",
  "H1vsA1" = "ΔINM1 Diauxic shift",
  "B2vsA2" = "Evolved Wild Type Post-diauxic growth",
  "C2vsA2" = "ΔGFD2 Post-diauxic growth",
  "D2vsA2" = "ΔYAP1 Post-diauxic growth",
  "E2vsA2" = "ΔSSN2 Post-diauxic growth",
  "F2vsA2" = "ΔSSN8 Post-diauxic growth",
  "G2vsA2" = "ΔSSN2/SSN8 Post-diauxic growth",
  "H2vsA2" = "ΔINM1 Post-diauxic growth"
)
strain_levels <- c(
  "Evolved Wild Type Diauxic shift", "Evolved Wild Type Post-diauxic growth",
  "ΔGFD2 Diauxic shift",      "ΔGFD2 Post-diauxic growth",
  "ΔYAP1 Diauxic shift",      "ΔYAP1 Post-diauxic growth",
  "ΔSSN2 Diauxic shift",      "ΔSSN2 Post-diauxic growth",
  "ΔSSN8 Diauxic shift",      "ΔSSN8 Post-diauxic growth",
  "ΔSSN2/SSN8 Diauxic shift", "ΔSSN2/SSN8 Post-diauxic growth",
  "ΔINM1 Diauxic shift",      "ΔINM1 Post-diauxic growth"
)

# Isolates significant proteins 
sig_proteins <- lapply(names(results_list), function(prot_comp){
  results_list[[prot_comp]] %>%
    rownames_to_column('PG.Genes') %>% 
    filter(adj.P.Val <= 0.05, abs(logFC) >= 1) %>%
    mutate(Comparison = prot_comp)
})
names(sig_proteins) <- names(results_list)

#Combines common proteins
common_sig_proteins <- bind_rows(sig_proteins) %>%
  group_by(PG.Genes) %>%
  filter(n_distinct(Comparison) >= 2) %>% #How many comparisons the proteins needs to be found in
  ungroup() %>%
  arrange(PG.Genes, Comparison)

#Creates summary that can be exported
common_sig_summary <- common_sig_proteins %>%
  left_join(
    imputed_analysisdata %>% select(PG.Genes, PG.ProteinGroups), #adds proteingroups
    by = "PG.Genes") %>% 
  group_by(PG.Genes) %>%
  summarise(protein_groups     = first(PG.ProteinGroups),
            n_comparisons = n_distinct(Comparison),
            comparisons = paste(comparison_to_strain[Comparison], collapse = ', '),
            directions = paste(ifelse(logFC > 0, 'Up', 'Down'), collapse = ', '),
            up_logFC      = paste(round(logFC[logFC > 0], 2), collapse = ', '),
            down_logFC    = paste(round(logFC[logFC < 0], 2), collapse = ', '),
            .groups = 'drop') %>%
  arrange(desc(n_comparisons))

#Data for plotting with extra sorting and filtering
common_sig_plot <- common_sig_proteins %>%
  group_by(PG.Genes) %>%
  filter(max(abs(logFC)) >= 2) %>%
  filter(n_distinct(Comparison) >= 2) %>%
  mutate(n_comparisons = n_distinct(Comparison)) %>%
  ungroup() %>%
  distinct(PG.Genes, n_comparisons) %>%
  slice_max(n_comparisons, n = 30) %>%
  semi_join(x = common_sig_proteins, by = "PG.Genes") %>%
  left_join(
    common_sig_proteins %>%
      group_by(PG.Genes) %>%
      summarise(n_comparisons = n_distinct(Comparison), .groups = "drop"),
    by = "PG.Genes"
  ) %>%
  left_join(
    mnar_lookup %>% select(PG.Genes, Comparison, MNAR_type),
    by = c("PG.Genes", "Comparison")
  ) %>%
  mutate(
    Strain = factor(
      gsub(" Diauxic shift| Post-diauxic growth", "", comparison_to_strain[Comparison]),
      levels = gsub(" Diauxic shift| Post-diauxic growth", "", strain_levels) %>% unique()
    ),
    Replicate = ifelse(grepl("1vs", Comparison), "Diauxic shift", "Post-diauxic growth"),
    MNAR_type = replace_na(MNAR_type, "MCAR")
  )

ggplot(common_sig_plot, aes(x = Strain, y = reorder(PG.Genes, n_comparisons))) +
  geom_tile(aes(fill = logFC), color = NA, linewidth = 0.6) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Log2 Fold Change",
    guide = guide_colorbar(
      barwidth    = 1,
      barheight   = 9,
      title.hjust = 0,
      title.theme = element_text(size = 20, face = "bold", margin = margin(b = 10))
    )
  ) +
  new_scale_fill() +
  geom_tile(
    data = subset(common_sig_plot, MNAR_type != "MCAR"),
    aes(fill = MNAR_type),
    color = NA, linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Only in Wild Type" = "#996e0c",
      "Only in mutant"    = "#009E73"
    ),
    name = ""
  ) +
  facet_wrap(~ Replicate, scales = "free_x") +
  labs(x = "Strain compared to Wild Type", y = "Protein") +
  theme_marcus() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text  = element_text(size = 14, face = "bold")
  )

ggsave('FinalFigures/proteincomparison.pdf',
  width  = 8.3,
  height = 8.3,
  device    = cairo_pdf)

# ----
# Biological process
# ----
# Function to get significant proteins with their BP for one comparison
get_bp_summary <- function(topTable_result, bp_lookup, fdr_cutoff = 0.05, fc_cutoff = 1) {
  
  # Get significant proteins
  sig_proteins <- topTable_result %>%
    rownames_to_column("PG.Genes") %>%
    filter(adj.P.Val <= fdr_cutoff, abs(logFC) >= fc_cutoff) %>%
    mutate(Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
    left_join(bp_lookup, by = "PG.Genes")
  
  # Split into long format - one row per protein-BP combination
  sig_long <- sig_proteins %>%
    filter(!is.na(PG.BiologicalProcess)) %>%
    separate_rows(PG.BiologicalProcess, sep = ",") %>%
    mutate(PG.BiologicalProcess = trimws(PG.BiologicalProcess))
  
  # Split into a named list, one element per BP term
  bp_split <- split(sig_long, sig_long$PG.BiologicalProcess)
  
  return(bp_split)
}

# First, create a lookup of protein -> biological process from your data
bp_lookup <- imputed_analysisdata %>%
  select(PG.Genes, PG.BiologicalProcess) %>%
  filter(!is.na(PG.BiologicalProcess))

# Apply to all comparisons
bp_results <- lapply(names(results_list), function(comp) {
  get_bp_summary(results_list[[comp]], bp_lookup)
})
names(bp_results) <- names(results_list)

# Flatten bp_results into one long dataframe  
bp_plot_data <- bind_rows(
  lapply(names(bp_results), function(comp) {
    bind_rows(bp_results[[comp]]) %>%
      mutate(Comparison = comp)
  })
) %>%
  group_by(Comparison, PG.BiologicalProcess, Direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(!PG.BiologicalProcess %in% c("biological_process", "Requires Ontology"))  # remove unannotated proteins

bp_plot_data_filtered <- bp_plot_data %>%
  group_by(PG.BiologicalProcess) %>%
  filter(
    n_distinct(Comparison) >= 2,
    sum(n) >= 10
  ) %>%
  ungroup() %>%
  mutate(
    Strain = factor(
      gsub(" Diauxic shift| Post-diauxic growth", "", comparison_to_strain[Comparison]),
      levels = gsub(" Diauxic shift| Post-diauxic growth", "", strain_levels) %>% unique()
    ),
    Replicate = ifelse(grepl("1vs", Comparison), "Diauxic shift", "Post-diauxic growth"),
    # Signed n: positive = upregulated, negative = downregulated
    signed_n = ifelse(Direction == "Upregulated", n, -n)
  )

top_up_bps <- bp_plot_data_filtered %>%
  filter(Direction == "Upregulated") %>%
  group_by(PG.BiologicalProcess) %>%
  summarise(total_n = sum(n)) %>%
  slice_max(total_n, n = 12, with_ties = FALSE) %>%
  pull(PG.BiologicalProcess)

top_down_bps <- bp_plot_data_filtered %>%
  filter(Direction == "Downregulated") %>%
  group_by(PG.BiologicalProcess) %>%
  summarise(total_n = sum(n)) %>%
  slice_max(total_n, n = 11, with_ties = FALSE) %>%
  pull(PG.BiologicalProcess)

bp_plot_data_extremes_split <- bp_plot_data_filtered %>%
  filter(
    (Direction == "Upregulated"   & PG.BiologicalProcess %in% top_up_bps) |
      (Direction == "Downregulated" & PG.BiologicalProcess %in% top_down_bps)
  )


bp_plot_data_extremes_split %>% 
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  ggplot(aes(x = Strain, 
             y = reorder(PG.BiologicalProcess, signed_n, sum), 
             fill = signed_n)) +          # <-- switch to signed_n
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(                   # <-- diverging scale
    low      = "darkblue",   # downregulated colour
    mid      = "white",
    high     = "firebrick",   # upregulated colour
    midpoint = 0,
    name     = "Protein count",
    labels   = abs 
  ) +
  facet_grid(Direction ~ Replicate, scales = "free_x") +
  scale_y_discrete(labels = label_wrap(50)) +
  labs(x = "Strain", y = "Biological Process") +
  theme_marcus() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text  = element_text(size = 14, face = "bold")
  )

ggsave('FinalFigures/bp_heatmap.pdf',
       width  = 10,
       height = 14,
       device = cairo_pdf)

#Creates a plot for the presentation with better readability
bp_plot_data_extremes_split %>% 
  mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))) %>%
  ggplot(aes(x = Strain, 
             y = reorder(PG.BiologicalProcess, signed_n, sum), 
             fill = signed_n)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "darkblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0,
    name     = "Protein count",
    labels   = abs 
  ) +
  facet_nested(. ~ Direction + Replicate, scales = "free_x") +
  scale_y_discrete(labels = label_wrap(50)) +
  labs(x = "Strain", y = "Biological Process") +
  theme_marcus() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text  = element_text(size = 14, face = "bold")
  )
ggsave('FinalFigures/bp_heatmappresentation.pdf',
       width  = 15,
       height = 6,
       device = cairo_pdf)

# ----
# Individual strains interesting proteins
# ----
get_interesting_proteins <- function(
    results_list,
    comparisons = c("B1vsA1", "B2vsA2", "C1vsA1", "C2vsA2", "D1vsA1", "D2vsA2",
                     "E1vsA1", "E2vsA2", "F1vsA1", "F2vsA2", "G1vsA1", "G2vsA2",
                     "H1vsA1", "H2vsA2"),
    direction = "both",
    n_up = 10,
    n_down = 20,
    fdr_cutoff = 0.05,
    lfc_cutoff = 1,
    analysisdata = NULL,
    cat_output = TRUE
) {
  
  results <- lapply(comparisons, function(comp) {
    
    df <- results_list[[comp]] %>%
      rownames_to_column("Protein") %>%
      mutate(neg_log10_fdr = -log10(adj.P.Val)) %>%
      filter(adj.P.Val <= fdr_cutoff) %>%
      mutate(score = abs(logFC) * neg_log10_fdr)
    
    filtered <- switch(direction,
                       "up"   = df %>% filter(logFC >=  lfc_cutoff & rank(-score) <= n_up),
                       "down" = df %>% filter(logFC <= -lfc_cutoff & rank(-score) <= n_down),
                       "both" = df %>% filter(
                         (logFC >=  lfc_cutoff & rank(-score) <= n_up) |
                           (logFC <= -lfc_cutoff & rank(-score) <= n_down)
                       )
    )
    
    if (!is.null(analysisdata)) {
      filtered <- filtered %>%
        left_join(
          analysisdata %>% select(PG.Genes, PG.ProteinGroups),
          by = c("Protein" = "PG.Genes")
        )
    }
    
    filtered %>% mutate(Comparison = comp)
  }) %>%
    setNames(comparisons)
  
  if (cat_output) {
    cat(paste(
      sapply(comparisons, function(x) paste(results[[x]]$Protein, collapse = ", ")),
      collapse = ", "
    ), "\n")
  }
  
  invisible(results)
}

res <- get_interesting_proteins(
  results_list  = results_list,
  comparisons = c("E1vsA1", "E2vsA2", "F1vsA1", "F2vsA2", "G1vsA1", "G2vsA2"),
  direction     = "up", 
)

# ----
# Isolate specific transmembrane proteins
# ----
transmembrane <- bp_lookup %>%
  filter(grepl("transmembrane transport", PG.BiologicalProcess)) %>%
  left_join(common_sig_summary, by = "PG.Genes") %>%
  mutate(
    up_logFC_sum = sapply(strsplit(up_logFC, ",\\s*"), function(x) sum(as.numeric(x), na.rm = TRUE))
  ) %>%
  arrange(desc(up_logFC_sum))
head(transmembrane)
