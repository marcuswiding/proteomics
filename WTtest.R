library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readr)
library(ggrepel)


data <- read_tsv('WTtest.tsv')
data_clean <- na.omit(data)
amountrows <- nrow(data)-nrow(data_clean)
print(amountrows)
nrow(data_clean)


#PCA plot
pcadata <- log(data_clean[ ,15:22]) #Isolates intensity data, logartihms it

pca <- prcomp(t(pcadata)) #Generates PCA using TRANSFORMED data (otherwise you get a point per protein = 4500 points)
summary(pca)
pca_scores <- as.data.frame(pca$x) #Dataframe with all PCA scores
pca_scores$Sample <- c(1,2,3,4,5,6,7,8)
pca_scores$Group <- c("A","A","A","A","B","B","B","B")

ggplot(pca_scores, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
  geom_point(size = 3) +   # points
  geom_text(vjust = -0.5) +                    # labels
  xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal()
ggsave("WTPCA.png", width = 10, height = 6, dpi = 300)

#Volcano plot
volcanodata <- log2(data_clean[,15:22])
log2FC <- rowMeans(volcanodata[,5:8]) - rowMeans(volcanodata[,1:4])

p_values <- apply(volcanodata, 1, function(row) {
  t.test(row[1:4], row[5:8], var.equal = FALSE)$p.value #var.equal = FALSE ensures Welsh t-test used (standard for proteomics), $p.value gives just p.value
})

#Adjust for FDR using Benjamini-Hochberg
fdr_values <- p.adjust(p_values, method = "BH")

results <- data.frame( 
  Protein_description = data_clean$PG.ProteinDescriptions,
  log2FC = log2FC,
  p_values =p_values,
  FDR = fdr_values
  )

#Takes negative log10 of FDR results
results$neg_log10_fdr <- -log10(results$FDR)


#Plotting:
# Thresholds
fdr_cutoff <- 0.05
fc_cutoff <- 1 # Absolute log2 fold change of 1

# Categorize proteins based on FDR instead of raw p-value
results <- results %>%
  mutate(
    Expression = case_when(
      log2FC >= fc_cutoff & FDR <= fdr_cutoff ~ "Upregulated",
      log2FC <= -fc_cutoff & FDR <= fdr_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

significant_proteins <- results %>%
  filter(Expression != "Not Significant")

# Draw the updated Volcano Plot
ggplot(results, aes(x = log2FC, y = neg_log10_fdr, color = Expression, label = Protein_description)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", 
                                "Not Significant" = "grey", 
                                "Upregulated" = "red")) +
  # Add threshold lines (y-intercept is now based on FDR)
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_cutoff), col = "black", linetype = "dashed") +
  labs(
    title = "Proteomics Volcano Plot coomparing WT OD0.5 and OD1 (FDR Corrected)",
    x = "Log2 Fold Change",
    y = "-Log10(FDR)"
  ) +
  # Protein LABELS ---
  geom_text_repel(data = top_proteins, 
                  aes(label = Protein_description), 
                  color = "black",        # Keep text black so it's readable
                  size = 3,               # Adjust text size
                  box.padding = 0.5,      # Add space around the text
                  max.overlaps = Inf) +   # Force it to draw all 10 labels
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )  
ggsave("WTvolcanoplot.png", width = 10, height = 6, dpi = 300)
