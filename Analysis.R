# ============================================
# Proteomics Differential Expression Analysis
# Author: Victoria Djana
# Description: Differential expression analysis
# using limma on fictional data
# ============================================

library(limma)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ── 1. LOAD DATA ──
data <- read.csv("data/proteomics_data.csv", row.names = 1)
cat("Dataset loaded:", nrow(data), "proteins ×", ncol(data), "samples\n")

# ── 2. DEFINE GROUPS ──
groups <- factor(c("CTRL", "CTRL", "CTRL", "TREAT", "TREAT", "TREAT"))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# ── 3. NORMALIZATION ──
data_log <- log2(data + 1) # Log2 transformation
cat("Log2 normalization applied\n")

# ── 4. DIFFERENTIAL EXPRESSION — limma ──
fit <- lmFit(data_log, design)
contrast <- makeContrasts(TREAT - CTRL, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, number = Inf, adjust = "BH")
results$protein <- rownames(results)

# Significant proteins (FDR < 0.05, |logFC| > 1)
sig <- results %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("\nSignificant proteins found:", nrow(sig), "\n")

# Save results
write.csv(results, "results/differential_expression.csv")
cat("Results saved.\n")

# ── 5. VOLCANO PLOT ──
results$significance <- "Not significant"
results$significance[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Up-regulated"
results$significance[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down-regulated"

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
geom_point(alpha = 0.7, size = 2) +
scale_color_manual(values = c(
"Not significant" = "grey",
"Up-regulated" = "#E74C3C",
"Down-regulated" = "#3498DB"
)) +
geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
labs(
title = "Volcano Plot — TREAT vs CTRL",
x = "log2 Fold Change",
y = "-log10 Adjusted p-value"
) +
theme_bw()

ggsave("results/volcano_plot.png", dpi = 150)
cat("Volcano plot saved.\n")

# ── 6. HEATMAP — TOP 20 PROTEINS ──
top20 <- head(results[order(results$adj.P.Val), ], 20)
heatmap_data <- as.matrix(data_log[top20$protein, ])

pheatmap(
heatmap_data,
scale = "row",
cluster_rows = TRUE,
cluster_cols = TRUE,
color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(50),
main = "Top 20 Differential Proteins",
filename = "results/heatmap_top20.png"
)
cat("Heatmap saved.\n")

cat("\nAnalysis complete. Results saved in results/\n")
