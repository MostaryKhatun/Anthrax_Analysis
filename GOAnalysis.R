# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("ggplot2")

# Load the libraries
library(clusterProfiler)
library(org.Hs.eg.db) # For human gene annotation
library(enrichplot)
library(ggplot2)

# Update packages
update.packages(ask = FALSE, checkBuilt = TRUE)

# Define the gene list
genes <- c("GULP1", "GATAD2A", "ACSS1", "LOC100506530", "CNTN4", 
           "LINC00485", "WNT6", "LOC101928461", "IL20RA", "VWC2", 
           "RNF128", "TTC6")
# Perform GO enrichment analysis
go_enrichment <- enrichGO(gene         = genes,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",  # Indicating the key type is SYMBOL
                          ont          = "CC",    # Choose BP, MF, CC, or ALL
                          pAdjustMethod = "BH",    # Benjamini-Hochberg correction for multiple testing
                          pvalueCutoff = 0.25,
                          qvalueCutoff = 0.25
)

# View and save GO enrichment results to CSV
go_enrichment_df <- as.data.frame(go_enrichment)
write.csv(go_enrichment_df, file = "GO_Enrichment_Results.csv", row.names = FALSE)

# Define the title
title <- "Cellular Component (CC)"

# Save Dotplot
dotplot_plot <- dotplot(go_enrichment, showCategory = 15) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
dotplot_plot
ggsave("Dotplot_GO_CC.jpeg", plot = dotplot_plot, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_CC.png", plot = dotplot_plot, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_CC.pdf", plot = dotplot_plot, width = 7, height = 8)

# Perform GO enrichment analysis for BP
go_enrichment_bp <- enrichGO(gene         = genes,
                             OrgDb        = org.Hs.eg.db,
                             keyType      = "SYMBOL",
                             ont          = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.25,
                             qvalueCutoff = 0.25
)

# Save results
go_enrichment_bp_df <- as.data.frame(go_enrichment_bp)
write.csv(go_enrichment_bp_df, file = "GO_Enrichment_BP.csv", row.names = FALSE)

# Title
title_bp <- "Biological Process (BP)"

# Dotplot
dotplot_bp <- dotplot(go_enrichment_bp, showCategory = 15) +
  ggtitle(title_bp) +
  theme(plot.title = element_text(hjust = 0.5))

dotplot_bp

# Save plots
ggsave("Dotplot_GO_BP.jpeg", plot = dotplot_bp, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_BP.png", plot = dotplot_bp, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_BP.pdf", plot = dotplot_bp, width = 7, height = 8)

# Perform GO enrichment analysis for MF
go_enrichment_mf <- enrichGO(gene         = genes,
                             OrgDb        = org.Hs.eg.db,
                             keyType      = "SYMBOL",
                             ont          = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.25,
                             qvalueCutoff = 0.25
)

# Save results
go_enrichment_mf_df <- as.data.frame(go_enrichment_mf)
write.csv(go_enrichment_mf_df, file = "GO_Enrichment_MF.csv", row.names = FALSE)

# Title
title_mf <- "Molecular Function (MF)"

# Dotplot
dotplot_mf <- dotplot(go_enrichment_mf, showCategory = 15) +
  ggtitle(title_mf) +
  theme(plot.title = element_text(hjust = 0.5))

dotplot_mf

# Save plots
ggsave("Dotplot_GO_MF.jpeg", plot = dotplot_mf, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_MF.png", plot = dotplot_mf, width = 7, height = 8, dpi = 1000)
ggsave("Dotplot_GO_MF.pdf", plot = dotplot_mf, width = 7, height = 8)

# Save Barplot
barplot_plot <- barplot(go_enrichment, showCategory = 15) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
barplot_plot
ggsave("Barplot_GO_MF.jpeg", plot = barplot_plot, width = 7, height = 8, dpi = 1000)
ggsave("Barplot_GO_MF.png", plot = barplot_plot, width = 7, height = 8, dpi = 1000)
ggsave("Barplot_GO_MF.pdf", plot = barplot_plot, width = 7, height = 8)

# Save Cnetplot
cnetplot_plot <- cnetplot(go_enrichment, showCategory = 10) +
  ggtitle("GO Term Network") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(20, 20, 20, 20)
  )
cnetplot_plot
ggsave("Cnetplot_GO_MF.jpeg", plot = cnetplot_plot, width = 12, height = 10, dpi = 1000)
ggsave("Cnetplot_GO_MF.png", plot = cnetplot_plot, width = 12, height = 10, dpi = 1000)
ggsave("Cnetplot_GO_MF.pdf", plot = cnetplot_plot, width = 12, height = 10)

# Calculate term similarity and Save Emapplot
go_enrichment_sim <- pairwise_termsim(go_enrichment)
emapplot_plot <- emapplot(go_enrichment_sim, showCategory = 10) +
  ggtitle("GO Term Enrichment Map") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin = margin(20, 20, 20, 20)
  )
go_enrichment_sim
ggsave("Emapplot_GO_MF.jpeg", plot = emapplot_plot, width = 12, height = 10, dpi = 1000)
ggsave("Emapplot_GO_MF.png", plot = emapplot_plot, width = 12, height = 10, dpi = 1000)
ggsave("Emapplot_GO_MF.pdf", plot = emapplot_plot, width = 12, height = 10)

