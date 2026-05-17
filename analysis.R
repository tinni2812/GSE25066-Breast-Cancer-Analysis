# -------------------------
# 1. Load DEGs table
# -------------------------
deg <- read.table("GSE25066.top.table.tsv", header=TRUE, sep="\t")
head(deg)  # Show first few rows




# -------------------------
# 2. Filter significant DEGs
# -------------------------
# Criteria: adj.P.Val < 0.05 & |logFC| >= 1
deg_filtered <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1, ]
head(deg_filtered)



# -------------------------
# 3. Volcano Plot
# -------------------------
# Assign colors based on up/down-regulation
deg$color <- "black"
deg$color[deg$logFC > 1 & deg$adj.P.Val < 0.05] <- "red"    # upregulated
deg$color[deg$logFC < -1 & deg$adj.P.Val < 0.05] <- "blue"  # downregulated

library(ggplot2)
ggplot(deg, aes(x=logFC, y=-log10(adj.P.Val), color=color)) +
  geom_point() +
  scale_color_identity() +
  theme_minimal() +
  labs(title="Volcano Plot of DEGs",          # plot title
       x="Log2 Fold Change",                 # x-axis label
       y="-Log10 Adjusted P-value")          # y-axis label





# -------------------------
# 4. Extract top 50 DEGs
# -------------------------
top50 <- head(deg_filtered[order(deg_filtered$adj.P.Val), ], 50)
head(top50)



# -------------------------
# 5. Prepare expression matrix for heatmap
# -------------------------
top50_genes <- top50$Gene.symbol

# Assuming geo_obj already loaded and expr_mat exists
expr_mat <- exprs(geo_obj[[1]])   # Full expression matrix
dim(expr_mat)

# Extract expression for top 50 genes
top50_ids <- top50$ID
expr_top50 <- expr_mat[rownames(expr_mat) %in% top50_ids, ]
dim(expr_top50)




# -------------------------
# 6. Plot Heatmap of Top 50 DEGs
# -------------------------
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# Map probe IDs to gene symbols
gene_names <- top50$Gene.symbol[match(rownames(expr_top50), top50$ID)]

# Replace missing symbols with probe IDs
gene_names[is.na(gene_names)] <- rownames(expr_top50)

# Assign gene symbols as row names
rownames(expr_top50) <- gene_names

# Remove any rows with NA values in expression data
expr_top50_clean <- expr_top50[complete.cases(expr_top50), ]

# Plot heatmap
pheatmap(expr_top50_clean,
         scale = "row",           # scale each gene row for better visualization
         show_colnames = FALSE,   # hide column/sample names
         fontsize_row = 5,        # small readable font for row names
         annotation_col = annotation_col,  # optional column annotation
         main = "Heatmap of Top 50 DEGs")  # heatmap title




# -------------------------
# 7. GO Enrichment Analysis
# -------------------------
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

# Clean gene list for enrichment
genes <- top50$Gene.symbol
genes <- na.omit(genes)
genes <- genes[genes != ""]
genes <- unique(genes)

# Map gene symbols to Entrez IDs
gene_ids <- bitr(genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# Perform GO Biological Process enrichment
go_result <- enrichGO(gene         = gene_ids$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05)

# View top results
head(as.data.frame(go_result))

# Load enrichplot for visualization
library(enrichplot)

# Horizontal bar plot of top 10 GO terms
barplot(go_result, showCategory = 10, title = "Top 10 GO Biological Processes")




# -------------------------
# 8. KEGG Pathway Enrichment
# -------------------------
# KEGG enrichment (can run independently)
kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

# Plot KEGG results
barplot(kegg_result, showCategory=10, title="Top 10 KEGG Pathways")
barplot(kegg_result, showCategory=10, title="Top 10 KEGG Pathways")
