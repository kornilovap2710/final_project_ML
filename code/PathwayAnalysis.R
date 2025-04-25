if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
if (!requireNamespace("fgsea", quietly = TRUE)) 
  BiocManager::install("fgsea")

library(fgsea)
library(ggplot2)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 2){
  stop("Wrong number of command line input parameters. Please check.")
}

file_path <- args[1]
# comparison <- args[2]
outdir <- args[2]

# Read data from the file
gene_data <- read.csv(file_path, header = TRUE, sep = '\t')

ranked_indices_desc <- order(gene_data$logFC, decreasing = TRUE)
ranked_genes_desc <- gene_data$GeneID[ranked_indices_desc]
log2foldchanges_desc <- gene_data$logFC[ranked_indices_desc]

stats_named_desc <- setNames(log2foldchanges_desc, ranked_genes_desc)
stats_named_desc

# msigdb_pathways <- gmtPathways("h.all.v7.1.symbols.gmt") # holemark set
kegg_pathways <- gmtPathways("../files/c2.cp.kegg.v7.1.symbols.gmt")
reactome_pathways <- gmtPathways("../files/c2.cp.reactome.v7.1.symbols.gmt")
combined_gene_sets <- c(kegg_pathways, reactome_pathways)

fgseaRes <- fgseaMultilevel(pathways = combined_gene_sets,
                  stats = stats_named_desc,
                  minSize = 15,
                  maxSize = 500,
                  scoreType = "std")

pathways <- fgseaRes$pathway
NES <- fgseaRes$NES
pval <- fgseaRes$pval
padj <- fgseaRes$padj
log2err <- fgseaRes$log2err
ES <- fgseaRes$ES
size <- fgseaRes$size

# Create a data frame for the bubble plot
bubble_data <- data.frame(pathways, pval, padj, log2err, ES, NES, size)

# top_10_pathways <- head(bubble_data[order(abs(bubble_data$NES), decreasing = TRUE), ], 10) # better to use FDR, if you want to use the enrichment score, you have to use the normalized ones

# ggplot(top_10_pathways, aes(x = NES, y = pathways, size = -log10(pval))) +
#   geom_point(aes(color = NES), alpha = 0.7) +
#   scale_size_continuous(range = c(3, 15)) +
#   labs(title = "Bubble Plot of Pathway Enrichment",
#        x = "Normalized Enrichment Score (NES)",
#        y = "Pathways",
#        size = "Significance (-log10 p-value)",
#        color = "NES")

# bubble_plot_file_name <- paste0("UC_vs_HC_average_6_top10_pathway_bubble.png")
# ggsave(file.path(outdir, bubble_plot_file_name), plot = last_plot(), width = 30, height = 30, units = "in")

pathways_file_name <- paste0("UC_vs_HC_only_GPCR_average_6_pathways.tsv")
write.table(bubble_data, file.path(outdir, pathways_file_name), sep = '\t', row.names = FALSE, quote = FALSE)
