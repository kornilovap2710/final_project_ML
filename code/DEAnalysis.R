# Install and load the limma package
if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma")

library(limma)
library(edgeR)

expression_data <- read.table("../results//UC_vs_HC_23_07_2024/UC_without_placebo_vs_HC_expression_datafile.tsv", header = TRUE, row.names = 1, sep = "\t")
sample_info <- read.table("../results/UC_vs_HC_23_07_2024/UC_without_placebo_vs_HC_metadata.tsv", header = TRUE, sep = "\t")

average_counts <- rowMeans(expression_data)
filtered_expression_data <- expression_data[average_counts > 4, ]

design <- model.matrix(~ sample_info$Condition)

fit <- lmFit(filtered_expression_data, design)

fit <- eBayes(fit)

results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")
write.table(results, file = "../results/UC_vs_HC_23_07_2024/UC_without_placebo_vs_HC_DEAresults_average_4.tsv", sep = '\t', quote = FALSE)

png("../results/UC_vs_HC_23_07_2024/UC_without_placebo_vs_HC_PCA_plot_average_4.png", width = 1024, height = 768)
pca_data <- prcomp(t(filtered_expression_data))
condition_numeric <- as.numeric(factor(sample_info$Condition))
plot(pca_data$x[, 1], pca_data$x[, 2], col = condition_numeric, pch = 16, cex = 2, main = "PCA Plot", xlab = "PC1", ylab = "PC2")
condition_numeric
legend("topright", legend = c("HC", "UC"), col = 1:7, pch = 16, title = "Condition")
dev.off()

# pdf("files/GPCR/UC_vs_HC_expression_density_plot_average_6.pdf")
# title <- "Expression density"
# plotDensities(filtered_expression_data, group = condition_numeric, main = title, legend = FALSE)
# legend("topright", legend = c("HC", "UC"), col = 1:2, pch = 16, title = "Condition")
# dev.off()
