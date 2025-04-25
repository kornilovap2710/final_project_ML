# Install and load the limma package
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")

library(limma)
library(edgeR)
library(stringr)
library(dplyr)


expression_data <- read.table(
  "data/UNIFI_for_Polina/filtered_data/UC_without_placebo_vs_HC_expression_file.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
sample_info <- read.table(
  "data/UNIFI_for_Polina/filtered_data/UC_without_placebo_vs_HC_metadata.tsv",
  header = TRUE,
  sep = "\t"
)

expression_data_remission <- read.table(
  "data/UNIFI_for_Polina/filtered_data/UC_placebo_vs_HC_remission_expression_file.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

sample_info_remission <- read.table(
  "data/UNIFI_for_Polina/filtered_data/UC_placebo_vs_HC_remission_metadata.tsv",
  header = TRUE,
  sep = "\t"
)


# sample_info_remission <- read.table()

average_counts <- rowMeans(expression_data)
filtered_expression_data <- expression_data[average_counts > 4, ]

filter_expression_data <- function(expression_data, threshold = 4) {
  average_counts <- rowMeans(expression_data)
  filtered_expression_data <- expression_data[average_counts > threshold, ]
  return(filtered_expression_data)
}

clean_condition_names <- function(sample_info) {
  sample_info$Condition <- str_trim(sample_info$Condition)
  sample_info$Condition <- gsub("[ /()\\-]", "_", sample_info$Condition)
  return(sample_info)
}

# Function to filter out healthy volunteers from sample metadata and expression data
filter_samples <- function(sample_info, filtered_expression_data, condition_exclude = "Healthy", week_filter = NULL) {
  # Remove specified condition
  if (!is.null(condition_exclude) && condition_exclude != "") {
    sample_info <- sample_info %>% filter(Condition != condition_exclude)
  }

  # Filter by week if specified
  if (!is.null(week_filter)) {
    sample_info <- sample_info %>% filter(Week %in% week_filter)
  }

  # Get the filtered sample IDs
  filtered_subjects <- sample_info$SampleID

  # Filter expression data based on selected subjects
  filtered_expression <- filtered_expression_data %>%
    select(all_of(filtered_subjects))

  return(list(filtered_sample_info = sample_info,
              filtered_expression = filtered_expression))
}
  


# 
# filtered_expression_data_rem <- filter_expression_data(expression_data_remission)
# sample_info_rem <- clean_condition_names(sample_info_remission)
# filtered_results_rem <- filter_samples(sample_info_rem, filtered_expression_data_rem,condition_exclude = "Healthy",week_filter = "WEEK I-8")
# filtered_sample_info_rem <- filtered_results_rem$filtered_sample_info
# filtered_expression_rem <- filtered_results_rem$filtered_expression

filtered_expression_data_UC_H <- filter_expression_data(expression_data)
sample_info_UC_H <- clean_condition_names(sample_info)


# Preparring metadata and exp data for UC all doses comparison baseline and Week 8
filtered_results_UC_baseline_Week8 <- filter_samples(sample_info_UC,filtered_expression_data_UC_H,condition_exclude = "Healthy",week_filter = NULL)
sample_info_UC_base_week8 <- filtered_results_UC_baseline_Week8$filtered_sample_info
expression_data_UC_base_week8 <- filtered_results_UC_baseline_Week8$filtered_expression

# # Drug dosages
# doses <- unique(sample_info_UC_H$Condition) %>%
#   filter(Condition != "Healthy")

# All doses + healthy
doses <- unique(sample_info_UC_H$Condition)
doses

results_list <- list()

for (dose in doses){
  print(paste("Processing dose:", dose))
  dose_metadata <- sample_info_UC_base_week8[sample_info_UC_base_week8$Condition == dose, ]
  dose_expression <- expression_data_UC_base_week8[,dose_metadata$SampleID]
  
  design <- model.matrix(~ dose_metadata$Week)
  colnames(design)

  fit <- lmFit(dose_expression, design)

  fit <- eBayes(fit)

  results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")

  result_name <- paste0(dose,"_bl_week8")
  results_list[[result_name]] <- results
  
  results_df <- as.data.frame(results)
  results_df$"Gene.symbol" <- rownames(results_df) 
  results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
  file_name <- paste0("results/DEGs/", dose,"bl_week8_DEGs.tsv")
  write.table(results_df, file = file_name, sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)
} 

for (i in names(results_list)){
  results <- results_list[[i]]
  filt_results <- results %>%
    filter((logFC > 1 | logFC < -1) & adj.P.Val < 0.05)
  assign(i, filt_results)
}

# All doses and healthy samples comparison

for (dose in doses){
  for (condition in comparison_conditions) {
    if (dose == condition) {
      next
    }
  print(paste("Processing dose:", dose))
  dose_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == dose, ]
  comparison_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == condition, ]
  
  # Combine dose and healthy metadata to keep both
  combined_metadata <- rbind(dose_metadata, comparison_metadata)
  dose_expression <- filtered_expression_data[, combined_metadata$SampleID]
  # Ensure the reference is set to 'Healthy' for dose vs healthy comparison
  if (condition == "Healthy") {
    combined_metadata$Condition <- factor(combined_metadata$Condition, levels = c("Healthy", dose))
  } 
  else {
    # For dose vs dose comparison, we can set the current dose as the reference
    combined_metadata$Condition <- factor(combined_metadata$Condition, levels = c(dose, condition))
  }
  
  design <- model.matrix(~ combined_metadata$Condition)
  colnames(design)
  
  fit <- lmFit(dose_expression, design)
  
  fit <- eBayes(fit)
  
  results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")
  result_name <- paste0(dose, "_vs_", condition)
  results_list[[result_name]] <- results
  
  results_df <- as.data.frame(results)
  results_df$"Gene.symbol" <- rownames(results_df) 
  results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
  file_name <- paste0("results/DEGs/", dose,"_vs_",condition,"_DEGs.tsv")
  # write.table(results_df, file = file_name, sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)
  } 
}
for (i in names(results_list)){
  results <- results_list[[i]]
  filt_results <- results %>%
    filter((logFC > 1 | logFC < -1) & adj.P.Val < 0.05)
  assign(i, filt_results)
}




# Initialize an empty list to store results
results_list <- list()

# Loop through each dose and make comparisons with Healthy and other doses
for (ref_dose in doses) {
  
  # **First comparison: Healthy vs ref_dose**
  if (ref_dose != "Healthy") {  # Skip if ref_dose is Healthy because we already do this in another part
    
    print(paste("Processing comparison:", "Healthy vs", ref_dose))
    
    # Subset metadata for the reference dose and Healthy condition
    ref_dose_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == ref_dose, ]
    healthy_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == "Healthy", ]
    
    # Combine metadata for the two conditions
    combined_metadata <- rbind(ref_dose_metadata, healthy_metadata)
    
    # Subset the expression data for the samples in both conditions
    dose_expression <- filtered_expression_data_UC_H[, combined_metadata$SampleID]
    
    # Set Healthy as the reference
    combined_metadata$Condition <- factor(combined_metadata$Condition, levels = c("Healthy", ref_dose))
    
    # Create the design matrix based on the combined metadata
    design <- model.matrix(~ combined_metadata$Condition)
    
    # Fit the linear model
    fit <- lmFit(dose_expression, design)
    
    # Apply empirical Bayes moderation
    fit <- eBayes(fit)
    
    # Extract results for the comparison (Healthy vs ref_dose)
    results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")
    
    # Store results in the results list
    result_name <- paste0("Healthy_vs_", ref_dose)
    results_list[[result_name]] <- results
    
    # Prepare results for output
    results_df <- as.data.frame(results)
    results_df$"Gene.symbol" <- rownames(results_df) 
    results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
    
    # # Save results to a file
    # file_name <- paste0("results/DEGs/Healthy_vs_", ref_dose, "_DEGs.tsv")
    # write.table(results_df, file = file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
  # **Second comparison: ref_dose vs each other dose**
  for (dose in doses) {
    if (dose == ref_dose || dose == "Healthy") {
      next  # Skip if dose is the same as ref_dose or if the dose is Healthy (already done)
    }
    
    print(paste("Processing comparison:", ref_dose, "vs", dose))
    
    # Subset metadata for the reference dose and the comparison dose
    ref_dose_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == ref_dose, ]
    comparison_metadata <- sample_info_UC_H[sample_info_UC_H$Condition == dose, ]
    
    # Combine metadata for the two doses
    combined_metadata <- rbind(ref_dose_metadata, comparison_metadata)
    
    # Subset the expression data for the samples in both doses
    dose_expression <- filtered_expression_data_UC_H[, combined_metadata$SampleID]
    
    # Set the current dose as the reference
    combined_metadata$Condition <- factor(combined_metadata$Condition, levels = c(ref_dose, dose))
    
    # Create the design matrix based on the combined metadata
    design <- model.matrix(~ combined_metadata$Condition)
    
    # Fit the linear model
    fit <- lmFit(dose_expression, design)
    
    # Apply empirical Bayes moderation
    fit <- eBayes(fit)
    
    # Extract results for the comparison (ref_dose vs dose)
    results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")
    
    # Store results in the results list
    result_name <- paste0(ref_dose, "_vs_", dose)
    results_list[[result_name]] <- results
    
    # Prepare results for output
    results_df <- as.data.frame(results)
    results_df$"Gene.symbol" <- rownames(results_df) 
    results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
    
    # # Save results to a file
    # file_name <- paste0("results/DEGs/", ref_dose, "_vs_", dose, "_DEGs.tsv")
    # write.table(results_df, file = file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


for (i in names(results_list)){
  results <- results_list[[i]]
  filt_results <- results %>%
    filter((logFC > 1 | logFC < -1) & adj.P.Val < 0.05)
  assign(i, filt_results)
}




  



# 
# perform_DE_analysis <- function(sample_info, expression_data, doses, condition) {
# 
#   results_list <- list()
#   
#   for (dose in doses) {
#     print(paste("Processing dose:", dose, "vs:", condition))
#     
#     # Filter metadata for dose and placebo samples
#     dose_metadata <- sample_info[sample_info$Condition == dose, ]
#     condition_metadata <- sample_info[sample_info$Condition == condition, ]
#     
#     # Combine metadata for dose and placebo
#     combined_metadata <- rbind(dose_metadata, condition_metadata)
#     
#     # Subset expression data for selected samples
#     dose_expression <- expression_data[, combined_metadata$SampleID]
#     
#     # Create design matrix (Placebo as reference)
#     design <- model.matrix(~ combined_metadata$Condition)
#     
#     # Fit the linear model
#     fit <- lmFit(dose_expression, design)
#     fit <- eBayes(fit)
#     
#     # Extract results for the dose effect
#     results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")
#     
#     # Store results in a named list
#     result_name <- paste0(dose, "_vs_",condition)
#     results_list[[result_name]] <- results
#     
#     # Save results as a file
#     results_df <- as.data.frame(results)
#     results_df$"Gene.symbol" <- rownames(results_df) 
#     results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
#     file_name <- paste0("results/DEGs/", dose, "_vs_",condition, "_DEGs.tsv")
#     write.table(results_df, file = file_name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   }
#   
#   return(results_list)
# }
# 
# results_DEGs_placebo <- perform_DE_analysis(filtered_sample_info,filtered_expression,doses,"Placebo")
# 





# 
# sample_info_new <- sample_info_rem %>% 
#   filter(Condition != "Healthy" & Week == "WEEK I-8")%>%
#   mutate(condition = ifelse(grepl("Ustekinumab", Condition), "UC", "Placebo"))
# 
# filt_data <- filter_samples(sample_info_new,filtered_expression_data_rem)
# 
# sample_info_r <- filt_data$filtered_sample_info
# exp_data_r <- filt_data$filtered_expression
# 
# design <- model.matrix(~ sample_info_r$Remission)
# colnames(design)
# 
# fit <- lmFit(exp_data_r, design)
# 
# fit <- eBayes(fit)
# 
# results <- topTable(fit, coef = 2, number = Inf, sort.by = "B")

# 
# results_df <- as.data.frame(results)
# results_df$"Gene.symbol" <- rownames(results_df) 
# results_df <- results_df[, c("Gene.symbol", setdiff(names(results_df), "Gene.symbol"))]
# file_name <- paste0("results/DEGs/UC_Placebo_Remission_Week8.tsv")
# write.table(results_df, file = file_name, sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)
# 
# file_p_value = paste0("results/DEGs/UC_Placebo_Remission_Week8_adj_p_value_0.05.tsv")
# file_p_value = paste0("results/DEGs/UC_Placebo_Remission_Week8_adj_p_value_0.05_logFC.tsv")
# results_sig_p_value <- results_df %>%
#   filter(adj.P.Val < 0.05 & (logFC > 1 | logFC < -1))
# write.table(results_sig_p_value,file = file_p_value,sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)

png(
  "~/Documents/Dissertation/results/UC_Placebo_remission_Week0_PCA_plot_average_4.png",
  width = 1024,
  height = 768
)

pca_data <- prcomp(t(exp_data_r))

condition_numeric <- as.numeric(factor(sample_info_r$Remission))

plot(
  pca_data$x[, 1],
  pca_data$x[, 2],
  col = condition_numeric,
  pch = 16,
  cex = 2,
  main = "PCA Plot",
  xlab = "PC1",
  ylab = "PC2"
)

condition_numeric

legend(
  "topright",
  legend = c("No remission", "Remission"),
  col = 1:2,
  pch = 16,
  title = "UC and Placebo PCA Week 0 - Remission"
)

dev.off()








result_files <- 
  c(
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0__vs_Ustekinumab_6_mg_kg__520_mg___I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0__vs_Ustekinumab_6_mg_kg__390_mg___I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0__vs_Ustekinumab_130_mg_IV__I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0__vs_Healthy_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0__vs_Placebo_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__260_mg___I_0_bl_week8_DEGs.tsv",
    "results/DEGs/Ustekinumab_130_mg_IV__I_0__vs_Healthy_DEGs.tsv",
    "results/DEGs/Ustekinumab_130_mg_IV__I_0__vs_Placebo_DEGs.tsv",
    "results/DEGs/Ustekinumab_130_mg_IV__I_0__vs_Ustekinumab_6_mg_kg__520_mg___I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_130_mg_IV__I_0__vs_Ustekinumab_6_mg_kg__390_mg___I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_130_mg_IV__I_0_bl_week8_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__390_mg___I_0_bl_week8_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__390_mg___I_0__vs_Ustekinumab_6_mg_kg__520_mg___I_0__DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__390_mg___I_0__vs_Healthy_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__520_mg___I_0_bl_week8_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__520_mg___I_0__vs_Placebo_DEGs.tsv",
    "results/DEGs/Ustekinumab_6_mg_kg__520_mg___I_0__vs_Healthy_DEGs.tsv",
    "results/DEGs/UC_Placebo_Remission_Week8_DEGs.tsv"
  )

# 1. Set your threshold values
logFC_threshold <- 1
pval_threshold <- 0.05

# 3. Create an empty list to store summary information
summary_list <- list()

# 4. Loop over each file to process them
for (file in result_files) {
  # 5. Extract the comparison name
  comparison_name <- gsub("_DEGs\\.tsv$", "", basename(file))
  
  # 6. Read the DEG results from the file
  results_df <- read.table(file, sep = "\t", header = TRUE)
  
  # 7. Apply filters based on logFC and adjusted p-value (adj.P.Val)
  significant_results <- results_df %>%
    filter((logFC > logFC_threshold | logFC < -logFC_threshold) & adj.P.Val < pval_threshold)
  
  # 8. Calculate the counts for the summary table
  total_DEGs <- nrow(results_df)  # Total number of DEGs (before filtering)
  significant_DEGs <- nrow(significant_results)  # Significant DEGs after filtering
  upregulated_DEGs <- sum(significant_results$logFC > logFC_threshold)  # Upregulated DEGs
  downregulated_DEGs <- sum(significant_results$logFC < -logFC_threshold)  # Downregulated DEGs
  
  # 9. Create a summary row for this comparison
  summary_row <- data.frame(
    Comparison = comparison_name,
    Total_DEGs = total_DEGs,
    Significant_DEGs = significant_DEGs,
    Upregulated_DEGs = upregulated_DEGs,
    Downregulated_DEGs = downregulated_DEGs
  )
  
  # 10. Add the summary row to the list
  summary_list[[comparison_name]] <- summary_row
}

# 11. Combine all summary rows into a single data frame
summary_df <- do.call(rbind, summary_list)

# 12. Save the final summary table to a file
write.table(summary_df, "results/DEGs/summary_DEGs_count_table.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# 13. View the summary table
head(summary_df)