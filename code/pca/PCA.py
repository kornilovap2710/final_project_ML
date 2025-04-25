import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# uniti_data_all = pd.read_csv('/Users/korni/Documents/Dissertation/data/UNITI_for_Polina/UNITI_for_Polina/clean_UNITI_biopsy_rma_intensities.csv')
# expression_data = uniti_data.iloc[:, 1:]

# # Calculate the mean and standard deviation for each gene
# means = expression_data.mean()
# std_devs = expression_data.std()

# print("Means:\n", means)
# print("\nStandard Deviations:\n", std_devs)
#
#
# scaler = StandardScaler()
# standardized_data = scaler.fit_transform(expression_data)
#
# # Create a DataFrame for the standardized data
# standardized_df = pd.DataFrame(standardized_data, columns=expression_data.columns)


# print("Standardized Data Means:\n", standardized_df.mean())
# print("\nStandardized Data Standard Deviations:\n", standardized_df.std())
pd.set_option('display.max_columns', None)
exp_data_UC_H = pd.read_csv("/Users/korni/Documents/Dissertation/results/UC_dose_130_without_placebo_vs_HC_expression_datafile.tsv", sep = '\t')
genes = pd.read_csv("/Users/korni/Documents/Dissertation/results/UC_dose_130_vs_H_DEGs.tsv", sep = '\t')
condition = pd.read_csv("/Users/korni/Documents/Dissertation/results/UC_dose_130_without_placebo_vs_HC_metadata.tsv", sep = '\t')

# # Filter the expression data to contain rows with DEGs only
# df = exp_data_UC_H[exp_data_UC_H["Gene.Symbol"].isin(genes["GeneID"].values)]
# df_transposed = df.melt(id_vars="Gene.Symbol",var_name="SampleID",value_name="Gene expression")
# merged_df = pd.merge(df_transposed,condition, on="SampleID")
# # print(merged_df.head())
#
# all_genes_transposed = exp_data_UC_H.melt(id_vars="Gene.Symbol",var_name="SampleID",value_name="Gene expression")
# merged_all_genes = pd.merge(all_genes_transposed,condition, on="SampleID")

df_t = exp_data_UC_H.set_index('Gene.Symbol').T.reset_index()
df_t.rename(columns={'index': 'SampleID'}, inplace=True)
print(df_t.head())
merged_df = pd.merge(df_t,condition, on="SampleID")
# print(merged_df.head())

y=merged_df["Condition"]
X = merged_df.drop(columns=["SampleID","Condition"])
print(X.var(axis=0))

# # Standardize the gene expression data
# scaler = StandardScaler()
# X_scaled = scaler.fit_transform(X)
# Perform PCA
pca = PCA(n_components=2)  # Adjust number of components if needed
principal_components = pca.fit_transform(X)

pca_df = pd.DataFrame(principal_components, columns=[f"PC{i+1}" for i in range(principal_components.shape[1])])
pca_df["Condition"] = y

print(pca_df.head())

explained_variance = pca.explained_variance_ratio_
print("Explained variance by each component:", explained_variance)


