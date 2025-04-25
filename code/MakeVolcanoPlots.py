from adjustText import adjust_text
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# original_input = "files/Marton_TNF/Ustekinumab_520mg_vs_placebo_T8/Ustekinumab_520mg_vs_placebo_T8_DEAresults.tsv"
original_input = "../results/DEGs/UC_Placebo_Remission_Week0.tsv"
# original_input = "files/Domenico_results/DE_res_UC_vs_HC.tsv"
checking_array = []

df = pd.read_csv(original_input, sep = '\t')
df['-log10(p-value)'] = -np.log10(df['adj.P.Val'])

top_genes = df.nlargest(50, 'logFC')
bottom_genes = df.nsmallest(50, 'logFC')

plt.figure(figsize = (20, 16))
plt.scatter(df['logFC'], df['-log10(p-value)'], c = 'gray', alpha=0.5)
#
# significant_genes = df[df['adj.P.Val'] < 0.05]
# plt.scatter(significant_genes['logFC'], significant_genes['-log10(p-value)'], c = 'red', alpha = 0.5)

significant_genes = df[(df['adj.P.Val'] < 0.05) & ((df['logFC'] > 1) | (df['logFC'] < -1))]
plt.scatter(significant_genes['logFC'], significant_genes['-log10(p-value)'], c='red', alpha=0.5)

texts = []
for i, row in top_genes.iterrows():
    if np.isfinite(row['logFC']) and np.isfinite(row['-log10(p-value)']):
        text = plt.text(row['logFC'], row['-log10(p-value)'], row['Gene.symbol'], fontsize=20, color='darkblue')
        texts.append(text)

# Adjust text labels to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle = '-', color = 'black'), force_text = 0.5, expand_points = (1.2, 1.2))

texts = []
for i, row in bottom_genes.iterrows():
    if np.isfinite(row['logFC']) and np.isfinite(row['-log10(p-value)']):
        text = plt.text(row['logFC'], row['-log10(p-value)'], row['Gene.symbol'], fontsize=20, color='darkblue')
        texts.append(text)

# Adjust text labels to avoid overlap
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black'), force_text=0.5, expand_points=(1.2, 1.2))

plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(ajd. p value)')
plt.title(f'Volcano plot for Remission vs No remission at Week 8 (Unifi - with placebo samples)')

plt.axvline(x=1, color='blue', linestyle='-', label='LogFC = 1')
plt.axvline(x=-1, color='blue', linestyle='-', label='LogFC = -1')

plt.axhline(y = -np.log10(0.05), color = 'blue', linestyle = '-', label = 'Significance Threshold (p=0.05)')

with open("../results/UC_Placebo_Remission_Week8_sign_genes.tsv", 'w') as out:

    out.write(f"GeneID\tlogFC\t-log10(p-value)\n")

    for i, row in significant_genes.iterrows():
        out.write(row['Gene.symbol'] + '\t' + format(row['logFC']) + '\t' + format(row['-log10(p-value)']) + '\n')        

plt.legend()
plt.savefig("../results/plots/UC_Placebo_Remission_Week8_volcano_plot.png")
plt.show()
