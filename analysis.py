import pandas as pd
import NaiveDE
import SpatialDE

# csv locations
meta_genes_csv = "raw_counts_SA501_exseq/exseq_SA501_meta_genes.csv"
meta_cells_csv = "raw_counts_SA501_exseq/exseq_SA501_meta_cells.csv"
counts_csv = "raw_counts_SA501_exseq/exseq_SA501_rawcounts_genes.csv"

# read files, note that counts are gene x cell (we transpose here)
meta_genes = pd.read_csv(meta_genes_csv)
meta_cells = pd.read_csv(meta_cells_csv)
counts = pd.read_csv(counts_csv, index_col=0).T

# add extra column to match with count's indices
index = list('C' + meta_cells['cell_id'].astype(str))
meta_cells = meta_cells.assign(cell_index=index)
counts = counts.reindex(index)

# use Anscombes approximation to variance stabilize Negative Binomial data
norm_expr = NaiveDE.stabilize(counts.T).T

# limma's removeBatchEffect function
resid_expr = NaiveDE.regress_out(meta_cells, norm_expr.T, 'np.log(total_counts)').T

# run SpatialDE requires:
#       X - spatial locations
#       resid_expr - normalized counts
# result will be a DataFrame with P-values and other relevant values for each gene:
#       g - The name of the gene
#       pval - The P-value for spatial differential expression
#       qval - Significance after correcting for multiple testing
#       l - A parameter indicating the distance scale a gene changes expression over
X = meta_cells[['center_x', 'center_y']]
results = SpatialDE.run(X, resid_expr)
results = results.sort_values('qval')
results.to_csv("results.csv", index=False)

#########################################################################

import numpy as np
import matplotlib.pyplot as plt

# visualize the data
plt.Figure(figsize=(6, 4))
plt.scatter(meta_cells['center_x'], meta_cells['center_y'], c='k')
plt.axis('equal')
plt.savefig("tissue.png")
plt.clf()

# get differentially expressed genes
alpha = 0.05
de_results = results.query('qval < 0.5')
de_genes = list(de_results.g)

# visualize differentially expressed genes, plot spatial location and color by expression level
show = len(de_genes)
plt.Figure(figsize=(10, 10))
for i, g in enumerate(de_genes[:show]):
    # plt.subplot(1, 3, i + 1)
    plt.scatter(meta_cells['center_x'], meta_cells['center_y'], c=norm_expr[g])
    plt.title(g)
    # plt.axis('equal')
    plt.colorbar(ticks=[])
    plt.savefig(f"tissue_{g}.png")
    plt.clf()

# visualize fraction of variance explained by spatial variation
plt.Figure(figsize=(10, 10))
plt.yscale('log')
plt.scatter(results['FSV'], results['qval'], c='black')
plt.axhline(0.05, c='black', lw=1, ls='--')
plt.gca().invert_yaxis()
plt.xlabel('Fraction spatial variance')
plt.ylabel('Adj. P-value')
plt.savefig("fsv.png")
plt.clf()

# run automatic expression histology (AEH) requires:
#       C: number of patterns
#       l: characteristic length scale for histological patterns.
# result will be two DataFrames with P-values and other relevant values for each gene:
#       histology_results - pattern membership information for each gene
#       patterns - realizations for the underlying expression for each histological pattern

idxs = np.array(de_results['l'].value_counts().index)
vals = np.array(de_results['l'].value_counts().values)
l = np.sum(idxs * vals) / np.sum(vals)
C = 3
histology_results, patterns = SpatialDE.aeh.spatial_patterns(X, resid_expr, de_results, C=C, l=l, verbosity=1)
histology_results.to_csv("histology.csv", index=False)
patterns.to_csv("patterns.csv", index=False)

# visualize expression in the tissue context
plt.Figure(figsize=(10, 10))
for i in range(3):
    # plt.subplot(1, 3, i + 1)
    plt.scatter(meta_cells['center_x'], meta_cells['center_y'], c=patterns[i])
    # plt.axis('equal')
    plt.title('Pattern {} - {} genes'.format(i, histology_results.query('pattern == @i').shape[0]))
    plt.colorbar(ticks=[])
    plt.savefig(f"pattern_{i}.png")
    plt.clf()
