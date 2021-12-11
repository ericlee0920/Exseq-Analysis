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
