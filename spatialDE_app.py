import pandas as pd
import NaiveDE
import SpatialDE


def run(counts_csv, metadata_csv, batch_correct=False):
    # read files, note that counts are gene x cell (we transpose here)
    meta_cells = pd.read_csv(metadata_csv)
    counts = pd.read_csv(counts_csv, index_col=0).T

    # use Anscombes approximation to variance stabilize Negative Binomial data
    norm_expr = NaiveDE.stabilize(counts.T).T

    # limma's removeBatchEffect function; when only one batch, comment out this line
    if batch_correct:
        norm_expr = NaiveDE.regress_out(meta_cells, norm_expr.T, 'np.log(total_counts)').T

    # run SpatialDE requires:
    #       X - spatial locations
    #       resid_expr - normalized counts
    # result will be a DataFrame with P-values and other relevant values for each gene:
    #       g - The name of the gene
    #       pval - The P-value for spatial differential expression
    #       qval - Significance after correcting for multiple testing
    #       l - A parameter indicating the distance scale a gene changes expression over

    X = meta_cells[['center_x', 'center_y']]
    results = SpatialDE.run(X, norm_expr)
    results = results.sort_values('qval')

    return results
