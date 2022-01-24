import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors


# csv locations
meta_cells_csv = "raw_counts_SA501_exseq/exseq_SA501_meta_cells.csv"
counts_csv = "raw_counts_SA501_exseq/exseq_SA501_rawcounts_genes.csv"

# read files, note that counts are gene x cell (we transpose here)
meta_cells = pd.read_csv(meta_cells_csv)
counts = pd.read_csv(counts_csv, index_col=0).T

# add extra column to match with count's indices
index = list('C' + meta_cells['cell_id'].astype(str))
meta_cells = meta_cells.assign(cell_index=index)
counts = counts.reindex(index)
pos = meta_cells[['center_x', 'center_y']]
pos.index = index
num_umi = counts.T.sum(axis=0)

# Create the Hotspot object and the neighborhood graph, calculate spatial variation
# NOTE: increasing n_neighbors will have affect on FDR
hs = hotspot.Hotspot(counts.T, model='danb', latent=pos)
hs.create_knn_graph(weighted_graph=False, n_neighbors=200,)
hs_results = hs.compute_autocorrelations(jobs=1)

# Select the genes with significant spatial autocorrelation
hs_genes = hs_results.index[hs_results.FDR < 0.05]

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=1)

# Create gene modules and plot module correlations
# NOTE: increasing min_gene_threshold will have affect on Modules
modules = hs.create_modules(min_gene_threshold=8, core_only=False, fdr_threshold=0.05)
modules.value_counts()
hs.results.join(hs.modules).to_csv("hotspot_results.csv")
hs.plot_local_correlations()
plt.savefig("hotspot_local_correlations.png")
plt.clf()

# Plot the module scores on top of positions, sort by Z score
module = 1
results = hs.results.join(hs.modules)
results = results.loc[results.Module == module]
genes = results.sort_values('Z', ascending=False).head(6).index

fig, axs = plt.subplots(2, 3, figsize=(11, 7.5))
for ax, gene in zip(axs.ravel(), genes):
    # NOTE: normalize data here
    expression = np.log2(hs.counts.loc[gene]/hs.umi_counts + 1)
    plt.sca(ax)
    sc = plt.scatter(x=hs.latent.iloc[:, 0], y=hs.latent.iloc[:, 1], s=2, c=expression,
                    edgecolors='none', cmap=matplotlib.cm.get_cmap('viridis'))
    plt.colorbar(sc)
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.title(gene)
plt.savefig(f"hotspot_modulescores_{module}.png")
plt.clf()

# Calculate module scores
module_scores = hs.calculate_module_scores()
module_scores.to_csv("hotspot_module_scores.csv")

# Plot the module scores on top of positions
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
for ax, module in zip(axs.ravel(), range(1, hs.modules.max()+1)):
    scores = hs.module_scores[module]
    plt.sca(ax)
    sc = plt.scatter(x=hs.latent.iloc[:, 0], y=hs.latent.iloc[:, 1], s=2, c=scores,
                    vmin=np.percentile(scores, 1), vmax=np.percentile(scores, 99),
                    edgecolors='none')
    plt.colorbar(sc)
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.xticks([])
    plt.yticks([])
    plt.title('Module {}'.format(module))
plt.savefig("hotspot_modules.png")
plt.clf()
