#mport pickle5
import scvelo as scv
scv.settings.set_figure_params('scvelo')
# load data
adata = scv.read('/work/project/xuanyu/project/MVD/velocyte/merged.loom')
adata.var_names_make_unique()
#Abundance of ['spliced', 'unspliced', 'ambiguous']: [0.8  0.15 0.05]


def convert_string_concise(original_string):
    # Split the string and remove 'x' using replace instead of rstrip
    identifier, sequence = original_string.split(':')
    sequence = sequence.replace('x', '')
    
    # Use string formatting to create the result
    result = f"{sequence}-1___{identifier}"
    
    return result

# Test the concise function
original = 'CTRL1:AAACGCTGTAACAAGTx'
converted = convert_string_concise(original)
print(converted)  # Output should be 'AAACGCTGTAACAAGT-1___CTRL1'

# change the index
CTRL1 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL1']
CTRL3 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL3']
CTRL4 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL4']
CTRL5 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL5']
CTRL6 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL6']
CTRL7 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL7']
CTRL9 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL9']
CTRL10 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'CTRL10']

FED2 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD2']
FED3 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD3']
FED4 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD4']
FED6 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD6']
FED7 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD7']
FED8 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD8']
FED9 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD9']
FED10 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD10']
FED14 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD14']
FED15 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD15']
FED16 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD16']
FED17 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD17']
FED22 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD22']
FED23 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD23']
FED24 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD24']

BD1 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'BD1']
BD2 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'BD2']
BD3 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'BD3']
BD4 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'BD4']
BD5 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'BD5']
BD6 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD11']
BD7 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD13']

RMVD1 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD1']
RMVD2 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV2']
RMVD3 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV3']
RMVD6 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV6']
RMVD7 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV7']
RMVD8 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV8']
RMVD10 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV10']
RMVD12 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD12']
RMVD14 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV14']
RMVD15 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'RMV15']
RMVD20 = [convert_string_concise(i) for i in adata.obs.index if i.split(':')[0] == 'DMVD20']
new_index = CTRL1 + CTRL3 + CTRL4 + CTRL5 + CTRL6 + CTRL7 + CTRL9 + CTRL10 + FED2 + FED3 + FED4 + FED6 + FED7 + FED8 + FED9 + FED10 + FED14 + FED15 +FED16 + FED17 + FED22 + FED23 + FED24 + BD1 + BD2 + BD3 + BD4 + BD5 + BD6 + BD7 + RMVD1 + RMVD2 + RMVD3 + RMVD6 + RMVD7 + RMVD8 + RMVD10 + RMVD12 + RMVD14 + RMVD15 + RMVD20

adata.obs = adata.obs.reindex(new_index) 
#save the adata
adata.write('adata_scvelo_raw407808cells.h5ad',compression='gzip')
# load adata
#adata = scv.read('adata_scvelo_raw.h5ad')

# load data; select gene and cells
import pickle5
csvDir='/work/project/xuanyu/project/MVD/scanpy/VIC_subset_56339cells'
file_adata_obs = open(csvDir+'/'+ 'adata.obs.obj', 'rb+')
adata_obs = pickle5.load(file_adata_obs)
file_adata_var = open(csvDir+'/'+ 'adata.var.obj', 'rb+')
adata_var = pickle5.load(file_adata_var)
file_adata_obsm = open(csvDir+'/'+ '/adata.obsm.obj', 'rb+')
adata_obsm = pickle5.load(file_adata_obsm)

intersect_index = adata.var.index.intersection(adata_var.index)
adata = adata[:,intersect_index]
adata = adata[adata_obs.index,:]
#vars  = adata.var.index.tolist()
adata.var.index.to_series().to_csv('intersection_vars.csv', header=False, index=False)
vars=adata.var.index
adata_var_convert = adata_var.loc[vars]


var_result = pd.concat([adata.var,adata_var_convert], axis=1)
# cleaning the data, i.e., throwing away everything not needed for velocity estimation
scv.utils.show_proportions(adata)
#Abundance of ['spliced', 'unspliced']: [0.84 0.16]

scv.utils.cleanup(adata, clean='all')
adata.obs = adata_obs
adata.var = var_result
adata.obsm = adata_obsm


#load the normalized expression data matrix
import pandas as pd
normalized_expr = pd.read_csv('/work/project/xuanyu/project/MVD/seurat/round2/VIC/subset_56339cells/seurat2scanpy/normalize_expr.data_intersect_genes.csv', sep=',', header=None, names=None, index_col=None)

import scipy.sparse as sp
# Convert the numpy array to a CSR matrix
adata_normalized_expr_csr_matrix = sp.csr_matrix(normalized_expr)
adata.layers['normalized_expr'] = adata_normalized_expr_csr_matrix

#visualize the expression of a gene
sc.pl.embedding(adata, basis='draw_graph_fa', color='magma', gene_symbols='MYH6', use_raw=None, layer='normalized_expr')

#subset
import pandas as pd
theSubsetCellId = pd.read_csv('/work/project/xuanyu/project/DMVD/subclustering/VIC/theSubsetCellId.tsv', sep=',', header=0, names=None, index_col=0) 
adata = adata[theSubsetCellId.index,:]

# Normalized X and spliced/unspliced count data；Don't worry if X is already processed; it will detect that automatically and not touch X in that case
scv.pp.normalize_per_cell(adata)

#The first and second order moments (basically mean and uncentered variance) are computed among nearest neighbors in PCA space
scv.pp.moments(adata, n_pcs=30, n_neighbors=30,use_rep="X_pca")


del adata.obsm['X_diffmap']
# run the dynamical model
scv.tl.recover_dynamics(adata,n_jobs=60)
scv.tl.velocity(adata, mode='dynamical')

# try different modes
scv.tl.velocity(adata, mode='steady_state')
scv.tl.velocity(adata, mode='stochastic')

# velocity_graph
scv.tl.velocity_graph(adata,n_jobs=1)


adata.uns['subcluster_colors'] = ["#5B74F4","#E82C61","#ff7f0e","#2ca02c","#7733B7","#8c565B","#E236AF","#8C8C00"]
# UMAP stream Plot
scv.pl.velocity_embedding_stream(adata, basis='umap',color='subcluster',min_mass=0.1)
#scvelo_UMAP_stream.pdf/svg

# UMAP grid Plot
scv.pl.velocity_embedding_grid(adata, basis='umap',color='subcluster',min_mass=0.1)
#scvelo_UMAP_grid.pdf

# FA
scv.pl.velocity_embedding_stream(adata, basis='draw_graph_fa',color='subcluster',min_mass=0.1, legend_loc ='right margin')
#scvelo_FA_stream.pdf

#Velocities in cycling progenitors
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

#Two more useful stats: - The speed or rate of differentiation is given by the length of the velocity vector. 
#The coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities) provides a measure of confidence
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

df = adata.obs.groupby('subcluster')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
>>> df
subcluster               VIC1      VIC2      VIC3      VIC4      VIC5      VIC6      VIC7      VIC8
velocity_length      5.278982  4.545827  4.897759  4.632742  4.857658  4.439327  4.553288  4.017638
velocity_confidence  0.701841  0.634118  0.651703  0.645642  0.645379  0.635073  0.627719  0.659393


scv.pl.velocity_graph(adata, threshold=.1)
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

# velocity pseudotime 
#Contrarily to diffusion pseudotime, it implicitly infers the root cells 
#and is based on the directed velocity graph instead of the similarity-based diffusion kernel.
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', size=80,basis='umap')

#run PAGA velocity-inferred directionality.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='subcluster')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(adata, basis='umap', size=50, alpha=.05, min_edge_width=3, node_size_scale=1.5)

# Latent time 
# only for dynamic model
scv.tl.latent_time(adata)

# UMAP layout
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,basis='umap')
#UMAP_latent_time.pdf

# Fa layout
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,basis='draw_graph_fa')
#FA_latent_time.pdf

# save
adata.write('adata_scvelo.h5ad',compression='gzip')

# Top-likelihood genes
#scVelo computes a likelihood for each gene and cell for a model-optimal latent time
#Driver genes display pronounced dynamic behavior and
#  are systematically detected via their characterization by high likelihoods in the dynamic model.
all_ranked_genes_series = adata.var['fit_likelihood'].sort_values(ascending=False).dropna()
all_ranked_genes_series.to_csv('all_ranked_genes.tsv',sep='\t',index=True)
all_ranked_genes = all_ranked_genes_series.index

# all_ranked genes heatmap
scv.pl.heatmap(adata, var_names=all_ranked_genes, sortby='latent_time', col_color='seurat_clusters',
n_convolve=100,color_map="viridis",row_cluster=False,font_scale=0.1,yticklabels=True)
#Heatmap_allrankedGenes.pdf

matplotlib_axis_heatmap = scv.pl.heatmap(adata, var_names=all_ranked_genes, sortby='latent_time', col_color='seurat_clusters',
n_convolve=100,color_map="viridis",row_cluster=False,font_scale=0.1,yticklabels=True,show=False)

# show ordered genes of the heatmap
list(matplotlib_axis_heatmap.data.index)
matplotlib_axis_heatmap.data.to_csv("allrankedGeneHeatmapData.tsv",sep='\t',header=True,index=True)


# top 100 genes heatmap
top_driver_genes = all_ranked_genes[:100]
scv.pl.heatmap(adata, var_names=top_driver_genes, sortby='latent_time', col_color='seurat_clusters',
n_convolve=100,color_map="viridis",row_cluster=False,font_scale=0.6,yticklabels=True)
#Heatmap_top100.pdf

# draw heatmap given a gene list
gene_list = ['DKK1','CXCL12','CXCR4','PDGFRB','CDH5','SEMA3G','ETV2', 'FLI1', 'SOX7', 'SOX17', 'TAL1','ERG']
scv.pl.heatmap(adata, var_names=gene_list, sortby='latent_time', col_color='seurat_clusters',
n_convolve=100,color_map="viridis",row_cluster=False,font_scale=0.6,yticklabels=True)
#candidate_gene_heatmap.pdf

# unspliced/spliced phase portrait plot given a gene list
var_names = ['DNASE1L3', 'LEPR', 'MIA','TIMP1','FGF7','CCDC3']
scv.pl.scatter(adata, var_names, frameon=True,color='seurat_clusters')
#candidate_gene_phase_portrait.pdf

# Gene expression dynamics along latent time for the driver genes given a gene list
scv.pl.scatter(adata, x='latent_time', y=var_names, frameon=True,color='seurat_clusters')
scv.pl.scatter(adata, x='latent_time', y=var_names, frameon=True,color='seurat_clusters',
add_polyfit='grey',linewidth=5,lowess=True)
# candidate_gene_expression_dynamics.pdf

#Cluster-specific top-likelihood genes
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.to_csv('cluster_specific_top_likelihood_genes.tsv',sep='\t')
df.head(5)
for cluster in ['EC0', 'EC1']:
    scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, frameon=False,color='seurat_clusters')

scv.pl.scatter(adata, x='latent_time', y=list(df['EC0'][:5]), ylabel='EC0', frameon=True, color='seurat_clusters')
