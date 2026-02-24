


import numpy as np
import scanpy as sc
import cellrank as cr

import scvelo as scv

scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2

#文字以字体形式保存
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']= 42
mpl.rcParams['savefig.dpi'] = 600

# load adata
adata = scv.read('/work/project/xuanyu/project/MVD/scvelo/VIC/adata_scvelo.h5ad')
import pandas as pd
df_pseudotime = pd.read_csv('monocle3_pseudotime.tsv', sep=',')
adata.obs['Pseudotime'] = df_pseudotime['Pseudotime']
sc.pl.violin(adata, keys='dpt_pseudotime', jitter=0.4, size=2, cut=0,groupby='subcluster')

scv.pl.proportions(adata)
adata

#-----------------------------------precalculated skip ------------------ 
#Preprocess the data
scv.pp.filter_and_normalize(
    adata, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False
)

sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

#Run scVelo
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")

#--------------------------------------------------------------------------


## --------- Combine RNA velocity with expression similarity in high dimensions
#Combine with gene expression similarity
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()


#use CellRank to compute a transition matrix using RNA velocity and gene expression similarity and how it can be visualized in low dimensions
#Set up the VelocityKernel
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix() #如果报错 使用FatNode服务器

# use Palantir pseudotime to compute a directed cell-cell transition matrix
pk = cr.kernels.PseudotimeKernel(adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix()


combined_kernel = 0.4 * vk + 0.4 * pk + 0.2 * ck
print(combined_kernel)

#Write results to file
adata.write("velocity_connectivity_kernel.h5ad", compression="gzip")


sc.pl.embedding(adata, basis="DM_EigenVectors_multiscaled", color="subcluster")

# Identify initial and terminal states
g = cr.estimators.GPCCA(combined_kernel)
g.fit(cluster_key="subcluster", n_states=[8, 20],n_cells=50) # n_states: Number of macrostates to compute.
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100, basis="DM_EigenVectors_multiscaled" )


#dentify terminal macrostates.
g.predict_terminal_states(allow_overlap=True)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, basis="draw_graph_fa")
g.plot_macrostates(which="terminal", discrete=False,basis="draw_graph_fa")
#Each cell is colored according to the terminal state it most likely belongs to; higher color intensity reflects greater confidence in the assignment.

#identify initial states 
g.predict_initial_states(allow_overlap=True)
g.plot_macrostates(which="initial", legend_loc="right", s=100,basis="DM_EigenVectors_multiscaled")


sc.pl.embedding(
    adata,
    basis="umap",
    color=["Ins1", "Ins2", "Gcg", "Ghrl", "Sox9", "Anxa2", "Bicc1"],
    size=50,
)

# set inital & terminal states manually
g.set_initial_states(states=["VIC3"],allow_overlap=True)
g.set_terminal_states(states=["VIC2", "VIC1", "VIC6_1", "VIC5_1"])
g.plot_macrostates(which="terminal", legend_loc="right", size=100, basis="DM_EigenVectors_multiscaled")


#Estimating Fate Probabilities and Driver Genes

#
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False,cmap='magma',size=30,basis="DM_EigenVectors_multiscaled")

# color intensity reflecting the degree of lineage-bias
cr.pl.circular_projection(adata, keys=["subcluster"], legend_loc="right",dpi=600,save='MS_circular_projection.jpg')


#-----------------------------------------------
#Uncover driver genes
#Correlate fate probabilities with gene expression
#Normally restrict the correlation-computation to the relevant clusters.
driver_clusters = ["VIC3","VIC2"]

delta_df = g.compute_lineage_drivers(
    lineages=["VIC2"], cluster_key="macrostates_fwd", clusters=driver_clusters,layer='MAGIC_imputed_data',use_raw=False
)
delta_df.head(10)

#Visualize putative driver genes
adata.obs["fate_probabilities_C1_2"] = g.fate_probabilities["C1_2"].X.flatten()

# show fate probabilities towards the population and some of the top-correlating genes
sc.pl.embedding(
    adata,
    basis="draw_graph_fa",
    color=["fate_probabilities_C1_2"] + list(delta_df.index[:8]),
    color_map="viridis",
    s=50,
    ncols=3,
    vmax="p96",
    layer='normalized_expr',use_raw=False
)
#------------------------------------------------------------------------------------------------

# compute driver genes for all trajectories,
#driver_df = g.compute_lineage_drivers()
# Save the DataFrame to a TSV file
#driver_df.to_csv('driver_genes_for_all_trajectories.tsv', sep='\t', index=True)

# define set of genes to annotate
C7_genes = ["MYH3", "MET", "RGS7"]
C1_2_genes = ["MYH6", "MYH7", "FGF12"]

genes_oi = {
    "C7": C7_genes,
    "C1_2": C1_2_genes,
}

# make sure all of these exist in AnnData
assert [
    gene in adata.var_names for genes in genes_oi.values() for gene in genes
], "Did not find all genes"

# compute mean gene expression across all cells
adata.var["mean expression"] = adata.layers['normalized_expr'].toarray().mean(axis=0)

# visualize in a scatter plot
g.plot_lineage_drivers_correlation(
    lineage_x="C7",
    lineage_y="C1_2",
    adjust_text=True,
    gene_sets=genes_oi,
    color="mean expression",
    legend_loc="none",
    figsize=(5, 5),
    dpi=150,
    fontsize=9,
    size=50,
)



sc.pl.violin(adata, keys=["dpt_pseudotime"], groupby="subcluster", rotation=90)

#----------------------------Visualize expression trends via line plots
#The first step in gene trend plotting is to initialite a model for GAM fitting.

model = cr.models.GAMR(adata, n_knots=5, smoothing_penalty=10)

#如果报错 conda install conda-forge::r-mgcv
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="normalized_expr",
    genes=["MYH3", "MYH6","MYH11","PDGFRB","POSTN"],
    same_plot=True,
    ncols=5,
    time_key="dpt_pseudotime",
    hide_cells=True,
    weight_threshold=(1e-3, 1e-3),
)

geneTypeDict={}
with open('geneSymbol_geneType.tsv','r') as FILE:
    for line in FILE:
        lineList = line.strip().split('\t')
        gene = lineList[0]
        geneTypeDict[gene]= lineList[1]

TFList=[]
with open('TcoF-DB_humamTF.tsv','r') as FILE:
    for line in FILE:
        lineList = line.strip().split('\t')
        gene = lineList[0]
        TFList.append(gene)


lineageName='C16'

#Correlates gene expression with lineage probabilities, for a given lineage and set of clusters. Often, it makes sense to restrict this to a set of clusters which are relevant for the specified lineages.
drivers_df = g.compute_lineage_drivers(lineages=lineageName,
                                       clusters=None,
                                       cluster_key="macrostates_fwd",
                                       layer='normalized_expr',
                                       use_raw=False)

geneTypeList=[]
for gene in drivers_df.index:
    if gene in geneTypeDict:
        geneTypeList.append(geneTypeDict[gene])
    else:
        geneTypeList.append('NA')

drivers_df= drivers_df.assign(geneType=geneTypeList)

TFinfoList=[]
for gene in drivers_df.index:
    if gene in TFList:
        TFinfoList.append('TRUE')
    else:
        TFinfoList.append('FALSE')

drivers_df= drivers_df.assign(TF=TFinfoList)



sorted_drivers_df = drivers_df.sort_values(axis=0, by= lineageName+'_qval',  ascending=True)
# Save the DataFrame to a TSV file
sorted_drivers_df.to_csv('lineage_'+lineageName+'_driver_genes_for_the_trajectories_raw.tsv', sep='\t', index=True)


# select the top 2000 protein coding genes
subset_drivers_df = sorted_drivers_df[(sorted_drivers_df['geneType'] == 'protein_coding')].head(2000)
subset_drivers_df= subset_drivers_df.sort_values(axis=0, by= lineageName+'_corr',  ascending=False)

# Save the DataFrame to a TSV file
subset_drivers_df.to_csv('lineage_'+lineageName+'_driver_genes_for_the_trajectories_top2000ProtGenes.tsv', sep='\t', index=True)

# select the top 50 TF
subset_drivers_df_TF = sorted_drivers_df[(sorted_drivers_df['geneType'] == 'protein_coding') & (sorted_drivers_df['TF'] == 'TRUE')].head(50)
subset_drivers_df_TF.to_csv('lineage_'+lineageName+'_driver_genes_for_the_trajectories_top50TFs.tsv', sep='\t', index=True)



# read
subset_drivers_df = pd.read_csv('lineage_VIC2_1_driver_genes_for_the_trajectories_raw.tsv', sep='\t',index_col=0)

# plot heatmap for the top 2000 protein coding genes
genes2plot=subset_drivers_df.index
lineageName='VIC2_1'
cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages=lineageName,
    cluster_key="macrostates_fwd",
    show_fate_probabilities=True,
    data_key="normalized_expr",
    genes=genes2plot,
    time_key="Pseudotime",
    figsize=(12, 10),
    show_all_genes=False,
    weight_threshold=(1e-3, 1e-3),
    n_jobs=1 #MemoryError for multicores
)

'C1_heatmap_top2000ProtGenes.pdf'

# plot heatmap for the top 50 TFs
genes2plot=subset_drivers_df_TF.index

cr.pl.heatmap(
    adata,
    model=model,  # use the model from before
    lineages=lineageName,
    cluster_key="macrostates_fwd",
    show_fate_probabilities=True,
    data_key="normalized_expr",
    genes=genes2plot,
    time_key="dpt_pseudotime",
    figsize=(12, 10),
    show_all_genes=True,
    weight_threshold=(1e-3, 1e-3),
    n_jobs=1 #MemoryError for multicores
)

'C1_heatmap_top50ProtGenes.pdf'


#Cluster gene expression trends
# delete adata.uns[lineageName_trend] if it exists
lineageName_trend = "lineage_"+ lineageName + "_trend"
del adata.uns[lineageName_trend]

genes2cluster=subset_drivers_df.index

cr.pl.cluster_trends(
    adata,
    model=model,  # use the model from before
    lineage=lineageName,
    data_key="normalized_expr",
    genes=genes2cluster,
    time_key="dpt_pseudotime",
    weight_threshold=(1e-3, 1e-3),
    n_jobs=1,
    random_state=0,
    clustering_kwargs={"resolution": 0.4, "random_state": 0},
    neighbors_kwargs={"random_state": 0},
)
'C7_cluster_trends.pdf'

# obtain some genes from the same cluster, sort by mean expression
lineageName_trend = "lineage_"+ lineageName + "_trend"
gdata = adata.uns[lineageName_trend].copy()
gdata.obs.to_csv('lineage_'+lineageName+'_cluster_trends.tsv', sep='\t', index=True)



gdata
cols = ["mean expression"]
gdata.obs = gdata.obs.merge(
    right=adata.var[cols], how="left", left_index=True, right_index=True
)
gdata

mid_peak_genes = (
    gdata[gdata.obs["clusters"] == "0"]
    .obs.sort_values("mean expression", ascending=False)
    .head(8)
    .index
)

# plot 
cr.pl.gene_trends(
    adata,
    model=model,
    lineages=lineageName,
    cell_color="subcluster",
    data_key="normalized_expr",
    genes=list(mid_peak_genes),
    same_plot=True,
    ncols=3,
    time_key="dpt_pseudotime",
    hide_cells=True,
    weight_threshold=(1e-3, 1e-3),
)



# save the adata
adata.write('adata_cellrank.h5ad',compression='gzip')






