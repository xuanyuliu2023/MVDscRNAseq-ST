conda activate hdWGCNA

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load the snRNA-seq dataset
seurat_obj <- get(load('/work/project/xuanyu/project/MVD/seurat/round2/VIC/subset_56339cells/seuObj.integrated_VIC_56339cells_subset.Rdata'))
DefaultAssay(seurat_obj) <- "RNA"

allGenes <- rownames(seurat_obj)
geneTypeInfo <- read.table(file='/work/project/xuanyu/resource/10XGenomics/geneSymbol_geneType.tsv',header=F,sep='\t',stringsAsFactors=F)
protein_coding_genes <- subset(geneTypeInfo,V2=="protein_coding")$V1
allProteinCodingGenes <- allGenes[allGenes %in% protein_coding_genes]

#only consider protein coding genes 
seurat_obj <- seurat_obj[allProteinCodingGenes,]

# Set up Seurat object for WGCNA
# Most of the information computed by hdWGCNA is stored in the Seurat object’s @misc slot.
# Selects the genes that will be used for WGCNA. The user can select genes using three different approaches using the gene_select paramete
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "hdWGCNA" # the name of the hdWGCNA experiment
)

#Construct metacells
#metacells are aggregates of small groups of similar cells originating from the same biological sample of origin. 
#The k-Nearest Neighbors (KNN) algorithm is used to identify groups of similar cells to aggregate, 
#and then the average or summed expression of these cells is computed,
# thus yielding a metacell gene expression matrix
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("subcluster", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter 20-75 # the larger the value of k, the fewer the obtained clusters.
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'subcluster' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

# Co-expression network analysis
# Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c('VIC1','VIC2','VIC3','VIC4','VIC5','VIC6','VIC7','VIC8'), # the name of the group of interest in the group.by column
  group.by='subcluster', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

#Select soft-power threshold
#This is an extremely important step in the scWGNCA pipleine
# Test different soft powers:

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # set this to FALSE since we did this above
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)
# Construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=7, # input soft power threshold 7 here based on PlotSoftPowers(
  setDatExpr=FALSE,
  detectCutHeight = 0.995,
  minModuleSize = 50, #minimal module size
  mergeCutHeight = 0.2
)
# visualize the WGCNA dendrogram，
# Each leaf on the dendrogram represents a single gene, 
# and the color at the bottom indicates the co-expression module assignment.
# Importantly, the “grey” module consists of genes that were not grouped into any co-expression module. 
# The grey module should be ignored for all downstream analysis and interpretation.

PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.copy2pdf(file='dendrogram.pdf')


plan("multisession", workers = 4)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- Seurat::ScaleData(
 seurat_obj,
 features = GetWGCNAGenes(seurat_obj),
 vars.to.regress = c('nCount_RNA', 'percent.mito')
)

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="orig.ident"
)

#The ME matrices are stored as a matrix where each row is a cell and each column is a module
# harmonized module eigengenes:
#Module eigengene is defined as the first principal component of the expression matrix of the corresponding module.
hMEs <- GetMEs(seurat_obj)
save(hMEs,file='hMEs.Rdata')

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

#Compute module connectivity
# “hub genes”, those which are highly connected within each module
# The conectivity is determined by the eigengene-based connectivity, also known as kME, of each gene
# compute intramodular connectivity:
seurat_obj <- ModuleConnectivity(seurat_obj) #slow

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)
pdf(file='KMEs.pdf')
p
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj)
write.table(modules,file='gene_modules_info.tsv',sep='\t',quote=F,row.names=F,col.names=T)

# get hub genes # set n_hubs=10000 to compute kME for all genes 
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100000)

#percentile rank
perc.rank <- function(x) trunc(rank(x))/length(x)
my.df <- data.frame(x=rnorm(200))
my.df <- within(my.df, xr <- perc.rank(x))

hub_df <- hub_df %>% group_by('module') %>% mutate(percRank=perc.rank(kME))
write.table(hub_df,file='hub_genes_modules_info.tsv',sep='\t',quote=F,row.names=F,col.names=T)


# compute gene scoring for the top 25 hub genes by kME for each module
# with AddModuleScore function
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)


# Basic Visualization

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs (harmonized module eigengenes)
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf("hME_UMAPplot.pdf")
wrap_plots(plot_list, ncol=4)
dev.off()

# get hMEs from seurat object
hMEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(hMEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)

mods <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12")
mods <- rev(c("M6","M1","M8","M11","M5","M4","M9","M2","M10","M12","M7","M3"))
# plot with Seurat's DotPlot function

DotPlot(seurat_obj, features=mods, group.by = 'subcluster',scale = T) + coord_flip() +
  RotatedAxis() + theme_bw() + scale_color_gradient2(high='red', mid='grey95', low='blue')
dev.copy2pdf(file="hME_score_Dotplot.pdf")

# violin plot
feature = 'M8'

VlnPlot(
  seurat_obj,
  features = feature,
  group.by = 'subcluster',
  pt.size = 0 # don't show actual data points
) + geom_boxplot(width=.25, fill='white') +
xlab('') + ylab('hME') + NoLegend() + theme_bw()
dev.copy2pdf(file=paste0(feature,'_hME_violinPlot.pdf'))

# Individual module network plots
library(igraph)
ModuleNetworkPlot(seurat_obj,n_hubs=25)

# hubgene network
#visualize the top 3 hub genes and 6 other genes per module
pdf(file='hubgene_network.pdf')
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=6,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

#Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)
write.table(umap_df,file='umap_df_hub_gene.tsv',sep='\t',quote=F,row.names=F,col.names=T)
 library(igraph)


plan("multisession", workers = 1)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)


ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=0 ,# how many hub genes to plot per module?
  keep_grey_edges=TRUE
)
dev.copy2pdf(file='Module_UMAP_Plot.pdf')

# save
saveRDS(seurat_obj, file='hdWGCNA_object.rds')

# load
seurat_obj <- readRDS(file='hdWGCNA_object.rds')

-----------
#-----------------------------------
#Differential module eigengene (DME) analysis
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# load
seurat_obj <- readRDS('hdWGCNA_object.rds')
# set random seed for reproducibility
set.seed(12345)
#DME analysis comparing two groups
group1_name <- 'FED'
group2_name <- 'CTRL'

group1 <- seurat_obj@meta.data %>% subset(group == group1_name) %>% rownames
group2 <- seurat_obj@meta.data %>% subset(group == group2_name) %>% rownames
DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='hdWGCNA'
)
DMEs

> DMEs$module
 [1] "M6"  "M7"  "M1"  "M8"  "M12" "M2"  "M3"  "M11" "M9"  "M4"  "M10" "M5"

> table(hub_df$module)

grey   M1   M2   M3   M4   M5   M6   M7   M8   M9  M10  M11  M12
   0  818  274  438  307   65  799  213   93   74   72  130   86

moduleSize.df <- data.frame(table(hub_df$module))
library(plyr)
DMEs$moduleSize <- mapvalues(rownames(DMEs),from=as.character(moduleSize.df$Var1),to=as.character(moduleSize.df$Freq))

write.table(DMEs,file=paste0(group1_name,"_vs_",group2_name,"_DMEs.tsv"),quote=F,sep='\t',row.names=T,col.names=T)

# vis

PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  mod_point_size =as.numeric(DMEs$moduleSize),
  wgcna_name = 'hdWGCNA',
  show_cutoff = TRUE
) + theme_bw()
ggsave(file=paste0(group1_name,"_vs_",group2_name,"_DMEsVolcano2.pdf"),units='in',width=6,height=4.81)

#plot with ggplot2
DEMs.df <- read.table(file=paste0(group1_name,"_vs_",group2_name,"_DMEs.tsv"),header=T,sep='\t')
ggplot(DEMs.df,aes(x=avg_log2FC,y=-log10(p_val_adj),size=moduleSize,color=module,label=module)) + geom_point() + 
theme_bw() + xlim(-1,1)
ggsave(file=paste0(group1_name,"_vs_",group2_name,"_DMEsVolcano.pdf"))
#One-versus-all DME analysis
group.by = 'subcluster'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'subcluster',
  wgcna_name = 'hdWGCNA'
)
write.table(DMEs,file="subcluster_DMEs.tsv",quote=F,sep='\t',row.names=T,col.names=T)
p <- PlotDMEsVolcano(
  seurat_obj,
  DMEs_all,
  wgcna_name = 'hdWGCNA',
  plot_labels=TRUE,
  show_cutoff=TRUE
)

# facet wrap by each cell type
pdf(file='subcluster_DMEs.pdf')
p + facet_wrap(~group, ncol=3)
dev.off()

#-------------------------
# Motif Analysis
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

# packages for TF motif analysis
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GeneOverlap)
library(ggseqlogo)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# load
seurat_obj <- readRDS('hdWGCNA_object.rds')

# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan with these settings for the human dataset
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86,
  wgcna_name = 'hdWGCNA'
)
dim(GetMotifMatrix(seurat_obj))


# TF target genes
target_genes <- GetMotifTargets(seurat_obj)

# check target genes for one TF:
head(target_genes$SOX9)

# overlap between modules & TF target genes:
seurat_obj<- OverlapModulesMotifs(seurat_obj)

# look at the overlap data
ModuleTFmotifOverlap.df<- GetMotifOverlap(seurat_obj)
write.table(ModuleTFmotifOverlap.df,file='ModuleTFmotifOverlap.tsv',sep='\t',col.names=T,row.names=F,quote=F)

# have bugs
MotifOverlapBarPlot(
  seurat_obj,
  module_names = 'M1',
  outdir = 'motifs/MotifOverlaps/',
  plot_size=c(5,6)
)

# Compute motif target gene expression score:
library(UCell)
seurat_obj <- MotifTargetScore(
  seurat_obj,
  method='UCell',
  maxRank=50000
)

#Plot overlap between motif of interest and modules:
df <- GetMotifOverlap(seurat_obj)

TF_name <- 'SOX9'


cur_df <- df %>% subset(tf == TF_name)
plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") + theme_bw() +
  ggtitle(paste0(TF_name," overlap")) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

pdf(paste0(TF_name,'_motif_overlap_or.pdf'))
p
dev.off()



# --------------Module Trait Correlation------------
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# load
seurat_obj <- readRDS('hdWGCNA_object.rds')
# convert variable
library(plyr)
#seurat_obj$msex <- as.factor(as.numeric(mapvalues(seurat_obj$group,from=c('female','male'),to=c(0,1))))
seurat_obj$mage <- as.numeric(seurat_obj$age)
seurat_obj$mBMI <- as.numeric(seurat_obj$BMI)
# list of traits to correlate
cur_traits <- c('nCount_RNA', 'nFeature_RNA', 'mage','msex')
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='seurat_clusters'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj)
#mt_cor is a list containing three items; cor which holds the correlation results,
#pval which holds the correlation p-values, and fdr which holds the FDR-corrected p-values.
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'black',
  high_color = 'firebrick3',
  mid_color = 'white',
  low_color = 'navy',
  plot_max = 0.2,
  combine=TRUE
)
