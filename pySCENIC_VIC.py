
#conda activate pyscenic

library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# expr_matrix.R
library(Seurat)
load('/work/project/xuanyu/project/MVD/seurat/round2/ALL/subset_78816cells/seuObj.integrated_MVD_ALL_subset_78816cells.Rdata')
expr_matrix <- as_matrix(seuObj.integrated@assays$RNA@counts)
dim(expr_matrix)
#cells <- row.names(subset(seuObj.integrated@meta.data, cellType =="CM"))
#expr_matrix <- expr_matrix[,cells]
#dim(expr_matrix)

allGenes <- rownames(seuObj.integrated)
geneTypeInfo <- read.table(file='/work/project/xuanyu/resource/10XGenomics/geneSymbol_geneType.tsv',header=F,sep='\t',stringsAsFactors=F)
protein_coding_genes <- subset(geneTypeInfo,V2=="protein_coding")$V1
allProteinCodingGenes <- allGenes[allGenes %in% protein_coding_genes]

#only consider protein coding genes 
expr_matrix <- expr_matrix[allProteinCodingGenes,]


header <- c('cell_id',colnames(expr_matrix))

write.table(paste(header, collapse='\t'),file='expr_matrix_RNA_counts.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(expr_matrix,file='expr_matrix_RNA_counts.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=TRUE)

 VIC1  VIC2  VIC3  VIC4  VIC5  VIC6  VIC7  VIC8
12833 11701  9186  7307  7034  6984  1095   199

VIC1.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC1')$cellID
VIC2.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC2')$cellID
VIC3.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC3')$cellID
VIC4.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC4')$cellID
VIC5.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC5')$cellID
VIC6.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC6')$cellID
VIC7.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC7')$cellID
VIC8.cells <- subset(seuObj.integrated@meta.data,subcluster=='VIC8')$cellID
VIC1.cells.sample <- VIC1.cells[sample(1:length(VIC1.cells),12833)]
VIC2.cells.sample <- VIC2.cells[sample(1:length(VIC2.cells),11701)]
VIC3.cells.sample <- VIC3.cells[sample(1:length(VIC3.cells),9186)]
VIC4.cells.sample <- VIC4.cells[sample(1:length(VIC4.cells),7307)]
VIC5.cells.sample <- VIC5.cells[sample(1:length(VIC5.cells),7034)]
VIC6.cells.sample <- VIC6.cells[sample(1:length(VIC6.cells),6984)]
VIC7.cells.sample <- VIC7.cells[sample(1:length(VIC7.cells),1095)]
VIC8.cells.sample <- VIC8.cells[sample(1:length(VIC8.cells),199)]
cells <- c(VIC1.cells.sample,VIC2.cells.sample,VIC3.cells.sample,VIC4.cells.sample,VIC5.cells.sample,VIC6.cells.sample,VIC7.cells.sample,VIC8.cells.sample)

expr_matrix <- expr_matrix[,cells]
dim(expr_matrix)
header <- c('cell_id',colnames(expr_matrix))
write.table(paste(header, collapse='\t'),file='expr_matrix_RNA_counts_56339cells.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(expr_matrix,file='expr_matrix_RNA_counts_56339cells.tsv',sep='\t',row.names=T,col.names=F,quote=F,append=TRUE)

# -------------step1_prep.sh -------------
Rscript ../expr_matrix.R /work/project/xuanyu/project/PEAT_latest/subclustering/Macrophage/seuObj.integrated.Macrophage.RData


#-----------------step2_pyscenic_run.sh 
#--------human--------------------------------------------
pyscenic grn expr_matrix_RNA_counts.tsv \
    /work/project/xuanyu/resource/SCENIC/TFs/allTFs_hg38.txt \
    -o adj.csv \
    --num_workers 20 \
    --transpose

pyscenic grn expr_matrix_RNA_counts_56339cells.tsv \
    /work/project/xuanyu/resource/SCENIC/TFs/allTFs_hg38.txt \
    -o adj_56339cells.csv \
    --num_workers 20 \
    --transpose


pyscenic ctx \
    --annotations_fname /work/project/xuanyu/resource/SCENIC/motif_to_TF_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname expr_matrix_RNA_counts.tsv \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 20 \
    --transpose \
    adj_56339cells.csv \
    /work/project/xuanyu/resource/SCENIC/cisTarget_databases/human/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    /work/project/xuanyu/resource/SCENIC/cisTarget_databases/human/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

python /work/project/xuanyu/script/pySCENIC/regulon_reformat.py

pyscenic aucell \
    expr_matrix_RNA_counts.tsv \
    reg.csv \
    --output aucell.csv \
    --num_workers 20 \
    --transpose



#-------------mouse---------------------------------------------
pyscenic grn expr_matrix_RNA_counts.tsv \
    /work/project/xuanyu/resource/SCENIC/TFs/allTFs_mm.txt \
    -o adj.csv \
    --num_workers 20 \
    --transpose

pyscenic ctx \
    --annotations_fname /work/project/xuanyu/resource/SCENIC/motif_to_TF_annotation/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname expr_matrix_RNA_counts.tsv \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 20 \
    --transpose \
    adj.csv \
    /work/project/xuanyu/resource/SCENIC/cisTarget_databases/mouse/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    /work/project/xuanyu/resource/SCENIC/cisTarget_databases/mouse/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

python /work/project/xuanyu/script/pySCENIC/regulon_reformat.py

pyscenic aucell \
    expr_matrix_RNA_counts.tsv \
    reg.csv \
    --output aucell.csv \
    --num_workers 20 \
    --transpose


#-------------regulonActivity_combinded.R
library(Seurat)
library(dplyr)
regulonAUC <- read.table(file='aucell.csv',sep=',',header=T,row.names=1,check.names=F)
load('/work/project/xuanyu/project/ISP/seurat/final/CM/seuObj.integrated_CM.Rdata')
#allcells <- colnames(seuObj.integrated)
#seuObj.integrated <- RenameCells(seuObj.integrated, new.names = sub("^_", "", allcells))
# create a v3 assay
seuObj.integrated[["regulonAUC"]] <- CreateAssayObject(counts = regulonAUC)

DefaultAssay(seuObj.integrated) <- 'regulonAUC'
Idents(seuObj.integrated) <- 'group'

ident_1 <- "CASE"
ident_2 <- "CTRL"

DEGs <- FindMarkers(object =seuObj.integrated, ident.1 = ident_1, ident.2 = ident_2, only.pos = FALSE, min.pct = 0.9,
 logfc.threshold = 0.8, test.use="wilcox")
DEGs$regulon <- row.names(DEGs)
DEGs$updown <- ifelse(DEGs$avg_log2FC>0,"up","down")                 
DEGs_sorted <- DEGs %>% arrange(desc(updown))
dim(DEGs_sorted)
write.table(DEGs_sorted,file=paste0(ident_1,"_vs_",ident_2,"_DE_regulons.tsv"),sep="\t",row.names=F,quote=F,col.names=T)

# Identify differential expressed genes across condition
Idents(seuObj.integrated) <- 'subcluster'
markers.all <- FindAllMarkers(object =seuObj.integrated, only.pos = TRUE, min.pct = 0.1,
 logfc.threshold =0.25, test.use="wilcox")
write.table(markers.all,file="regulons.subcluster.tsv",sep="\t",row.names=F,quote=F,col.names=T)

head(subset(markers.all,cluster == 'Ad4'),n=20)

save(seuObj.integrated, file="seuObj.integrated.regulonAUC.RData")
#--------------- visualization---------------

library(Seurat)
library(ggplot2)
library(viridis)
library(Nebulosa)
library(ggExtra)
suppressPackageStartupMessages(library(escape))
library(ggthemes)
suppressPackageStartupMessages(library(scRNAtoolVis))
library(reshape2)
library(plyr)

# change the current plan to access parallelization
plan("multisession", workers = 10)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)


markers.all_grouped <- markers.all %>%
  group_by(cluster) %>%
  slice(1:5) %>%
  ungroup()

# dot plot
#features=markers.all_grouped$gene
features=rev(markers.all$gene)
jjDotPlot(object = seuObj.integrated,
          gene = features,
          dot.col = c('blue','white','red'),
          id='subcluster',
          ytree=F,
          legend.position='top',
          x.text.angle = 45,x.text.hjust=0.5,x.text.vjust=0.5)
ggsave(file='top5Regulons_subcluster_jjDotplot.pdf')


# plot regulon activity score
features <- DEGs_sorted$regulon
cluster.order=c('CM1 (CTRL)','CM1 (CASE)','CM2 (CTRL)','CM2 (CASE)')
jjDotPlot(object = seuObj.integrated,
          gene = rev(features),
          dot.col = c('blue3','white','red3'),
          id = 'subcluster',
          split.by = 'group',
          split.by.aesGroup = T,ytree=F,legend.position='right',
          slot="data",
          assay="regulonAUC",
          cluster.order=cluster.order,
          x.text.angle = 45,x.text.hjust=0.5,x.text.vjust=0.5) + coord_flip()
ggsave(file="DR_split_jjDotplot.pdf")

# plot TF expression
features <- gsub("\\(\\+\\)$", "", DEGs_sorted$regulon)
jjDotPlot(object = seuObj.integrated,
          gene = rev(features),
          dot.col = c('blue3','white','red3'),
          id = 'subcluster',
          split.by = 'group',
          split.by.aesGroup = T,ytree=F,legend.position='right',
          slot="data",
          assay="RNA",
          dot.min=3,
          dot.max=7,
          cluster.order=cluster.order,
          x.text.angle = 45,x.text.hjust=0.5,x.text.vjust=0.5) + coord_flip()

ggsave(file="TF_split_jjDotplot.pdf")





# grouped heatmap

DefaultAssay(seuObj.integrated) <- 'regulonAUC'
Idents(seuObj.integrated) <- 'group'
vars_to_regress <- c("nCount_RNA","nFeature_RNA","percent.mito","S.Score","G2M.Score")
seuObj.integrated <- ScaleData(object = seuObj.integrated,, features=rownames(seuObj.integrated),vars.to.regress = vars_to_regress,
  verbose = T,model.use = "linear", do.scale = TRUE,do.center = TRUE)
features <- DEGs_sorted$regulon
SCpubr::do_ExpressionHeatmap(sample = seuObj.integrated,
                                  assay = 'regulonAUC',
                                  features = features,
                                  group.by = c("group","subcluster"),
                                  cluster=FALSE,
                                  features.order=rev(features),
                                  enforce_symmetry = TRUE,
                                  flip=TRUE,use_viridis =F,
                                  max.cutoff=0.25,
                                  min.cutoff=-0.25,
                                  viridis.palette = "A",viridis.direction = 1,
                                  slot = "scale.data")

ggsave(file='grouped_heatmap_Regulon.pdf')

#---------------------- genes -------------------------------------------
#for human 
#geneTypeInfo <- read.table(file='/work/project/xuanyu/resource/10XGenomics/geneSymbol_geneType.tsv',header=F,sep='\t',stringsAsFactors=F)
DefaultAssay(seuObj.integrated) <- 'RNA'
geneTypeInfo <- read.table(file='/work/project/xuanyu/resource/10XGenomics/geneSymbol_geneType.tsv',header=F,sep='\t',stringsAsFactors=F)
protein_coding_genes <- subset(geneTypeInfo,V2=="protein_coding")$V1
allGenes <- rownames(seuObj.integrated)
allProteinCodingGenes <- allGenes[allGenes %in% protein_coding_genes]

vars_to_regress <- c("nCount_RNA","nFeature_RNA","percent.mito","S.Score","G2M.Score")
seuObj.integrated <- ScaleData(object = seuObj.integrated,assay='RNA', features=allProteinCodingGenes,vars.to.regress = vars_to_regress,
  verbose = T,model.use = "linear", do.scale = TRUE,do.center = TRUE)

genes <- gsub("\\(\\+\\)$", "", rownames(DEGs_sorted))

# grouped heatmap
SCpubr::do_ExpressionHeatmap(sample = seuObj.integrated,
                                  assay = 'RNA',
                                  features = genes,
                                  group.by = c("group","subcluster"),
                                  cluster=FALSE,
                                  features.order=rev(genes),
                                  enforce_symmetry = TRUE,
                                  flip=TRUE,use_viridis =F,
                                  viridis.palette = "A",viridis.direction = 1,
                                  slot = "scale.data")
ggsave(file='grouped_heatmap_TFs.pdf')




#-------Prepare input file for cytoscape ------------------
#-----single regulon -------
regulon_info <- read.table(file='regulon.tsv',header=F,sep='\t')
iRegulon <- "BHLHE40(+)"
TF <- gsub("\\(\\+\\)$", "", iRegulon)
regulon <- subset(regulon_info,V1%in%iRegulon)$V2
regulonVec<- strsplit(regulon,";")[[1]]
outdf <- data.frame(source=rep(TF,length(regulonVec)),target=regulonVec)
write.table(outdf,file=paste0(TF,'_regulon_cytoscape_input.tsv'),quote=F,sep='\t',col.names=T,row.names=F)

#---------multiple regulons----------------
iRegulons <- c("Ar(+)","Ebf3(+)","E2f1(+)")
TFs <- gsub("\\(\\+\\)$", "", iRegulons)
regulons <- subset(regulon_info,V1%in%iRegulons)$V2
regulonVec<- unlist(strsplit(regulons,";"))
vecLenth<- as.numeric(lapply(strsplit(regulons,";"),length))
# Repeat each element of 'a' according to the values in 'b'
repeated_a <- as.character(unlist(sapply(seq_along(TFs), function(i) rep(TFs[i], vecLenth[i]))))
out2df <- data.frame(source=repeated_a,target=regulonVec)                           
data.frame(source=rep(TFs,length(regulonVec)),target=regulonVec)
write.table(out2df,file=paste0("TFs",'_regulon_cytoscape_input.tsv'),quote=F,sep='\t',col.names=T,row.names=F)


