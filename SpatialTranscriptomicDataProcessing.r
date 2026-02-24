library(Seurat)
library(future)
library(plyr)
plan("multicore", workers = 10)
options(future.globals.maxSize = 200 * 1024 ^ 3)
plan()



data_dir <- '/work/project/xuanyu/project/MVD/ST/spaceranger/s87393/outs/'
STdata.s87393 <- Load10X_Spatial(data.dir = data_dir)
STdata.s87393@meta.data$sampleID <-rep('s87393',ncol(STdata.s87393))

data_dir <- '/work/project/xuanyu/project/MVD/ST/spaceranger/s72791/outs/'
STdata.s72791 <- Load10X_Spatial(data.dir = data_dir)
STdata.s72791@meta.data$sampleID <-rep('s72791',ncol(STdata.s72791))


data_dir <- '/work/project/xuanyu/project/MVD/ST/spaceranger/s88613/outs/'
STdata.s88613 <- Load10X_Spatial(data.dir = data_dir)
STdata.s88613@meta.data$sampleID <-rep('s88613',ncol(STdata.s88613))


data_dir <- '/work/project/xuanyu/project/MVD/ST/spaceranger/s93863/outs/'
STdata.s93863 <- Load10X_Spatial(data.dir = data_dir)
STdata.s93863@meta.data$sampleID <-rep('s93863',ncol(STdata.s93863))


seuObj.list <- c(s87393=STdata.s87393,s88613=STdata.s88613,s72791=STdata.s72791,s93863=STdata.s93863)

# ---------- Gene expression data processing ------------------
# SCTransform is designed to eliminate the influence of sequencing depth. The nCount_RNA is used to construct the sct mode, so you don't need to put it into vars.to.regress.
seuObj.list <- lapply(X = seuObj.list, function(i){SCTransform(i, method = "glmGamPoi",assay='Spatial',vst.flavor = "v2")})
# Typically use 3,000 or more features for analysis downstream of sctransform
features <- SelectIntegrationFeatures(object.list = seuObj.list, nfeatures = 3000, assay =c('SCT','SCT','SCT','SCT'))
# Run the PrepSCTIntegration() function prior to identifying anchors
seuObj.list <- PrepSCTIntegration(object.list = seuObj.list, anchor.features = features, assay =c('SCT','SCT','SCT','SCT'))

# integration
GEX.anchors <- FindIntegrationAnchors(object.list = seuObj.list, normalization.method = "SCT", anchor.features = features)
STobj.GEXintegrated <- IntegrateData(anchorset = GEX.anchors, normalization.method = "SCT",new.assay.name = "GEXintegrated",dims = 1:20)

#Perform dimensionality reduction by PCA and UMAP embedding
STobj.GEXintegrated <- RunPCA(STobj.GEXintegrated, verbose = FALSE, assay = "GEXintegrated",npcs = 50)
STobj.GEXintegrated <- RunUMAP(STobj.GEXintegrated, reduction = "pca", assay = "GEXintegrated", dims = 1:20)
STobj.GEXintegrated <- FindNeighbors(STobj.GEXintegrated, reduction = "pca", dims = 1:20)


STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.4)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.6)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 0.8)
STobj.GEXintegrated <- FindClusters(STobj.GEXintegrated, resolution = 1)


# Identify markers
STobj.GEXintegrated <- PrepSCTFindMarkers(object =STobj.GEXintegrated, assay = "SCT", verbose = TRUE)
markers.all <- FindAllMarkers(object =STobj.GEXintegrated, assay = "SCT", only.pos = TRUE, min.pct = 0.25,
 logfc.threshold = 0.25, test.use="wilcox")
write.table(markers.all,file="markers.all.tsv",sep="\t",row.names=F,quote=F,col.names=T)
library(dplyr)
represent_markers <- markers.all %>%  group_by(cluster) %>%  slice(1:5) %>% select(gene)
head(subset(markers.all,cluster == 'sc0'),n=10)

save(STobj.GEXintegrated,file='STobj.GEXintegrated.ALL.RData')



STobj.GEXintegrated@meta.data$seurat_cluster <- paste0('sc',STobj.GEXintegrated@meta.data$GEXintegrated_snn_res.0.6)
STobj.GEXintegrated@meta.data$seurat_cluster <- factor(STobj.GEXintegrated@meta.data$seurat_cluster,levels=c('sc0','sc1','sc2','sc3','sc4','sc5','sc6','sc7','sc8','sc9','sc10','sc11'))
Idents(STobj.GEXintegrated) <- 'seurat_cluster'

# set colors for SpatialDimPlot
# color pallet
library(ggplot2)
categary.color.pallet <- c("#5B74F4","#E82C61","#ff7f0e","#2ca02c","#7733B7",
                           "#8c565B","#E236AF","#8C8C00","#700B35","#06BACE",
                           "#193C3E","#992756","#DADD00","#F76AA7","#C4796A",
                           "#246A73","#FF6360","#4DDD30","#FFA449","#74A6E8",
                           "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
                           "#63C73D","#800026","#4567EA","#234556","#1f7234")

Idents(STobj.GEXintegrated) = STobj.GEXintegrated$seurat_cluster 
colors = categary.color.pallet
library(dplyr)
names(colors) <- Idents(STobj.GEXintegrated) %>% levels()

#STobj.GEXintegrated$group <- mapvalues(STobj.GEXintegrated$sampleID,from=c('','AD2287','YZ72','YZ74'),to=c('diseased','diseased','healthy','healthy'))
# visualization
# Dimplot
library(ggplot2)

# UMAP categary data plot
categary='GEXintegrated_snn_res.0.6'
categary='seurat_cluster'

DimPlot(STobj.GEXintegrated, label = T,raster=FALSE, repel = TRUE, pt.size=1.5,reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'_UMAP_DimPlot_labelT.jpg'),dpi = 800,width=5.88,height=4.76,unit='in')

DimPlot(STobj.GEXintegrated, label = F,raster=FALSE, repel = TRUE, pt.size=1.5, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'_UMAP_DimPlot_labelF.jpg'),dpi = 800,width=5.88,height=4.76,unit='in')


# SpatialDimPlot for all slices
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, image.alpha =0.2, alpha = c(1, 1), pt.size.factor =1.4,
stroke = 0.2,cols=colors,crop=F)
ggsave(file="ALL_slice_Spatial_cluster_DimPlot_with_image.pdf")
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, image.alpha =0, alpha = c(1, 1), pt.size.factor =1.4,
stroke = 0.2,cols=colors,crop=F)
ggsave(file="ALL_slice_Spatial_cluster_DimPlot_without_image.jpg",dpi=800)

# SpatialDimPlot for one slice
s87393 s88613 s72791

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1', image.alpha =0.5, alpha = c(1, 1),pt.size.factor =5,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s87393_Spatial_cluster_DimPlot.pdf")
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1', image.alpha =0, alpha = c(1, 1),pt.size.factor =1.3,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s87393_Spatial_cluster_DimPlot.jpg",dpi = 800)

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.2', image.alpha =0.5, alpha = c(1, 1),pt.size.factor =1.3,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s88613_Spatial_cluster_DimPlot.pdf")
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.2', image.alpha =0, alpha = c(1, 1),pt.size.factor =1.3,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s88613_Spatial_cluster_DimPlot.jpg",dpi = 800)

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.3', image.alpha =0.5, alpha = c(1, 1),pt.size.factor =1.3,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s72791_Spatial_cluster_DimPlot.pdf")
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.3', image.alpha =0, alpha = c(1, 1),pt.size.factor =1.3,
stroke = 0.5,cols=colors,crop=F)
ggsave(file="s72791_Spatial_cluster_DimPlot.jpg",dpi = 800)




# highlight a specific cluster
SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1', image.alpha =0, alpha = c(1, 1),pt.size.factor =2,
stroke = 1, cells.highlight = CellsByIdentities(object = STobj.GEXintegrated, idents = c('sc0', 'sc1')), facet.highlight = TRUE, ncol = 4)

clustersTohighlight <- levels(Idents(STobj.GEXintegrated))

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1', image.alpha =0, alpha = c(1, 1),pt.size.factor =2,
stroke = 1, cells.highlight = CellsByIdentities(object = STobj.GEXintegrated, idents =clustersTohighlight), facet.highlight = TRUE, ncol = 4)
ggsave(file="slice1_STcluster_highlight.jpg",dpi=600,width=17.1,height=11.9,unit='in')

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.2', image.alpha =0, alpha = c(1, 1),pt.size.factor =1.5,
stroke = 1,cells.highlight = CellsByIdentities(object = STobj.GEXintegrated, idents =clustersTohighlight), facet.highlight = TRUE, ncol = 4)
ggsave(file="slice1.2_STcluster_highlight.jpg",dpi=600,width=17.1,height=11.9,unit='in')

SpatialDimPlot(STobj.GEXintegrated, label = FALSE, images='slice1.3', image.alpha =0, alpha = c(1, 1),pt.size.factor =1.5,
stroke = 1, cells.highlight = CellsByIdentities(object = STobj.GEXintegrated, idents =clustersTohighlight), facet.highlight = TRUE, ncol = 4)
ggsave(file="slice1.3_STcluster_highlight.jpg",dpi=600,width=17.1,height=11.9,unit='in')


#Gene expression visualization
DefaultAssay(STobj.GEXintegrated) <- 'SCT'
feature='CHAD'
feature='FMOD'
feature='SFRP1'
feature='CNMD'
feature='IGFBP4'
feature='CRLF1'
feature='CTHRC1'
feature='MYH11'
feature='nCount_Spatial'
feature='nCount_SCT'
feature='nFeature_Spatial'
feature='nFeature_SCT'
feature='RGS5'

# all slices
features =c('FMOD','SFRP1','CNMD','PMEPA1','MYH11','RGS5','NR4A2','RARRES2','PCOLCE2','ISG15','MFAP5','NTM','COL1A1','nCount_SCT','nFeature_SCT')
for (thefeature in features) {
print(thefeature)
SpatialFeaturePlot(STobj.GEXintegrated, image.alpha= 0, features = thefeature, crop= F, pt.size.factor = 1.4,alpha = c(1, 1),stroke=0.25)
ggsave(file=paste0(thefeature,'_spatialFearturePlot_rainbow.jpg'),dpi = 800,width=14.7,height=8.8,units='in')
}

# one slice
SpatialFeaturePlot(STobj.GEXintegrated, image.alpha= 0, images='slice1.2', features = feature, crop= F, pt.size.factor = 1.4,alpha = c(1, 1),stroke=0.25)

# spatial feature plot
library(viridis)
feature='FMOD'
feature='SFRP1'
feature='IGFBP4'
feature='BGN'
feature='ACTA2'
feature='ISG15'
feature='RARRES2'
feature='PCOLCE2'
feature='CTHRC1'
feature='PRELP'
feature='nCount_Spatial'
feature='nFeature_Spatial'
feature='PPP1R1B'
feature='GFAP'
feature='LYZ'
feature='MYH11'
feature='ACTA2'
DefaultAssay(STobj.GEXintegrated) <- 'SCSE.REACTOME'
feature='REACTOME-SIGNALING-BY-TGF-BETA-RECEPTOR-COMPLEX'
#magma color palette
for (imageID in c('slice1','slice1.2','slice1.3','slice1.4')) { 
SpatialFeaturePlot(STobj.GEXintegrated, image.alpha= 0,features = feature,images=imageID,crop=F,
 pt.size.factor = 1.3,alpha = c(1, 1)) +  scale_fill_viridis(option="magma")
ggsave(file=paste0(feature,'_',imageID,'_spatialFearturePlot.alpha0.1.magma.jpg'),dpi = 800)
}

#default color palette
for (imageID in c('slice1','slice1.2','slice1.3','slice1.4')) { 
SpatialFeaturePlot(STobj.GEXintegrated, image.alpha= 0,features = feature,images=imageID,crop=F,
 pt.size.factor = 1.3,alpha = c(1, 1)) 
ggsave(file=paste0(feature,'_',imageID,'_spatialFearturePlot.alpha0.1.defaultColor.jpg'),dpi = 800)
}



# Integration with single-cell data
#load('/work/project/xuanyu/project/MVD/ST/seurat/STobj.GEXintegrated.Rdata')
load('/work/project/xuanyu/project/MVD/seurat/round2/ALL/seuObj.integrated_MVD_ALL_final.RData')
#seuObj.integrated <- SCTransform(seuObj.integrated, vars.to.regress = c("nCount_RNA","percent.mito","S.Score", "G2M.Score"))
save(seuObj.integrated,file='/work/project/xuanyu/project/MVD/seurat/round2/ALL/seuObj.integrated_MVD_ALL_with_SCT_final.RData')

DefaultAssay(seuObj.integrated) <- 'SCT'
anchors <- FindTransferAnchors(reference = seuObj.integrated, query =STobj.GEXintegrated, normalization.method = "SCT",
reference.assay ="SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seuObj.integrated$subcluster2, prediction.assay = TRUE,
    weight.reduction = STobj.GEXintegrated[["pca"]], dims = 1:30)
STobj.GEXintegrated[["predictions"]] <- predictions.assay

#Now we get prediction scores for each spot for each subpopulation.
DefaultAssay(STobj.GEXintegrated) <- "predictions"

features <- row.names(STobj.GEXintegrated[["predictions"]]@data)
for (theFeature in features){
SpatialPlot(object = STobj.GEXintegrated, features = theFeature,crop=F, ncol = 4,alpha = c(1, 1),image.alpha = 0,pt.size.factor = 1.3,stroke = 0.25)
ggsave(file=paste0(theFeature,'.spatialPlot.jpg'))
}


# dim density_plot_for_sparse_marker
DefaultAssay(STobj.GEXintegrated) <- 'SCT'

library(Nebulosa)
feature='FMOD'
plot_density(STobj.GEXintegrated, feature, pal="magma",size = 1) +theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.marker_density_plot.jpg'),width=4.89,height=4.17,unit='in',dpi = 800)


# module score
features <- list(c('DCN','BGN','LUM','ACAN'))
STobj.GEXintegrated <- AddModuleScore(STobj.GEXintegrated,feature=features,ctrl = 3,name='PG_features')
SpatialPlot(object = STobj.GEXintegrated, features = 'PG_features1',crop=F, ncol = 4,alpha = c(1, 1),image.alpha = 0,pt.size.factor = 1.3,stroke = 0.25)

# prep SCSE input data
normData <- t(as.matrix(STobj.GEXintegrated@assays$SCT@data))
write.table( normData, "ALL_SCTexpr.tsv", sep="\t", row.names = TRUE, col.names=NA)
mv /work/project/xuanyu/project/MVD/ST/seurat/ALL_SCTexpr.tsv /work/project/xuanyu/software/SCSE/SingleCellSignatureExplorer/SingleCellSignatureScorer/data
#-------------------------Add pathway enrichment score assay -------------------------------------------------------

REACTOME_count <- read.table(file="REACTOME_v2023.2_HUMAN_ALL_SCTexpr.tsv",sep="\t",header=T,row.names=1)
REACTOME_count <-t(REACTOME_count)
# create a v3 assay
STobj.GEXintegrated[["SCSE.REACTOME"]] <- CreateAssayObject(counts = REACTOME_count)


KEGG_count <- read.table(file="KEGGlegacy_v2023.2_HUMAN_ALL_SCTexpr.tsv",sep="\t",header=T,row.names=1)
KEGG_count <-t(KEGG_count)
# create a v3 assay
STobj.GEXintegrated[["SCSE.KEGG"]] <- CreateAssayObject(counts = KEGG_count)


GOBP_count <- read.table(file="GOBP_v2023.2_HUMAN_ALL_SCTexpr.tsv",sep="\t",header=T,row.names=1)
GOBP_count <-t(GOBP_count)
# create a v3 assay
STobj.GEXintegrated[["SCSE.GOBP"]] <- CreateAssayObject(counts = GOBP_count)


#save 
save(STobj.GEXintegrated,file='STobj.GEXintegrated_new.Rdata')
