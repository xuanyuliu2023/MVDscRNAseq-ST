library(Seurat)
library(future)
library(plyr)

# change the current plan to access parallelization
plan("multisession", workers = 20)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

seuObj.list.DMVD <- get(load('seuObj.list.DMVD.scrublet.Rdata'))
seuObj.list.BD <- get(load('seuObj.list.scrublet.BD.Rdata'))
seuObj.list.RMV <- get(load('seuObj.list.scrublet.RMV.Rdata'))
seuObj.list <- c(seuObj.list.DMVD,seuObj.list.BD,seuObj.list.RMV)

#sample names & aggr data path
sampleNames <- c("K707CTRL","K709CTRL","K710CTRL","K711day3","K718day3","K719day3", "K694day7","K695day7","K698day7")


# filter cells
#seuObj.list <- lapply(seuObj.list,function(x){AddMetaData(x, row.names(x@meta.data), col.name = 'cellID')})
#seuObj.list <- lapply(seuObj.list,function(x){subset(x, cellID %in% cellsToKeep)})


# integration
## ----------- harmony integration for small datasets
library(harmony)
library(scCustomize)
seuObj.integrated <- Merge_Seurat_List(
  seuObj.list,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject"
)

# log-normalization
seuObj.integrated <- NormalizeData(seuObj.integrated,normalization.method = "LogNormalize",
    scale.factor = 10000)
seuObj.integrated <- FindVariableFeatures(seuObj.integrated,selection.method = "vst", nfeatures = 2000, verbose = FALSE)


if (species == "human") {    
    seuObj.integrated <- CellCycleScoring(seuObj.integrated, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes,set.ident = F)
    } else {
    load('/work/project/xuanyu/resource/human_mouse_genesymbol_convert/cc.genes.mm.Rdata')
    seuObj.integrated <- CellCycleScoring(seuObj.integrated, s.features = cc.genes.mm$s.genes, g2m.features = cc.genes.mm$g2m.genes,set.ident = F)
}


vars_to_regress <- c("nCount_RNA","nFeature_RNA","percent.mito","S.Score","G2M.Score")
seuObj.integrated <- ScaleData(object = seuObj.integrated,vars.to.regress = vars_to_regress,
  verbose = T,model.use = "linear", do.scale = TRUE,do.center = TRUE)

# or sctransform normalization
seuObj.integrated <- SCTransform(seuObj.integrated, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mito","S.Score","G2M.Score"), verbose = FALSE)

# PCA Harmony
seuObj.integrated <- RunPCA(seuObj.integrated, verbose = FALSE,npcs = 50)
ElbowPlot(seuObj.integrated, ndims = 50, reduction = "pca")
seuObj.integrated <- RunHarmony(seuObj.integrated, group.by.vars="orig.ident")
seuObj.integrated <- RunUMAP(seuObj.integrated , reduction = "harmony", dims = 1:30)

# clustering
seuObj.integrated <- FindNeighbors(seuObj.integrated, reduction = "harmony", dims = 1:30)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.2)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.4)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.6)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.8)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 1.0)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 1.2)



# add meta info

library(plyr)

sampleInfo.df <- read.table(file="/work/project/xuanyu/project/MVD/seurat/sampleInfo.tsv",sep='\t',header=T,stringsAsFactors=F)

seuObj.integrated@meta.data$group <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$group)
seuObj.integrated@meta.data$group <- factor(seuObj.integrated@meta.data$group,levels=c('CTRL','FED','BD','RMVD'))
seuObj.integrated@meta.data$group2 <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$group2)
seuObj.integrated@meta.data$sex <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$sex)
seuObj.integrated@meta.data$age <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$age)
seuObj.integrated@meta.data$BMI <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$BMI)
seuObj.integrated@meta.data$newID <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleInfo.df$oldID,to=sampleInfo.df$newID)
seuObj.integrated@meta.data$newID <- factor(seuObj.integrated@meta.data$newID,levels=c("CTRL1","CTRL3","CTRL4","CTRL5","CTRL6","CTRL7","CTRL9","CTRL10","FED2","FED3","FED4",
"FED6","FED7","FED8","FED9","FED10","FED14","FED15","FED16","FED17","FED22","FED23","FED24","BD1","BD2","BD3","BD4","BD5","BD6","BD7",
"RMVD1","RMVD2","RMVD3","RMVD6","RMVD7","RMVD8","RMVD10","RMVD12","RMVD14","RMVD15","RMVD20"))

# set Idents
seuObj.integrated$seurat_clusters <- paste0('c',seuObj.integrated@meta.data$RNA_snn_res.0.6)
Idents(seuObj.integrated) <- 'seurat_clusters'


# Identify differential expressed genes across condition
markers.all <- FindAllMarkers(object =seuObj.integrated, assay = "RNA", only.pos = TRUE, min.pct = 0.25,
 logfc.threshold = 0.25, test.use="bimod")
write.table(markers.all,file="markers.all.tsv",sep="\t",row.names=F,quote=F,col.names=T)
head(subset(markers.all,cluster == 'c12'),n=20)


#Identify conserved cell type markers across conditions
Idents(seuObj.integrated) <- 'lineage'
adipocyte.markers <- FindConservedMarkers(seuObj.integrated, assay = "SCT", ident.1 = "Adipocyte", grouping.var = "group",verbose = FALSE)
head(adipocyte.markers)

# Phylogenetic Analysis of Identity Classes
library(ape)
seuObj.integrated <- BuildClusterTree(object = seuObj.integrated,assay='RNA',
features=seuObj.integrated@assays$RNA@var.features)
plot(Tool(object = seuObj.integrated, slot = 'BuildClusterTree'))


# add meta info
sample_info = read.table(file="/work/project/xuanyu/project/LDS/sampleInfo.tsv",header=T,sep='\t',stringsAsFactors=F)
library(plyr)
sampleNames <- c("YZ17","YZ19","YZ20","YZ23","YZ38","AD1400","AD1449","AD1501","AD2437","AD793","AD852","AD1061","AD1346","AD2017","AD2215",
"AD2271","AD997","AD1037","AD1201","AD1569","AD1669","AD1915","AD2654")

seuObj.integrated@meta.data$group <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sample_info$sampleID,to=sample_info$group)
seuObj.integrated@meta.data$sex <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sample_info$sampleID,to=sample_info$sex)
seuObj.integrated@meta.data$mutatedGene <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sample_info$sampleID,to=sample_info$mutatedGene)
seuObj.integrated@meta.data$age <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sample_info$sampleID,to=sample_info$age)


# save
save(seuObj.integrated,file='seuObj.integrated.ALL.RData')


#discard
cells2discard <- row.names(subset(meta.data.table,RNA_snn_res.0.8 %in% c('20','22')))
save(cells2discard,file='cells2discard.Rdata')


#----------------------subset----------------------------
seuObj.integrated.CM <- subset(seuObj.integrated, cellType=="CM")
seuObj.integrated.CM@meta.data$subcluster <- factor(seuObj.integrated.CM@meta.data$subcluster,levels=c('CM1','CM2'))

seuObj.integrated.CM <- FindVariableFeatures(seuObj.integrated.CM,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
vars_to_regress <- c("nCount_RNA","nFeature_RNA","percent.mito","S.Score","G2M.Score")
seuObj.integrated.CM <- ScaleData(object = seuObj.integrated.CM,vars.to.regress = vars_to_regress,
  verbose = T,model.use = "linear", do.scale = TRUE,do.center = TRUE)
# PCA Harmony
seuObj.integrated.CM <- RunPCA(seuObj.integrated.CM, verbose = FALSE,npcs = 50)
seuObj.integrated.CM <- RunHarmony(seuObj.integrated.CM, group.by.vars="orig.ident")
seuObj.integrated.CM <- RunUMAP(seuObj.integrated.CM, reduction = "harmony", dims = 1:20)

# clustering
seuObj.integrated.CM <- FindNeighbors(seuObj.integrated.CM, reduction = "harmony", dims = 1:20)
seuObj.integrated.CM <- FindClusters(seuObj.integrated.CM, resolution = 0.2)

save(seuObj.integrated.CM, file='seuObj.integrated.CM.Rdata')