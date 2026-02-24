library(Seurat)
library(future)
library(plyr)
plan("multicore", workers = 40)
options(future.globals.maxSize = 200 * 1024 ^ 3)
plan()

species <- 'human' #or mouse


cellrangerDir <- "/work/project/xuanyu/project/LDS/cellranger/"
sampleNames <- c("YZ17","YZ19","YZ20","YZ23","YZ38","AD1400","AD1449","AD1501","AD2437","AD793","AD852","AD1061","AD1346","AD2017","AD2215",
"AD2271","AD997","AD1037","AD1201","AD1569","AD1669","AD1915","AD2654")

seuObj.list <- list()

#import data & create seurate object
for(i in 1:length(sampleNames)) {
    print(i)
    dataDir<- paste0(cellrangerDir,sampleNames[i],"/outs/filtered_feature_bc_matrix")
    print(dataDir)
    counts <- Read10X(data.dir = dataDir)
    seuObj.list[[i]] <- CreateSeuratObject(counts, project = sampleNames[i], assay = "RNA",min.cells = 3, min.features = 0,
    names.field = 2,names.delim = "-", meta.data = NULL)
    #seuObj.list[[i]] <- RenameCells(seuObj.list[[i]], add.cell.id = sampleNames[i])
    seuObj.list[[i]] <- RenameCells(seuObj.list[[i]], new.names = paste0(colnames(seuObj.list[[i]]),"___",sampleNames[i]))
    seuObj.list[[i]]@meta.data$orig.ident <- mapvalues(seuObj.list[[i]]@meta.data$orig.ident,from="1",to=sampleNames[i])
}
names(seuObj.list) <- sampleNames

# calculate percent.mito (%)
for (i in 1:length(x = seuObj.list)) {
        if (species == "human") {
        mito.features <- grep(pattern = "^MT-", x = rownames(x = seuObj.list[[i]]), value = TRUE)
        } else {
        mito.features <- grep(pattern = "^mt-", x = rownames(x = seuObj.list[[i]]), value = TRUE)
        }
        seuObj.list[[i]] <- PercentageFeatureSet(seuObj.list[[i]],features=mito.features,col.name="percent.mito")
}

# calculate percent.ribo (%)
for (i in 1:length(x = seuObj.list)) {
        if (species == "human") {
        riboProtGenes <- read.table(file='/work/project/xuanyu/resource/ribosomalProteinGenes/ribosomalProteinGenes.human.tsv',header=F,stringsAsFactors=F)$V1
        } else {
        riboProtGenes <- read.table(file='/work/project/xuanyu/resource/ribosomalProteinGenes/ribosomalProteinGenes.nouse.tsv',header=F,stringsAsFactors=F)$V1
        }
        ribo.features <- riboProtGenes[riboProtGenes %in% rownames(x = seuObj.list[[i]])]
        seuObj.list[[i]] <- PercentageFeatureSet(seuObj.list[[i]],features=ribo.features,col.name="percent.ribo")
}

# red blood cell score
for (i in 1:length(x = seuObj.list)) {
    if (species == "human") {
    haemoglobinGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2","HBM", "HBQ1","HBZ")
    } else {
    haemoglobinGenes <- c("Hba-x","Hba-a1","Hba-a2","Hbb-bt","Hbb-bs", "Hbb-bh2", "Hbb-bh1", "Hbb-y","Hbq1a","Hbq1b")
    }
    haemoglobinGenes <- CaseMatch(search =haemoglobinGenes, match = rownames(x = seuObj.list[[i]]))
    seuObj.list[[i]] <- PercentageFeatureSet(seuObj.list[[i]],features=haemoglobinGenes,col.name="percent.rbc")
}




save(seuObj.list,file='seuObj.list.23samps.raw.Rdata')
#-------------------------------------------------------------QC--------------------------------------
seuObj.list.filtered <- seuObj.list
# plot density plot
library(MASS)
library(viridis)
library(scales)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(base_size = 16))

########## get_density funtion calculate the density, the smaller the n, the bigger the grid
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# IQR based cell filtering (half the IQR equals the median absolute deviation (MAD).)
outlierFinder <- function(data){
    lowerq = quantile(data)[2]
    upperq = quantile(data)[4]
    iqr = upperq - lowerq
    mild.threshold.upper = (iqr * 1.5) + upperq
    mild.threshold.lower = lowerq - (iqr * 1.5)
    print(as.numeric(mild.threshold.upper))
    print(as.numeric(mild.threshold.lower))
}
#-------------------------------------------------------------------------
sampleNames <- c("CTRL9","RMV2","RMV3","RMV6","RMV7","RMV8","RMV10","RMV14","RMV15")

for (i in sampleNames) {
print(i)
sampleID = i
sampleIndex = which(sampleID == sampleNames)
print(sampleIndex)
# density Plot
nCount_RNA <- ggplot(seuObj.list[[sampleIndex]]@meta.data,aes(x=nCount_RNA)) + geom_density(fill="#2E9FDF",alpha=.5) + theme_bw()
nFeature_RNA <- ggplot(seuObj.list[[sampleIndex]]@meta.data,aes(x=nFeature_RNA)) + geom_density(fill="#2E5400",alpha=.5) + theme_bw()
percent.mito <- ggplot(seuObj.list[[sampleIndex]]@meta.data,aes(x=percent.mito)) + geom_density(fill="#FF6666",alpha=.5) + theme_bw()
percent.ribo <- ggplot(seuObj.list[[sampleIndex]]@meta.data,aes(x=percent.ribo)) + geom_density(fill="#E7B800",alpha=.5) + theme_bw()
figure <- ggarrange(nCount_RNA, nFeature_RNA, percent.mito, percent.ribo,
                    labels = c("a", "b", "c","d"),
                    ncol = 2, nrow = 2)
#add title
annotate_figure(figure, top = text_grob(sampleID,
               color = "red", face = "bold", size = 14))
ggsave(file=paste0(sampleID,"_combined_QC_metrics_densityPlot.pdf"),width=9.92,height=5.82,units="in")
}


# interactive plot
#ggplotly(nCount_RNA)
#ggplotly(nFeature_RNA)
#ggplotly(percent.ribo)

#--------------------------------------------outlier QC for each sample ------------------------------
sampleNames <- c("CTRL9","RMV2","RMV3","RMV6","RMV7","RMV8","RMV10","RMV14","RMV15")

sampleID = 'RMV15'
sampleIndex = which(sampleID == sampleNames)
# 2D plot
meta.data <- seuObj.list[[sampleIndex]]@meta.data
dat <- meta.data[,c('nCount_RNA','percent.mito','nFeature_RNA','percent.ribo')]

# calculate density
dat$density1 <- get_density(dat$nCount_RNA,dat$percent.mito,n=100)
dat$density2 <- get_density(dat$nFeature_RNA,dat$percent.mito,n=100)

# set the cutoff
outlierFinder(seuObj.list[[sampleIndex]]@meta.data$percent.mito)
outlierFinder(seuObj.list[[sampleIndex]]@meta.data$nCount_RNA)
outlierFinder(seuObj.list[[sampleIndex]]@meta.data$nFeature_RNA)
outlierFinder(seuObj.list[[sampleIndex]]@meta.data$percent.ribo)

#-----------------------------
percent.mito_upperCutoff=7.6
percent.mito_lowerCutoff=0
nCount_RNA_upperCutoff=15238
nCount_RNA_lowerCutoff=1500
nFeature_RNA_upperCutoff=4647
nFeature_RNA_lowerCutoff=500
percent.ribo_upperCutoff=26
percent.ribo_lowerCutoff=0

# plot percent.mito_nCountRNA with lines
percent.mito_nCountRNA <- ggplot(dat,aes(x=nCount_RNA,y=percent.mito)) + geom_point(aes(x=nCount_RNA,y=percent.mito, color = density1),size=1)  +
geom_density_2d() +scale_color_viridis() + labs(x="nCount RNA",y='percent.mito') + theme(legend.position='none') +
scale_y_continuous(minor_breaks = seq(0, max(dat$percent.mito),1)) +
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
annotation_logticks(sides = "lb") + ggtitle(sampleID) +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
geom_vline(xintercept=nCount_RNA_lowerCutoff, linetype="dashed", color = "red") +
geom_vline(xintercept=nCount_RNA_upperCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=percent.mito_lowerCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=percent.mito_upperCutoff, linetype="dashed", color = "red")
#ggsave(file=paste0(sampleID,'_percent.mito_nCountRNA.addlines.pdf'),height=5,width=5,units='in')

# plot nFeature_RNA_nCountRNA with lines
nFeature_RNA_nCountRNA <- ggplot(dat,aes(x=nCount_RNA,y=nFeature_RNA)) + geom_point(aes(x=nCount_RNA,y=nFeature_RNA, color = density1),size=1)  +
geom_density_2d() +scale_color_viridis() + labs(x="nCount RNA",y='nFeature_RNA') + theme(legend.position='none') +
scale_y_continuous(minor_breaks = seq(0, max(dat$nFeature_RNA),1)) +
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
annotation_logticks(sides = "lb") + ggtitle(sampleID) +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
geom_vline(xintercept=nCount_RNA_lowerCutoff, linetype="dashed", color = "red") +
geom_vline(xintercept=nCount_RNA_upperCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=nFeature_RNA_lowerCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=nFeature_RNA_upperCutoff, linetype="dashed", color = "red")
#ggsave(file=paste0(sampleID,'_nFeature_RNA_nCountRNA.addlines.pdf'),height=5,width=5,units='in')

# plot percent.ribo_nCountRNA with lines
percent.ribo_nCountRNA <- ggplot(dat,aes(x=nCount_RNA,y=percent.ribo)) + geom_point(aes(x=nCount_RNA,y=percent.ribo, color = density1),size=1)  +
geom_density_2d() +scale_color_viridis() + labs(x="nCount RNA",y='percent.ribo') + theme(legend.position='none') +
scale_y_continuous(minor_breaks = seq(0, max(dat$percent.ribo),1)) +
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
annotation_logticks(sides = "lb") + ggtitle(sampleID) +
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
geom_vline(xintercept=nCount_RNA_lowerCutoff, linetype="dashed", color = "red") +
geom_vline(xintercept=nCount_RNA_upperCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=percent.ribo_lowerCutoff, linetype="dashed", color = "red") +
geom_hline(yintercept=percent.ribo_upperCutoff, linetype="dashed", color = "red")
#ggsave(file=paste0(sampleID,'_percent.ribo_nCountRNA.addlines.pdf'),height=5,width=5,units='in')
figure <- ggarrange(percent.mito_nCountRNA, nFeature_RNA_nCountRNA, percent.ribo_nCountRNA,
                    labels = c("a", "b", "c"),
                    ncol = 2, nrow = 2)
#add title
annotate_figure(figure, top = text_grob(sampleID,
               color = "red", face = "bold", size = 14))
ggsave(file=paste0(sampleID,"_combined_QC_metrics_2DPlot.pdf"),width=12.1,height=7.42,units="in")


# filter out low quality cells
print(paste0(sampleID," cell number before filtering ",ncol(seuObj.list[[sampleIndex]])))
seuObj.list.filtered[[sampleIndex]] <- subset(x =seuObj.list[[sampleIndex]], subset = nFeature_RNA < nFeature_RNA_upperCutoff & nFeature_RNA > nFeature_RNA_lowerCutoff
              & nCount_RNA < nCount_RNA_upperCutoff & nCount_RNA > nCount_RNA_lowerCutoff
              & percent.mito < percent.mito_upperCutoff & percent.mito > percent.mito_lowerCutoff & percent.ribo < percent.ribo_upperCutoff & percent.ribo > percent.ribo_lowerCutoff)
print(paste0(sampleID," cell number after filtering ",ncol(seuObj.list.filtered[[sampleIndex]])))


# save
save(seuObj.list.filtered,file='seuObj.list.outlierfilter.BD.Rdata')

#### seu2cellRanger
dir.create('./seu2cellRanger')
setwd('./seu2cellRanger')
library(Matrix)
# seurat to cellRanger output
seu2cellRanger <- function(sample){
  sampleIndex = which(sample == sampleNames)
  sparse.gbm <- seuObj.list.filtered[[sampleIndex]]@assays$RNA@counts
  writeMM(obj = sparse.gbm, file=paste0(sample,".matrix.mtx"))
  write(x = rownames(sparse.gbm), file=paste0(sample,".genes.tsv"))
  write(x = colnames(sparse.gbm), file=paste0(sample,".barcodes.tsv"))
}
lapply(sampleNames,seu2cellRanger)

##### ------------------------------------------
/work/project/xuanyu/software/python/python3/Python-3.7.1/python
###python script to remove doublets with scrublet
sampleNames = ["CTRL9","RMV2","RMV3","RMV6","RMV7","RMV8","RMV10","RMV14","RMV15"]

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

#import data
for sample in sampleNames:
  counts_matrix = scipy.io.mmread(sample+'.matrix.mtx').T.tocsc()
  genes = np.array(scr.load_genes(sample+'.genes.tsv', delimiter='\t', column=0))
  print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
  print('Number of genes in gene list: {}'.format(len(genes)))
  scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06,sim_doublet_ratio=2)
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                          min_cells=3,
                                                          min_gene_variability_pctl=85,
                                                          n_prin_comps=30)
  thresh=0.35
  scrub.call_doublets(threshold=thresh)
  predicted_doublets=doublet_scores >thresh
  predicted_doublets.tofile(sample+'.predicted_doublets.csv',sep=',',format='%s')
  sum(predicted_doublets)

#------------------------
### rmDouplet
while read id;do echo cat ${id}.predicted_doublets.csv\|sed \'s\/\,\/\\n\/g\' \>${id}.predicted_doublets.tsv;done < id.tsv|bash

# remove predicted doublets from seurat object
library(Seurat)
load('seuObj.list.outlierfilter.Rdata')
sampleNames <- c("CTRL9","RMV2","RMV3","RMV6","RMV7","RMV8","RMV10","RMV14","RMV15")

for (sample in sampleNames) {
sampleIndex = which(sample == sampleNames)
doublets <- read.table(file=paste0("seu2cellRanger","/",sample,".predicted_doublets.tsv"),header=F)$V1
seuObj.list.filtered[[sampleIndex]]$doublet <- doublets
seuObj.list.filtered[[sampleIndex]] <- subset(x = seuObj.list.filtered[[sampleIndex]], doublet == "False")
}
save(seuObj.list.filtered,file="seuObj.list.scrublet.Rdata")