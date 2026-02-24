
conda activate Seurat4.3.0
library(Seurat)
library(DESeq2)
load('/work/project/xuanyu/project/MVD/seurat/round2/seuObj.integrated_MVD_ALL.Rdata')


# 体量的稀疏矩阵转换为稠密矩阵，使用C++函数 as_matrix
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

sampleNames <- c("CTRL1","CTRL3","CTRL4","CTRL5","CTRL6","CTRL7","CTRL9","CTRL10",
                 "FED2","FED3","FED4","FED6","FED7","FED8","FED9","FED10","FED14",
                 "FED15","FED16","FED17","FED22","FED23","FED24","BD1","BD2","BD3",
                 "BD4","BD5","BD6","BD7","RMVD1","RMVD2","RMVD3","RMVD6","RMVD7","RMVD8"
                 ,"RMVD10","RMVD12","RMVD14","RMVD15","RMVD20")



#----------------------- concise code --------------------------------------
# Assuming that 'seuObj.integrated' is already defined and 'RNA' is an appropriate slot in 'seuObj.integrated'
counts <- as.data.frame(as_matrix(seuObj.integrated@assays$RNA@counts))
meta.data.table <- seuObj.integrated@meta.data


library("BiocParallel")
register(MulticoreParam(10))



#----------------------- indicate celltype ------START HERE!!!
celltype = "DC"

> table(meta.data.table$cellType2,meta.data.table$group)

               CTRL    FED     BD   RMVD
  DC            711    298    113    512
  Macrophage   4574   2490   1119   3619
  Mast_cell      22     26      1     56
  NK              1     12      1     84
  Pericyte       33    205     35    356
  SMC           467   1494     89   1680
  TC           1228    602     84    637
  valvularEC    143    356    208    429
  vascularEC     41     37      2    712
  VIC         45409 113323  50772  72190

# Using a vectorized approach to create a data frame with the total expression for each sample
Results <- sapply(
  sampleNames,
  function(sample) {
    # Ensure that the subsetting results in a data frame and use drop = FALSE
    subset_indices <- row.names(subset(meta.data.table, newID == sample & cellType2 == celltype, drop = FALSE))
    # Only attempt rowSums if there are more than 2 indices to sum over
    if (length(subset_indices) >= 2) {
      rowSums(counts[, subset_indices, drop = FALSE])
    } else {
      # Return NA or some other indicator if the condition is not met
      NA
    }
  }
)

if (type(Results)=='double') {
  df <- Results
  } else {
df <- as.data.frame(Results[!is.na(Results)])
  }


# Remove samples that do not meet the condition (where results are NA)

valid_samples <- colnames(df)
print(valid_samples)
print(paste0("There are ",length(valid_samples),"/",length(sampleNames)," valid samples, which are listed above"))

# only protein-coding genes
#geneType <- read.table(file='/work/project/xuanyu/resource/10XGenomics/geneSymbol_geneType.tsv',sep='\t',header=F)
#protein_coding_genes <- subset(geneType, V2 == 'protein_coding')$V1
#df <- df[row.names(df) %in% protein_coding_genes,]

library(plyr)
subgroup <- mapvalues(colnames(df),from= as.character(meta.data.table$newID),to= as.character(meta.data.table$group))
sex <- mapvalues(colnames(df),from= as.character(meta.data.table$newID),to= as.character(meta.data.table$sex))
group2 <- mapvalues(colnames(df),from= as.character(meta.data.table$newID),to= as.character(meta.data.table$group2))

coldata <- data.frame(condition =subgroup, row.names=colnames(df))
coldata$sex <- sex
coldata$group2 <- group2
save(coldata,file=paste0(celltype,'_coldata.Rdata'))



dds_all_groups <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design= ~ sex + condition)

# set the reference level
dds_all_groups$condition <- factor(dds_all_groups$condition, levels = c('CTRL','FED','BD','RMVD'))
dds_all_groups <- DESeq(dds_all_groups)

#plotCounts(dds_all_groups, gene='CTHRC1', intgroup="condition")

# Extracting normalized counts
nomalized_counts_all_groups <- counts(dds_all_groups, normalized=TRUE)
save(nomalized_counts_all_groups,file=paste0(celltype,"_nomalized_counts.Rdata"))

# Extracting normalized and transformed values  
vsd_all_groups <- vst(dds_all_groups, blind=FALSE) #  normalized counts that were transformed using Variance Stabilizing Transformation (VST); ’
rld_all_groups <- rlog(dds_all_groups, blind=FALSE) # regularized log2 transcformed counts to the log2 scale
#head(assay(vsd_all_groups), 3)
save(dds_all_groups,file=paste0(celltype,"_dds_all_groups.Rdata"))
save(vsd_all_groups,file=paste0(celltype,"_vsd_all_groups.Rdata"))
save(rld_all_groups,file=paste0(celltype,"_rld_all_groups.Rdata"))


#subset the data; indicate group info ---------------------------
CTRL_group = 'CTRL'
CASE_group = 'FED'

CTRL_group = 'CTRL'
CASE_group = 'BD'

CTRL_group = 'CTRL'
CASE_group = 'RMVD'

CTRL_group = 'FED'
CASE_group = 'BD'

CTRL_group = 'FED'
CASE_group = 'RMVD'

CTRL_group = 'BD'
CASE_group = 'RMVD'

coldata_subset = subset(coldata,condition %in% c(CTRL_group,CASE_group))
print("The number of samples of each group:")
table(coldata_subset$condition)
sampleNames_subset =rownames(coldata_subset)
df_subset = df[,sampleNames_subset]
save(coldata_subset,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,'_coldata_subset.Rdata'))



#By adding variables to the design, one can control for additional variation in the counts. For example, if the condition samples are balanced across experimental batches, by including the batch factor to the design, one can increase the sensitivity for finding differences due to condition.

dds <- DESeqDataSetFromMatrix(countData = df_subset,
                              colData = coldata_subset,
                              design= ~ sex + condition)


#Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the reference level
dds$condition <- factor(dds$condition, levels = c(CTRL_group,CASE_group))

dds <- DESeq(dds)
res <- results(dds, name=paste0("condition_",CASE_group,"_vs_",CTRL_group), alpha=0.05)

summary(res)
res.df <- data.frame(res)
res.df$gene <- row.names(res.df)
res.df <- res.df[order(res.df$pvalue),]

print(paste0("The number of DEGs with padj < 0.05 is ",sum(res$padj < 0.05, na.rm=TRUE)))

res.sig.df <- subset(res.df,padj < 0.05 & abs(log2FoldChange) >0.25)
print(paste0("The number of DEGs with padj < 0.05 and log2FC >0.25 is ",nrow(res.sig.df)))

write.table(res.df,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,".pseudobulk_results.raw.tsv"),sep= '\t',row.names=F,col.names=T,quote=F)
write.table(res.sig.df,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,".pseudobulk_results.padj0.05_log2FC0.25.tsv"),sep= '\t',row.names=F,col.names=T,quote=F)
save(dds,res,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,'.dds.res.Rdata'))

# Extracting normalized counts
nomalized_counts <- counts(dds, normalized=TRUE)
save(nomalized_counts,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_nomalized_counts.Rdata"))

# Extracting normalized and transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
#head(assay(vsd), 3)
save(vsd,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_vsd.Rdata"))
save(rld,file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_rld.Rdata"))



#---------------
# gene expression plot; multiple groups
library(reshape2)
library(ggplot2)
library(DESeq2)
library(ggrepel)

gene = 'CCDC3'
celltype = "ValvularEC"

vsd = get(load(paste0(celltype,"_vsd_all_groups.Rdata")))
coldata = get(load(paste0(celltype,'_coldata.Rdata')))

normalizedCountMatrix <- assay(vsd)

normalizedCountMatrix.subset <-normalizedCountMatrix[gene,]
normalizedCountMatrix.subset <- melt(normalizedCountMatrix.subset)
normalizedCountMatrix.subset$sample <- rownames(normalizedCountMatrix.subset)
normalizedCountMatrix.subset$condition <- factor(coldata$condition,levels=c('CTRL','FED','BD','RMVD'))

myplot <- ggplot(normalizedCountMatrix.subset,aes(x=condition,y=value,color=condition,fill=condition))
myplot +  geom_boxplot(outlier.colour = NA, width=0.2,notch=FALSE,alpha=0.5) + geom_point(position = position_jitter(width = 0.1),alpha=1,size=1)+
labs(x='',y='Normalized counts transformed using VST') + theme_bw() + scale_color_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"))  + scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
theme(legend.position="none") +
geom_text_repel(data=normalizedCountMatrix.subset,aes(x=condition,y=value,label=sample),point.padding = unit(0.3,"lines"),max.overlaps=20) +
 ggtitle(paste0(celltype,"  ",gene))
ggsave(file=paste0(celltype,'_',gene,'.pseudoBulk.withlabel.pdf'),width=6.83,height=3.46)

myplot +  geom_boxplot(outlier.colour = NA, width=0.2,notch=FALSE,alpha=0.5) + geom_point(position = position_jitter(width = 0.1),alpha=1,size=1)+
labs(x='',y='Normalized aggregated counts') + theme_bw() + scale_color_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"))  + scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")) +
theme(legend.position="none") +  ggtitle(paste0(celltype,"  ",gene))
ggsave(file=paste0(celltype,'_',gene,'.pseudoBulk.withoutlabel.pdf'),width=6.83,height=3.46)


----------
# gene expression plot; Two groups
library(reshape2)
library(ggplot2)
library(DESeq2)
library(ggrepel)
gene = 'CCDC3'

celltype = "valvularEC"
CTRL_group = 'CTRL'
CASE_group = 'FED'


load(paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_vsd.Rdata"))
load(paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,'_coldata_subset.Rdata'))
#load(paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_nomalized_counts.Rdata"))
#load(paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,"_rld.Rdata"))

normalizedCountMatrix <- assay(vsd)
#normalizedCountMatrix <- assay(rld)
#normalizedCountMatrix <- nomalized_counts
normalizedCountMatrix.subset <-normalizedCountMatrix[gene,]
normalizedCountMatrix.subset <- melt(normalizedCountMatrix.subset)
normalizedCountMatrix.subset$sample <- rownames(normalizedCountMatrix.subset)
normalizedCountMatrix.subset$condition <- factor(coldata_subset$condition,levels=c('CTRL','FED'))

myplot <- ggplot(normalizedCountMatrix.subset,aes(x=condition,y=value,color=condition,fill=condition))
myplot +  geom_boxplot(outlier.colour = NA, width=0.2,notch=FALSE,alpha=0.5) + geom_point(position = position_jitter(width = 0.1),alpha=1,size=1)+
labs(x='',y='Normalized counts transformed using VST') + theme_bw() + scale_color_manual(values=c("#5B74F4","#E82C61"))  + scale_fill_manual(values=c("#5B74F4","#E82C61")) +
theme(legend.position="none") +
geom_text_repel(data=normalizedCountMatrix.subset,aes(x=condition,y=value,label=sample),point.padding = unit(0.3,"lines")) +
 ggtitle(paste0(celltype,"  ",gene))

ggsave(file=paste0(celltype,"_condition_",CASE_group,"_vs_",CTRL_group,'_',gene,'.pseudoBulk.pdf'),width=2.83,height=3.97)


