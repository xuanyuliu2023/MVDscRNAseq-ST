# 🧬 MVDscRNAseq-ST: Single-Cell & Spatial Transcriptomics Pipeline for Human Mitral Valve Disease

> **Repository**: [https://github.com/xuanyuliu2023/MVDscRNAseq-ST](https://github.com/xuanyuliu2023/MVDscRNAseq-ST)  
> **Last Updated**: July 2026  
> **License**: MIT

---

## 📌 Overview

This repository provides a **complete analysis pipeline** for processing **single-cell RNA sequencing (scRNA-seq)** and **spatial transcriptomic data** from **human mitral valve tissues**. The workflow covers quality control, preprocessing, clustering, trajectory inference, gene regulatory network analysis, and co-expression network construction. All scripts are organized in a logical order to ensure reproducibility.

---

## 📁 Repository Structure
MVDscRNAseq-ST/
├── scRNASeq_QC.r # Step 1: scRNA-seq quality control
├── scRNAseq_preprocessing.R # Step 2: scRNA-seq preprocessing
├── SpatialTranscriptomicDataProcessing.r # Step 3: Spatial data preprocessing & visualization
├── pseudobulk_celltype.r # Step 4: Pseudobulk analysis by cell type
├── hdWGCNA.r # Step 5: Co-expression network analysis
├── pySCENIC_VIC.py # Step 6: Regulatory network inference
├── scvelo.py # Step 7: RNA velocity analysis
└── cellrank.py # Step 8: Trajectory inference using CellRank


---

## 🔄 Overall Pipeline Sequence

The analysis follows a **strictly sequential workflow** from raw data to biological interpretation. Each step depends on the output of the previous one:

### Step 1: Quality Control (`scRNASeq_QC.r`)
- **Input**: Raw UMI counts matrix and metadata
- **Process**: Filtering, doublet detection, normalization, PCA, clustering, UMAP visualization
- **Output**: QC-filtered Seurat object

### Step 2: scRNA-seq Preprocessing (`scRNAseq_preprocessing.R`)
- **Input**: QC-filtered Seurat object (from Step 1)
- **Process**: Batch correction, integration, cell cycle regression, normalization, dimensionality reduction, clustering refinement
- **Output**: Preprocessed, integrated AnnData/Seurat object (`.h5ad`)

### Step 3: Spatial Transcriptomics Processing (`SpatialTranscriptomicDataProcessing.r`)
- **Input**: Spatial transcriptomics data (e.g., Visium or Slide-seq) and preprocessed scRNA-seq reference (from Step 2)
- **Process**: Spot-level QC, spatial normalization, deconvolution using scRNA-seq reference, spatial clustering
- **Output**: Spatial Anndata/RDS object

### Step 4: Pseudobulk Analysis (`pseudobulk_celltype.r`)
- **Input**: Preprocessed scRNA-seq data (from Step 2)
- **Process**: Aggregate counts per cell type, differential expression analysis, cell-type-specific marker identification
- **Output**: Pseudobulk expression matrix and DEG tables

### Step 5: Co-expression Network Analysis (`hdWGCNA.r`)
- **Input**: Pseudobulk expression matrix (from Step 4)
- **Process**: Weighted gene co-expression network analysis per cell type, module detection, module-trait correlation
- **Output**: Gene modules and hub genes

### Step 6: Regulatory Network Inference (`pySCENIC_VIC.py`)
- **Input**: Preprocessed scRNA-seq expression matrix (from Step 2)
- **Process**: Transcription factor motif analysis, regulon activity scoring, cell-type-specific regulon identification
- **Output**: Regulon activity matrix and regulatory network

### Step 7: RNA Velocity (`scvelo.py`)
- **Input**: Spliced/unspliced count matrices (from Step 2)
- **Process**: Velocity estimation, latent time inference, velocity stream visualization
- **Output**: Velocity vectors and cell state transition graph

### Step 8: Trajectory Inference (`cellrank.py`)
- **Input**: Velocity results (from Step 7) and preprocessed scRNA-seq data (from Step 2)
- **Process**: Compute transition probabilities, identify root and terminal states, fate mapping, pseudotime estimation
- **Output**: CellRank object with fate probabilities and pseudotime

---

## ⚙️ Prerequisites

- **R ≥ 4.2** with packages: `Seurat`, `SpatialExperiment`, `hdWGCNA`, `monocle3`
- **Python ≥ 3.8** with packages: `scanpy`, `cellrank`, `scvelo`, `pyscenic`

---

## 📊 Key Outputs

- **Step 1**: QC-filtered Seurat object, QC plots
- **Step 2**: Integrated, batch-corrected AnnData object, refined clusters
- **Step 3**: Spatial feature maps, deconvolution proportions
- **Step 4**: Pseudobulk DEGs, cell-type-specific markers
- **Step 5**: Gene modules, module eigengenes, hub genes
- **Step 6**: Regulon activity scores, regulatory network visualizations
- **Step 7**: Velocity streams, latent time estimates
- **Step 8**: Cell fate probabilities, pseudotime trajectories

---

