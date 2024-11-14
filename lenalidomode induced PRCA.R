library(SeuratWrappers)
library(stringr) 
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(liger)
library(cowplot)
library(limma)
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
####PART1
####upload all sample

bmncHD4_file <- str_c(workpath, './20221118/SC1072_RSEC_MolsPerCell.csv')
bmncHD2_file <- str_c(workpath, "./20221118/SC1402_SampleTag10_hs_RSEC_MolsPerCell.csv")
bmncHD3_file <- str_c(workpath, "./20221118/SC1402_SampleTag11_hs_RSEC_MolsPerCell.csv")
bmncHD4_file <- str_c(workpath, './20221118/SC1402_SampleTag09_hs_RSEC_MolsPerCell.csv')
bmncPRCA_P3_file <- str_c(workpath, "./20221118/PRCA_PRCA_P3_RSEC_MolsPerCell.csv")
bmncPRCA_P4_file <- str_c(workpath, "./20221118/PRCA_PRCA_P4_RSEC_MolsPerCell.csv")


##HD1
bmncHD1<- read.table(bmncHD1_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncHD1)=bmncHD1[,1]
bmncHD1 <- bmncHD1[,-1]
bmncHD1 <- t(bmncHD1)

bmncHD1 <- as.matrix(bmncHD1)
bmncHD1 <- as(bmncHD1,"sparseMatrix")
bmncHD1.object <- CreateSeuratObject(counts = bmncHD1, project = "HD1", min.cells = 3, min.features = 200)
bmncHD1.object
bmncHD1.object[['percent_mt']] <- PercentageFeatureSet(bmncHD1.object,pattern = '^MT.')
bmncHD1.object[['id']] <- 'bmncHD1'
head(bmncHD1.object)
head(bmncHD1.object@meta.data, 5)
VlnPlot(bmncHD1.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncHD1.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncHD1.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncHD1.object <- subset(bmncHD1.object, subset = nFeature_RNA > 300 & nFeature_RNA < 4300 & percent_mt < 40)
bmncHD1.object <- NormalizeData(bmncHD1.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncHD1.object <- FindVariableFeatures(bmncHD1.object, selection.method = "vst", nfeatures = 2000)

##HD2
bmncHD2<- read.table(bmncHD2_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncHD2)=bmncHD2[,1]
bmncHD2 <- bmncHD2[,-1]
bmncHD2 <- t(bmncHD2)
bmncHD2 <- as.matrix(bmncHD2)
bmncHD2 <- as(bmncHD2,"sparseMatrix")
bmncHD2.object <- CreateSeuratObject(counts = bmncHD2, project = "HD4", min.cells = 3, min.features = 200)
bmncHD2.object
bmncHD2.object[['percent_mt']] <- PercentageFeatureSet(bmncHD2.object,pattern = '^MT.')
bmncHD2.object[['id']] <- 'bmncHD2'
head(bmncHD2.object)
head(bmncHD2.object@meta.data, 5)
VlnPlot(bmncHD2.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncHD2.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncHD2.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncHD2.object <- subset(bmncHD2.object, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent_mt < 40)
bmncHD2.object <- NormalizeData(bmncHD2.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncHD2.object <- FindVariableFeatures(bmncHD2.object, selection.method = "vst", nfeatures = 2000)


##HD3
bmncHD3<- read.table(bmncHD3_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncHD3)=bmncHD3[,1]
bmncHD3 <- bmncHD3[,-1]
bmncHD3 <- t(bmncHD3)
bmncHD3 <- as.matrix(bmncHD3)
bmncHD3 <- as(bmncHD3,"sparseMatrix")
bmncHD3.object <- CreateSeuratObject(counts = bmncHD3, project = "HD4", min.cells = 3, min.features = 200)
bmncHD3.object
#QC
##HD4QC
bmncHD3.object[['percent_mt']] <- PercentageFeatureSet(bmncHD3.object,pattern = '^MT.')
bmncHD3.object[['id']] <- 'bmncHD3'
head(bmncHD3.object)
head(bmncHD3.object@meta.data, 5)
VlnPlot(bmncHD3.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncHD3.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncHD3.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncHD3.object <- subset(bmncHD3.object, subset = nFeature_RNA > 300 & nFeature_RNA < 3500 & percent_mt < 40)
bmncHD3.object <- NormalizeData(bmncHD3.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncHD3.object <- FindVariableFeatures(bmncHD3.object, selection.method = "vst", nfeatures = 2000)

##HD4
bmncHD4<- read.table(bmncHD4_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncHD4)=bmncHD4[,1]
bmncHD4 <- bmncHD4[,-1]
bmncHD4 <- t(bmncHD4)

bmncHD4 <- as.matrix(bmncHD4)
bmncHD4 <- as(bmncHD4,"sparseMatrix")
bmncHD4.object <- CreateSeuratObject(counts = bmncHD4, project = "HD4", min.cells = 3, min.features = 200)
bmncHD4.object
bmncHD4.object[['percent_mt']] <- PercentageFeatureSet(bmncHD4.object,pattern = '^MT.')
bmncHD4.object[['id']] <- 'bmncHD4'
head(bmncHD4.object)
head(bmncHD4.object@meta.data, 5)
VlnPlot(bmncHD4.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncHD4.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncHD4.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncHD4.object <- subset(bmncHD4.object, subset = nFeature_RNA > 300 & nFeature_RNA < 4300 & percent_mt < 40)
bmncHD4.object <- NormalizeData(bmncHD4.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncHD4.object <- FindVariableFeatures(bmncHD4.object, selection.method = "vst", nfeatures = 2000)

##PRCAP3
bmncPRCAP3<- read.table(bmncPRCAP3_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncPRCAP3)=bmncPRCAP3[,1]
bmncPRCAP3 <- bmncPRCAP3[,-1]
bmncPRCAP3 <- t(bmncPRCAP3)

bmncPRCAP3 <- as.matrix(bmncPRCAP3)
bmncPRCAP3 <- as(bmncPRCAP3,"sparseMatrix")
bmncPRCAP3.object <- CreateSeuratObject(counts = bmncPRCAP3, project = "PRCA_P3", min.cells = 3, min.features = 200)
bmncPRCAP3.object
bmncPRCAP3.object[['percent_mt']] <- PercentageFeatureSet(bmncPRCAP3.object,pattern = '^MT.')
bmncPRCAP3.object[['id']] <- 'bmncPRCAP3'
head(bmncPRCAP3.object)
head(bmncPRCAP3.object@meta.data, 5)
VlnPlot(bmncPRCAP3.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncPRCAP3.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncPRCAP3.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncPRCAP3.object <- subset(bmncPRCAP3.object, subset = nFeature_RNA > 300 & nFeature_RNA < 4300 & percent_mt < 40)
bmncPRCAP3.object <- NormalizeData(bmncPRCAP3.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncPRCAP3.object <- FindVariableFeatures(bmncPRCAP3.object, selection.method = "vst", nfeatures = 2000)

##PRCAP4
bmncPRCAP4<- read.table(bmncPRCAP4_file, sep=",",header=TRUE,stringsAsFactors = FALSE,check.names = T)
rownames(bmncPRCAP4)=bmncPRCAP4[,1]
bmncPRCAP4 <- bmncPRCAP4[,-1]
bmncPRCAP4 <- t(bmncPRCAP4)

bmncPRCAP4 <- as.matrix(bmncPRCAP4)
bmncPRCAP4 <- as(bmncPRCAP4,"sparseMatrix")
bmncPRCAP4.object <- CreateSeuratObject(counts = bmncPRCAP4, project = "PRCA_P4", min.cells = 3, min.features = 200)
bmncPRCAP4.object
bmncPRCAP4.object[['percent_mt']] <- PercentageFeatureSet(bmncPRCAP4.object,pattern = '^MT.')
bmncPRCAP4.object[['id']] <- 'bmncPRCAP4'
head(bmncPRCAP4.object)
head(bmncPRCAP4.object@meta.data, 5)
VlnPlot(bmncPRCAP4.object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt",  ncol = 3))
plot1 <- FeatureScatter(bmncPRCAP4.object, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(bmncPRCAP4.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bmncPRCAP4.object <- subset(bmncPRCAP4.object, subset = nFeature_RNA > 300 & nFeature_RNA < 4300 & percent_mt < 40)
bmncPRCAP4.object <- NormalizeData(bmncPRCAP4.object, normalization.method = "LogNormalize", scale.factor = 10000)
bmncPRCAP4.object <- FindVariableFeatures(bmncPRCAP4.object, selection.method = "vst", nfeatures = 2000)
#integrat all sample
merge.object <- merge(bmncPRCAP4.object,y = bmncPRCAP3.object)
merge.object <- merge(merge.object,y = bmncHD4.object)
merge.object <- merge(merge.object,y = bmncHD3.object)
merge.object <- merge(merge.object,y = bmncHD2.object)
merge.object <- merge(merge.object,y = bmncHD1.object)

# integrat analysis
merge.object <- merge(bmncHD1.object,y = bmncHD2.object)
merge.object <- merge(some,y = merge.object)
merge.object <- NormalizeData(merge.object,normalization.method = "LogNormalize",scale.factor = 10000)
merge.object <- FindVariableFeatures(merge.object,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(merge.object), 10)
plot1 <- VariableFeaturePlot(merge.object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
CombinePlots(plots = list(plot1+plot2))
all.genes <- rownames(merge.object)
merge.object <- ScaleData(merge.object, features = all.genes)
merge.object <- RunPCA(merge.object, features = VariableFeatures(object = merge.object))
print(merge.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merge.object, dims = 1:2, reduction = "pca")
DimPlot(merge.object, reduction = "pca",split.by = 'id')
DimHeatmap(merge.object, dims = 1:15, cells = 500, balanced = TRUE)
merge.object <- JackStraw(merge.object, num.replicate = 100)
merge.object <- ScoreJackStraw(merge.object, dims = 1:20)
JackStrawPlot(merge.object, dims = 1:20)
ElbowPlot(merge.object)
merge.object <- FindNeighbors(merge.object, dims = 1:20)
merge.object <- FindClusters(merge.object, resolution = 0.75)
head(Idents(merge.object), 5) 
set.seed(123)
merge.object <- RunTSNE(merge.object, dims = 1:20)
DimPlot(merge.object, reduction = "tsne",label = T)
DimPlot(merge.object, reduction = "tsne", label = TRUE,split.by = 'id')
merge.object <- RunUMAP(merge.object, dims = 1:20)
DimPlot(merge.object, reduction = "umap")
DimPlot(merge.object, reduction = "umap", label = TRUE,split.by = 'id')
DimPlot(merge.object, reduction = "umap", label = TRUE,split.by = 'id')




###annotation
merge.object@meta.data$group <- 'HD'
merge.object@meta.data[merge.object@meta.data$id == 'PRCA_P3',"group"] <- 'PRCA_P3'
merge.object@meta.data[merge.object@meta.data$id == 'PRCA_P4',"group"] <- 'PRCA_P4'

A <- FindAllMarkers(merge.object)
write.csv(A,'./20221118/ALLMARKER.csv')

FeaturePlot(merge.object,features = c('HBB','GYPA','TFRC','CD19','CD79A'),split.by = 'group',label = T)
FeaturePlot(merge.object,features = c('FCGR3A','CD3D','CD4','CD8A'),split.by = 'group',label = T)
FeaturePlot(merge.object,features = c('CD14','CD19','CD34'),split.by = 'group',label = T)
merge.object@meta.data$seurat_annotation <- 'undefined'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '14',"seurat_annotation"] <- 'plasma cell'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '10',"seurat_annotation"] <- 'ncMono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '7',"seurat_annotation"] <- 'ery'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '12',"seurat_annotation"] <- 'ery'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '11',"seurat_annotation"] <- 'HSC'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '17',"seurat_annotation"] <- 'pDC'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '9',"seurat_annotation"] <- 'cDC'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '8',"seurat_annotation"] <- 'premono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '1',"seurat_annotation"] <- 'cMono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '4',"seurat_annotation"] <- 'cMono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '6',"seurat_annotation"] <- 'cMono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '20',"seurat_annotation"] <- 'cMono'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '21',"seurat_annotation"] <- 'B'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '3',"seurat_annotation"] <- 'B'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '15',"seurat_annotation"] <- 'B'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '19',"seurat_annotation"] <- 'B'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '0',"seurat_annotation"] <- 'CD4T'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '13',"seurat_annotation"] <- 'CD4T'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '2',"seurat_annotation"] <- 'CD8T'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '5',"seurat_annotation"] <- 'NK'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '16',"seurat_annotation"] <- 'B'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '16',"seurat_annotation"] <- 'proB'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '16',"seurat_annotation"] <- 'proB'
merge.object@meta.data[merge.object@meta.data$seurat_clusters == '23',"seurat_annotation"] <- 'eosinophil'
highCells=colnames(subset(x = merge.object, subset = C1QA > 0.2, slot = 'counts'))
highORlow=ifelse(colnames(merge.object) %in% highCells,'high','low')
table(highORlow)
merge.object@meta.data$C1QAhighORlow=highORlow
merge.object@meta.data[merge.object@meta.data$C1QAhighORlow == 'high' & merge.object@meta.data$seurat_clusters == '4',"seurat_annotation"] <- 'ncMonocomp'
merge.object@meta.data[merge.object@meta.data$C1QAhighORlow == 'high' & merge.object@meta.data$seurat_clusters == '10',"seurat_annotation"] <- 'ncMonocomp'
merge.object@meta.data[merge.object@meta.data$C1QAhighORlow == 'high' & merge.object@meta.data$seurat_clusters == '9',"seurat_annotation"] <- 'ncMonocomp'
merge.object@meta.data[merge.object@meta.data$C1QAhighORlow == 'high' & merge.object@meta.data$seurat_clusters == '2',"seurat_annotation"] <- 'ncMonocomp'


##find diPRCA_P4erential genes
LIST <- SplitObject(merge.object,split.by = 'seurat_annotation')

OTHERmarker <- FindMarkers(LIST$ery,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ERY PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$premono,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/premono PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$HSC,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/hsc PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$GMP,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/gmp PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$cDC,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cdc PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$mono,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/mono PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$`plasma cell`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/plasmacell PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$ncMcomp,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ncm comp PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$CD8T,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cd8t PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$cell1,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cell1 PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$ncM,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ncM PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$proB,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/PROB comp PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$pDC,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/PDC PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$CD4T,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/CD4T PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$NK,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/NK PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$B,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/B comp PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$earlyB,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/earlyB PRCA_P4VSHD.csv')
OTHERmarker <- FindMarkers(LIST$NaiveB,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/NaiveB PRCA_P4VSHD.csv')




OTHERmarker <- FindMarkers(LIST$ery,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ERY PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$premono,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/premono PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$HSC,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/hsc PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$GMP,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/gmp PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$cDC,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cdc PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$mono,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/mono PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$`plasma cell`,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/plasmacell PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$ncMcomp,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ncm comp PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$CD8T,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cd8t PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$cell1,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/cell1 PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$ncM,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/ncM PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$proB,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/PROB comp PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$pDC,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/PDC PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$CD4T,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/CD4T PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$NK,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/NK PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$B,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/B comp PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$earlyB,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/earlyB PRCA_P3VSHD.csv')
OTHERmarker <- FindMarkers(LIST$NaiveB,ident.1 = 'PRCA_P3',ident.2 = 'HD',group.by = 'group',min.pct = 0.25)
write.csv(OTHERmarker,file= './20221118/NaiveB PRCA_P3VSHD.csv')


########PART2
#ERY analysis
ery <- subset(merge.object,seurat_annotation %in% 'ery')
DimPlot(ery,group.by = 'id',reduction = 'tsne')
DimPlot(ery,label = T,split.by = 'id',reduction = 'tsne')
FeaturePlot(ery,features = c('TFRC','HBA1','HBB','CA1'),ncol = 2,cols = c('#16168B',"#3f9baf","#7bc344","#f3e407" ,"#ed2c24" ,"#93181e"),pt.size = 0.1)
FeaturePlot(ery,features = c('GYPA','GATA1','SLC4A1','ALAS2'),ncol = 2,cols = c('#16168B',"#7bc344","#f3e407",'#F2711B' ,"#ed2c24" ,"#93181e"),pt.size = 0.1)
FeaturePlot(ery,features = c('GATA2','KIT','PLEK'),cols = c('#16168B',"#3f9baf","#7bc344","#f3e407" ,"#ed2c24" ,"#93181e"),ncol = 2,pt.size = 0.1)
FeaturePlot(ery,features = c('CD34','STAT5A','SOX4','ENO1'),cols = c('#16168B',"#7bc344","#f3e407" ,"#ed2c24" ,"#93181e"),ncol = 2,pt.size = 0.1)
DimPlot(ery,group.by = 'group',cols = c('#PRCA_P4000000','#PRCA_P4000000','#8D0000'),pt.size = 0.3)
DimPlot(ery,group.by = 'group',cols = c('#PRCA_P4000000','#E63636','#PRCA_P4000000'),pt.size = 0.3)
DimPlot(ery,group.by = 'group',cols = c('black','black','black'),pt.size = 0.3)

highCells=colnames(subset(x = ery, subset = GYPA > 0.1, slot = 'counts'))
highORlow=ifelse(colnames(ery) %in% highCells,'high','low')
table(highORlow)
ery@meta.data$CD235highORlow=highORlow

ery@meta.data$state <- 'early'
ery@meta.data[ery@meta.data$CD235highORlow == 'high',"state"] <- 'late'
VlnPlot(ery,features = c('TFRC','GYPA','HBB','GATA1','HLA.A'),split.by = 'group',group.by = 'state')
ERYLIST <- SplitObject(ery,split.by = 'state')
table(ERYLIST$early$id)
table(ERYLIST$late$id)

OTHERmarker <- FindMarkers(ery,ident.1 = 'early',group.by = 'state')
write.csv(OTHERmarker,'./20221118/early marker.csv')
erylist <- SplitObject(ery,split.by = 'state')
VlnPlot(ery,features = c('CCNE1','CDC27','CDK1','CENPE','CENPU','MCM4','MCM10','PCNA','POLE2','TFDP1','TOP2','HMMR','TPX2','CENPF','TOP2A'),group.by = 'state',split.by = 'group',pt.size = 0,cols = c('#2E8B57','#483D8B','#9370DB'),ncol = 5)

OTHERmarker <- FindMarkers(erylist$early,ident.1 ='PRCA_P4',ident.2 = 'HD',group.by = 'group')
write.csv(OTHERmarker,'./20221118/earlyery PRCA_P4vshd.csv')
OTHERmarker <- FindMarkers(erylist$late,ident.1 ='PRCA_P4',ident.2 = 'HD',group.by = 'group')
write.csv(OTHERmarker,'./20221118/late ery PRCA_P4vshd.csv')

OTHERmarker <- FindMarkers(erylist$early,ident.1 ='HD',group.by = 'group')
write.csv(OTHERmarker,'./20221118/earlyery HDVSPRCA.csv')
OTHERmarker <- FindMarkers(erylist$late,ident.1 ='HD',ident.2 = 'PRCA_P4',group.by = 'group')
write.csv(OTHERmarker,'./20221118/late ery HDVSPRCA_P4.csv')
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(R.utils)
options(clusterProfiler.download.method = "wininet")
R.utils::setOption("clusterProfiler.download.method",'auto')
genelist_input <- read.table(file = './20221118/ery early PRCA_p4vshd up.txt', header = T, sep='\t')
genelist_input <- read.table(file = './20221118/ery late PRCA_p4vshd up.txt', header = T, sep='\t')
genelist_input <- read.table(file = './20221118/hspc PRCA_p4vshd up.txt', header = T, sep='\t')
genelist_input <- read.table(file = './20221118/ery early PRCA_p4vshd down.txt', header = T, sep='\t')
genelist_input <- read.table(file = './20221118/ery late PRCA_p4vshd down.txt', header = T, sep='\t')
genelist_input <- read.table(file = './20221118/hspc PRCA_p4vshd down.txt', header = T, sep='\t')

genename <- as.character(genelist_input[,1]) 
gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
gene=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
head(gene_map)
head(genelist_input)
colnames(gene_map)[1]<-"Gene"
colnames(genelist_input)[1]<-"Gene"
aaa<-inner_join(gene_map,genelist_input,by = "Gene")
head(aaa)
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)

geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoPRCA_P4=1)
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoPRCA_P4=1)
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoPRCA_P4=1)

write.table (Go_gseresult, file ="./20221118/ery early PRCA_p4vshd up GO.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="./20221118/ery early PRCA_p4vshd up KEGG .csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="./20221118/ery early PRCA_p4vshd up GSEA.csv", sep =",", row.names =TRUE)

ery <- FindNeighbors(ery, dims = 1:20)
ery <- FindClusters(ery, resolution = 0.75)#?ı?resolution??Ӱ?????ֵ?Ⱥ????��??0.1-1 Խ??ȺԽ??
head(Idents(ery), 5) 
set.seed(123)
ery <- RunTSNE(ery, dims = 1:20)
DimPlot(ery, reduction = "tsne",group.by = 'id',label = T)
DimPlot(ery, reduction = "tsne", label = TRUE,split.by = 'id')
ery <- RunUMAP(ery, dims = 1:20)
DimPlot(ery, reduction = "umap",label = T)
DimPlot(ery, reduction = "umap", label = TRUE,split.by = 'id',group.by = 'seurat_annotation')
DimPlot(ery, reduction = "umap", label = TRUE,split.by = 'id')
VlnPlot(ery,features = c('CD34','GATA2','CD38','CSF2RB','GPI','ANK1','ATP5IF1','AHSP','TFRC','GATA1','CENPE','CCNB1','H4C3','HMGB2','TUBB','ALAS2','GYPA','HBD','HBB','CA1'),pt.size = 0)
FeaturePlot(ery,features = c('CD34','GATA2','CD38','CSF2RB','GPI','ANK1','ATP5IF1','AHSP','TFRC','GATA1','CENPE','CCNB1','H4C3','HMGB2','TUBB','ALAS2','GYPA','HBD','HBB','CA1'),label = T)
head(ery$erystate)
ery@meta.data$erystate <- '3-Poly/Ortho'

ery@meta.data[ery@meta.data$seurat_clusters == '6',"erystate"] <- '1-MEP/BFU/CFU'
ery@meta.data[ery@meta.data$seurat_clusters == '9',"erystate"] <- '1-MEP/BFU/CFU'
ery@meta.data[ery@meta.data$seurat_clusters == '3',"erystate"] <- '1-MEP/BFU/CFU'
ery@meta.data[ery@meta.data$seurat_clusters == '2',"erystate"] <- '2-Pro/Baso'
ery@meta.data[ery@meta.data$seurat_clusters == '0',"erystate"] <- '2-Pro/Baso'
table(ery$group)
prop.table(table(ery$erystate))
table(ery$erystate, ery$group)
Cellratio <- prop.table(table(ery$erystate, ery$group), margin = 2)#??????????????ͬϸ??Ⱥ????
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Group',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
VlnPlot(ery, features = c('HLA.A','HLA.B','HLA.C','B2M'),group.by = 'erystate',pt.size = 0,split.by = 'group',cols = c('#078992','#f8766d','#c10037'))

VlnPlot(ery,features = c('TAPBP','TAP1','TAP2','PDIA3','ERAP1','PSMB1','PSMB2','PSMB5','PSMB8','PSMB9','HLA.A','HLA.B','HLA.C','B2M'),pt.size = 0,group.by = 'erystate',split.by = 'group',ncol = 5,cols = c('#078992','#f8766d','#c10037'))

#######part3 Tcell

TCELL <- subset(merge.object,seurat_annotation %in% c('CD8T','CD4T'))

TCELL <- NormalizeData(TCELL,normalization.method = "LogNormalize",scale.factor = 10000)
TCELL <- FindVariableFeatures(TCELL,selection.method = 'vst',nfeatures = 2000)
top10 <- head(VariableFeatures(TCELL), 10)
plot1 <- VariableFeaturePlot(TCELL)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
CombinePlots(plots = list(plot1+plot2))
all.genes <- rownames(TCELL)
TCELL <- ScaleData(TCELL, features = all.genes)
TCELL <- RunPCA(TCELL, features = VariableFeatures(object = TCELL))
print(TCELL[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TCELL, dims = 1:2, reduction = "pca")
DimPlot(TCELL, reduction = "pca",split.by = 'id')
DimHeatmap(TCELL, dims = 1:15, cells = 500, balanced = TRUE)
TCELL <- JackStraw(TCELL, num.replicate = 100)
TCELL <- ScoreJackStraw(TCELL, dims = 1:20)
JackStrawPlot(TCELL, dims = 1:20)
ElbowPlot(TCELL)
TCELL <- FindNeighbors(TCELL, dims = 1:20)
TCELL <- FindClusters(TCELL, resolution = 1.5)#?ı?resolution??Ӱ?????ֵ?Ⱥ????��??0.1-1 Խ??ȺԽ??
head(Idents(TCELL), 5) 
set.seed(123)
TCELL <- RunTSNE(TCELL, dims = 1:20)
DimPlot(TCELL, reduction = "tsne",group.by = 'seurat_annotation',label = T)
DimPlot(TCELL, reduction = "tsne", label = TRUE,split.by = 'id')
TCELL <- RunUMAP(TCELL, dims = 1:20)
DimPlot(TCELL, reduction = "umap",label = T)

FeaturePlot(TCELL,features = c('FCGR3A','CD3E','CD4','CD8A'),split.by = 'tgroup',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('FCGR3A','CD3E','CD4','CD8A'),split.by = 'tgroup',ncol = 1)
FeaturePlot(TCELL,features = c('CD4','SELL','LEF1','S100A4','ITGB1'),split.by = 'tgroup',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('NCAM1','FCGR3A'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1)

FeaturePlot(TCELL,features = c('FOXP3','GZMK','GZMB','GZMH'),split.by = 'tgroup',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('FOXP3','GZMK','GZMB','GZMH'),split.by = 'tgroup',ncol = 1)
VlnPlot(TCELL,features = c('FCGR3A','CD8A','GZMK','GZMB','GZMH'),split.by = 'tgroup',ncol = 1)
VlnPlot(TCELL,features = c('FCGR3A'),split.by = 'tgroup',ncol = 1)

TCELL@meta.data$TCELLsubtype <- 'other'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '13',"TCELLsubtype"] <- 'CD4T-Treg'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '7',"TCELLsubtype"] <- 'CD8T-memory'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '9',"TCELLsubtype"] <- 'CD8T-memory'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '20',"TCELLsubtype"] <- 'CD8T-memory'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '10',"TCELLsubtype"] <- 'CD8T-memory'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '6',"TCELLsubtype"] <- 'CD8T-cytotoxic'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '8',"TCELLsubtype"] <- 'CD8T-cytotoxic'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '11',"TCELLsubtype"] <- 'CD8T-cytotoxic'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '18',"TCELLsubtype"] <- 'CD8T-cytotoxic'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '1',"TCELLsubtype"] <- 'NK'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '12',"TCELLsubtype"] <- 'NK'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '0',"TCELLsubtype"] <- 'CD4T-help'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '4',"TCELLsubtype"] <- 'CD4T-help'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '14',"TCELLsubtype"] <- 'CD4T-help'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '20',"TCELLsubtype"] <- 'CD4T-help'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '3',"TCELLsubtype"] <- 'CD4T-naive'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '3',"TCELLsubtype"] <- 'CD4T-naive'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '2',"TCELLsubtype"] <- 'CD4T-naive'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '5',"TCELLsubtype"] <- 'CD4T-naive'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '15',"TCELLsubtype"] <- 'NK2'

highCells=colnames(subset(x = TCELL, subset = FCGR3A > 0.1, slot = 'counts'))
highORlow=ifelse(colnames(TCELL) %in% highCells,'high','low')
table(highORlow)
TCELL@meta.data$CD16highORlow=highORlow
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '8' & TCELL@meta.data$CD16highORlow == 'high',"TCELLsubtype"] <- 'NK'
TCELL@meta.data[TCELL@meta.data$CD16highORlow == 'low'& TCELL@meta.data$TCELLsubtype == 'CD8' ,"TCELLsubtype"] <- 'CD8T-cytotoxic'

highCells=colnames(subset(x = TCELL, subset = CD8A > 0.8, slot = 'counts'))
highORlow=ifelse(colnames(TCELL) %in% highCells,'high','low')
table(highORlow)
TCELL@meta.data$CD8highORlow=highORlow


highCells=colnames(subset(x = TCELL, subset = CD4 > 0.1, slot = 'counts'))
highORlow=ifelse(colnames(TCELL) %in% highCells,'high','low')
table(highORlow)
TCELL@meta.data$CD4highORlow=highORlow


TCELL@meta.data[TCELL@meta.data$seurat_clusters == '2' & TCELL@meta.data$CD8highORlow == 'high',"TCELLsubtype"] <- 'CD8T-naive'
TCELL@meta.data[TCELL@meta.data$seurat_clusters == '3' & TCELL@meta.data$CD16highORlow == 'high',"TCELLsubtype"] <- 'CD8T-naive'
FeaturePlot(TCELL,features = c('FCGR3A','CD3E','CD4','CD8A'),split.by = 'id',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('FCGR3A','CD3E','CD4','CD8A'),split.by = 'TCELLsubtype',group.by = 'id',ncol = 1)
FeaturePlot(TCELL,features = c('CD4','SELL','LEF1','S100A4','ITGB1'),split.by = 'tgroup',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('CD4','SELL','LEF1','S100A4','ITGB1'),split.by = 'tgroup',ncol = 1)
VlnPlot(TCELL,features = c('CD28'),split.by = 'tgroup',ncol = 1)

FeaturePlot(TCELL,features = c('FOXP3','GZMK','GZMB','GZMH'),split.by = 'tgroup',cols = c('GREY' ,"#5599PRCA_P4",'ORANGE',"#DE1F1F"),label = T)
VlnPlot(TCELL,features = c('FOXP3','GZMK','GZMB','GZMH'),split.by = 'tgroup',ncol = 1)
VlnPlot(TCELL,features = c('FCGR3A','CD8A','GZMK','GZMB','GZMH'),split.by = 'tgroup',ncol = 1)
VlnPlot(TCELL,features = c('FCGR3A'),split.by = 'tgroup',ncol = 1)

VlnPlot(TCELL,features = c('CD4'),split.by = 'TCELLsubtype',group.by = 'id',ncol = 1)
VlnPlot(TCELL,features = c('CD8A'),split.by = 'TCELLsubtype',group.by = 'id',ncol = 1)
VlnPlot(TCELL,features = c('SELL'),split.by = 'TCELLsubtype',group.by = 'id',ncol = 1)
VlnPlot(TCELL,features = c('LEF1'),split.by = 'TCELLsubtype',group.by = 'id',ncol = 1)

TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD4T-naive' & TCELL@meta.data$CD8highORlow == 'high',"TCELLsubtype"] <- 'CD8Tnaive'
TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD8Tnaive',"TCELLsubtype" ] <- 'CD8T-naive'
TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD8T-ePRCA_P4ecor' & TCELL@meta.data$CD16highORlow == 'high',"TCELLsubtype"] <- 'NK'
TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD8T-memory' & TCELL@meta.data$id == 'hdsmall' & TCELL@meta.data$CD16highORlow == 'high',"TCELLsubtype"] <- 'NK'
TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD4T-naive' & TCELL@meta.data$id == 'gdg_lh' & TCELL@meta.data$CD4highORlow == 'low',"TCELLsubtype"] <- 'CD8-naive'
TCELL@meta.data[TCELL@meta.data$TCELLsubtype == 'CD8T-ePRCA_P4ecor', "TCELLsubtype"] <- 'CD8T-cytotoxic'
TCELLlist <- SplitObject(TCELL,split.by = 'TCELLsubtype')
DimPlot(TCELL,group.by = 'TCELLsubtype',label = T)
DimPlot(TCELL,group.by = 'TCELLsubtype')
DimPlot(TCELL,group.by = 'id')

OTHERmarker <- FindMarkers(TCELLlist$`CD8T-cytotoxic`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/PRCA_P4vdhd cd8 cytotxic.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD8T-cytotoxic`,ident.1 = 'PRCA',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/lhvdhd cd8 cytotxic.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD8T-cytotoxic`,ident.1 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/hd  vs other cd8 cytotxic.csv')

OTHERmarker <- FindMarkers(TCELLlist$`CD8T-memory`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/PRCA_P4vdhd cd8 memory.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD8T-memory`,ident.1 = 'PRCA',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/lhvdhd cd8 memory.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD8T-memory`,ident.1 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/hd  vs other cd8 memory.csv')

OTHERmarker <- FindMarkers(TCELLlist$`CD4T-naive`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/PRCA_P4vdhd CD4T-naive.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-naive`,ident.1 = 'PRCA',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/lhvdhd CD4T-naive.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-naive`,ident.1 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/hd  vs other CD4T-naive.csv')

OTHERmarker <- FindMarkers(TCELLlist$`CD4T-help`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/PRCA_P4vdhd CD4T-help.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-help`,ident.1 = 'PRCA',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/lhvdhd CD4T-help.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-help`,ident.1 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/hd  vs other CD4T-help.csv')

OTHERmarker <- FindMarkers(TCELLlist$`CD4T-Treg`,ident.1 = 'PRCA_P4',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/PRCA_P4vdhd CD4T-Treg.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-Treg`,ident.1 = 'PRCA',ident.2 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/lhvdhd CD4T-Treg.csv')
OTHERmarker <- FindMarkers(TCELLlist$`CD4T-Treg`,ident.1 = 'HD',group.by = 'tgroup')
write.csv(OTHERmarker,'./20221118/TCELL/hd  vs other CD4T-Treg.csv')



VlnPlot(TCELL,features = c('CCR7','TCF7','LEF1','SELL'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('NKG7','CCL4','CST7','PRF1','GZMA','GZMB','IFNG','CCL3'),pt.size = 0,split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 2,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('PRF1','IFNG','GNLY','NKG7','GZMB','GZMA','GZMH','KLRK1','KLRB1','KLRD1','CTSW','CST7'),pt.size = 0,split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 2,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))

VlnPlot(TCELL,features = c('PDCD1','TIGIT','LAG3','HAVCR2','CTLA4'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 2,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))

VlnPlot(TCELL,features = c('MKI67','TYMS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 2,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))

VlnPlot(cd8tx,features = c('CD69','NKG7'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 2,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('CD69','GZMB','IFNG','PRF1'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0.005)
VlnPlot(TCELL,features = c('CD69'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('GIMAP7','FGFBP2','LGALS1','PRF1','CD52','GZMH','GNLY'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0.005)
VlnPlot(TCELL,features = c('FASLG','FAS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0.005)
VlnPlot(TCELL,features = c('IL2','TNF'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0.005)
VlnPlot(TCELL,features = c('CMC1','CCL5','CST7','CCL4'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('NKG7','GZMB','GZMA','KLRD1'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('GZMB','GZMA','GZMH','GZMM','GZMK'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('FCGR3A'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('KLRD1','KLRB1','TIGIT'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('LAG3','IL2RA','CD28','ICOS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('LAG3','IL2RA','CD28','ICOS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('GZMB','GZMA','GZMH','GZMM','GZMK'),group.by = 'tgroup',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('CD28'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('NFKB1','NFKBIA','SQSTM1'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)

VlnPlot(merge.object,features = c('NFKB1','NFKBIA','SQSTM1'),split.by = 'group',group.by = 'seurat_annotation',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(merge.object,features = c('NFKB1'),split.by = 'group',group.by = 'seurat_annotation',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
table(merge.object$group)


VlnPlot(TCELL,features = c('LAG3','IL2RA','CD28','ICOS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('LAG3','IL2RA','CD28','ICOS'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'))
VlnPlot(TCELL,features = c('GZMB','GZMA','GZMH','GZMM','GZMK'),group.by = 'tgroup',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('CD28'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)
VlnPlot(TCELL,features = c('NFKB1','NFKBIA','SQSTM1'),split.by = 'tgroup',group.by = 'TCELLsubtype',ncol = 1,cols = c('#2B4ACA','#PRCA_P47777','#D11414'),pt.size = 0)


#####PART4
##ITALK 分析
library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)
library(SeuratWrappers)
library(stringr) 
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(liger)
library(cowplot)
library(limma)

iTalk_data <- as.data.frame(t(merge.object@assays$RNA@counts))

iTalk_data$cell_type <- merge.object@meta.data$seurat_annotation

iTalk_data$compare_group <- merge.object@meta.data$group

unique(iTalk_data$cell_type)# "cd56_nk" "cd14_monocytes" "b_cells" "cytotoxic_t" "regulatory_t" "memory_t" "naive_t"
unique(iTalk_data$compare_group)
# "group1" "group2" "group3"
my10colors <-c('#ff9a36','#FFAAC6','#E1B675','#FFD700','#8FBC8F','#C7EB52','#D45E5E','#B54CB9','#AC8F14','#F7E1C7','#D8D822','#7E95F2','#25aff5','#1E7BB0','#DCC4FA','#C66B03','grey','#9EF6FF','#6B8E23','grey','pink','orange','blue')
cell_types <- c('CD8T','CD4T',"B", "cMono", "cDC","HSC", "plasma cell","NK","GMP","ery",'CD8T','ncMonocomp','pDC','ncM')

highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
write.csv(highly_exprs_genes,file = './20221118/ITALK/PRCAFFgene.csv')
highly_exprs_genes <- read.csv('./20221118/ITALK/PRCAFFgene.csv')
rownames(highly_exprs_genes)<-highly_exprs_genes[,1]

highly_exprs_genes <- highly_exprs_genes[,-1]

comm_list <- c('growth factor','other','cytokine','checkpoint')

head(comm_list)

cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)
cell_col
cell_types
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res,file= str_c('./20221118/ITALK/all.csv'))

LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)


iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),]
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(iTalk_res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:20],link.arr.width=iTalk_res$cell_to_mean_exprs[1:20])

a <- cell_types[-1]
for (i in seq_along(a)) { 
  j=i 
  
  cell1= cell_types[i]
  deg_t<-DEG(iTalk_data %>% filter(cell_type== str_c(cell1)),method='Wilcox',contrast=c('PRCA_P4','HD'))
  while (j<=16) {
    cell2= cell_type[j]
    deg_nk<-DEG(iTalk_data %>% filter(cell_type==str_c(cell2)),method='Wilcox',contrast=c('PRCA_P4','HD'))
    
    res<-NULL
    res.omit<-NULL
    for(comm_type in comm_list){
      res_cat<-FindLR(deg_nk, deg_t,  datatype='DEG',comm_type=comm_type)
      #res_cat<-FindLR(deg_t,  datatype='DEG',comm_type=comm_type)
      res<-rbind(res,res_cat)
    }
    ʾ
    write.csv(res,file= str_c('E:/20221118/ITALK/',str_c(cell1),'_',str_c(cell2),'wilcox.csv'))
    res.omit <- na.omit(res)
    write.csv(res.omit,file= str_c('E:/20221118/ITALK/',str_c(cell1),'_',str_c(cell2),'omit new wilcox.csv'))
    
    
    j = j +1 
  }
}
    hh <- read.csv('./20221118/ITALK/CD8Tvs ery.csv')
    LRPlot(hh.omit,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
    hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$receptor.on))
    hh$order2=factor(rev(as.integer(rownames(hh))),labels = rev(hh$ligand.from))
    ggplot(hh,aes(x=cell_to,y=order))+
      geom_point(aes(size=`Log10.cell_from_q.value`,
                     color=`cell_from_logFC`))+
      theme_bw()+
      theme(panel.grid = element_blank(),
            axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
      scale_color_gradientn(colours =c('#0000CD','#456DEF','#7C9BFF','#FFA500','#FF8C00','#FF7F50','#FF4500','red'))+
      labs(x=NULL,y=NULL)
    