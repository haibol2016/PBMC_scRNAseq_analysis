## using R 4.1
library(SoupX)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(metap)
library(celldex)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(BiocParallel)
library(DoubletFinder)
library(SummarizedExperiment)

plan("multicore", workers = 4)
set.seed(12345)
options(future.globals.maxSize = 4000 * 1024 ^ 2)

# Load the pbmc dataset

out <- "05152022.10xPBMC.scRNA-seq.out"
if (!dir.exists(out)) {
    dir.create(out)
}

library_id <- c("pbmc_10k_v3", "pbmc_10k_v3.1", 
                "pbmc_1k_v3",  "pbmc_5k_v3", "pbmc_5k_v3_ng")
cellranger_outs <- file.path("../", library_id, "outs")
names(cellranger_outs) <- gsub("_", "", library_id)

sc <- mapply(function(.x, .y) {
    sc <- load10X(.x)
    sc <- autoEstCont(sc, doPlot = FALSE)
    out <- adjustCounts(sc)
    cntSoggy <- rowSums(sc$toc > 0)
    cntStrained <- rowSums(out > 0)
    print(tail(sort((
        cntSoggy - cntStrained
    ) / cntSoggy), n = 10))
    print(tail(sort(
        rowSums(sc$toc > out) / rowSums(sc$toc > 0)
    ), n = 10))
    colnames(out) <- paste0(.y, "_", colnames(out))
    out
}, cellranger_outs, names(cellranger_outs))

sc = do.call(cbind, sc)
sc = CreateSeuratObject(
    sc,
    project = names(cellranger_outs),
    min.cells = 3,
    min.features = 200
)

saveRDS(sc, file = file.path(out, "SoupX.corrected.scRNA-seq.Seurat.RDS"))
sc <- readRDS(file.path(out, "SoupX.corrected.scRNA-seq.Seurat.RDS"))
# The [[ operator can add columns to object metadata.
# This is a great place to stash QC stats
sc[["percent_mt"]] <- PercentageFeatureSet(sc,
                                           pattern = "MT-.+-ENS.+")
sc[["percent_ribo"]] <- PercentageFeatureSet(sc,
                                             "(M)?RP[SL]\\d+([^K]+)?-ENS.+$")
sc[["percent_hb"]] <-  PercentageFeatureSet(sc, "HB[BA]-ENS.+")
sc[["percent_plat"]] <-
    PercentageFeatureSet(sc, "(PECAM1|PF4)-ENS.+$")

# Visualize QC metrics as a violin plot
pdf(file.path(out, "Fig 1.Basic QC plots.pdf"),
    height = 10,
    width = 15)
VlnPlot(
    sc,
    features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent_mt",
        "percent_ribo",
        "percent_hb",
        "percent_plat"
    ),
    split.by = 'orig.ident',
    ncol = 3
)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships,
# but can be used for anything calculated by the object, i.e. columns in object
# metadata, PC scores etc.
pdf(file.path(
    out,
    paste(
        "Fig 2.Scatter plot showing relationship between",
        "UMIcounts and mt_perc.pdf"
    )
),
height = 5, width = 10)
plot1 <- FeatureScatter(sc,
                        feature1 = "nCount_RNA",
                        feature2 = "percent_mt",
                        group.by = 'orig.ident')
plot2 <- FeatureScatter(sc,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA",
                        group.by = 'orig.ident')
plot1 + plot2
dev.off()

## mt perncentage per library
pdf(
    file.path(out, "Fig 3.Ridge plot showing MT percentage per library.pdf"),
    width = 8,
    height = 5
)
RidgePlot(object = sc,
          features = "percent_mt",
          group.by = 'orig.ident',)
dev.off()

# gene expression
selected_c <- WhichCells(sc, expression = nFeature_RNA > 200 &
                             percent_mt < 12.5 & percent_ribo > 5)
sc <- subset(sc, cells = selected_c)

# 25437 genes X 41208 cells
cDat <-  as.matrix(GetAssayData(object = sc, slot = 'counts'))

selected_f <- rownames(sc)[rowSums(matrix(cDat != 0,
                                          nrow = nrow(cDat))) > 
                               ncol(cDat) * 1 / 1000]
sc <- subset(sc, features = selected_f)

# Filter genes
# As the level of expression of mitochondrial and MALAT1 genes are judged as
# mainly technical, it can be wise to remove them from the dataset before
# any further analysis.
sc <-
    sc[!grepl("MALAT1", rownames(sc)) &
           !grepl("MT-.+-ENS.+", rownames(sc)) &
           !grepl("HB[BA]-ENS.+", rownames(sc)),]  # 19070 genes

# cells per library
# table(sc@meta.data$orig.ident)

# pbmc10kv3 pbmc10kv3.1    pbmc1kv2    pbmc1kv3    pbmc4kv2    pbmc5kv3 
# 9921        9061         977         936        4504        4256 
# pbmc5kv3ng    pbmc8kv2 
# 4363        8560 


pdf(
    file.path(out, "Fig 4.Basic QC plots post filtering.pdf"),
    height = 10,
    width = 15
)
VlnPlot(
    sc,
    features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent_mt",
        "percent_ribo",
        "percent_hb",
        "percent_plat"
    ),
    split.by = 'orig.ident',
    ncol = 3
)
dev.off()

s_genes <-  cc.genes$s.genes
s_genes <-
    rownames(sc)[gsub("-ENS.*", "", rownames(sc)) %in% s_genes]

g2m_genes <-  cc.genes$g2m.genes
g2m_genes <-
    rownames(sc)[gsub("-ENS.*", "", rownames(sc)) %in% g2m_genes]


# split the dataset into a list of two seurat objects (stim and CTRL)
pbmc.list <- SplitObject(object = sc, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
pbmc.list <- lapply(
    X = pbmc.list,
    FUN = function(x) {
        x <- NormalizeData(x)
        x <- CellCycleScoring(object = x,
                              g2m.features = g2m_genes,
                              s.features = s_genes)
        x <-
            FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    }
)


pbmc.list <- mapply(function(x, y) {
    x <- ScaleData(
        x,
        vars.to.regress = c("nFeature_RNA", "percent_mt"),
        verbose = FALSE
    )
    x <- RunPCA(x, verbose = FALSE)
    pdf(file.path(
        out,
        paste0("Fig 5.", y,
               " Elbow plot showing PC significance.pdf")
    ),
    height = 10,
    width = 20)
    print(ElbowPlot(x, ndims = 50, reduction = "pca"))
    dev.off()
    x
}, pbmc.list, names(pbmc.list), SIMPLIFY = FALSE)



save.image(file = file.path(out, "0515.RData"))
load(file.path(out, "0515.RData"))


pbmc.list <- mapply(function(x, y) {
    x <- RunUMAP(x,
                 dims = 1:y,
                 verbose = TRUE,
                 reduction = "pca",
    )
}, pbmc.list, c(20, 20, 18, 21, 22), SIMPLIFY = FALSE)


pK <- mapply(function(x, y) {
    sweep.res <- paramSweep_v3(x)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pdf(file.path(out, paste0("Fig.6", y, ".pK.determination.pdf")))
    barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2)
    dev.off()
}, pbmc.list, names(pbmc.list), SIMPLIFY = FALSE)

pbmc.list <-  mapply(function(x, y, z, j) {
    nExp <- round(ncol(x) * y)
    x <- doubletFinder_v3(
        x,
        pN = 0.25,
        pK = z,
        nExp = nExp,
        PCs = 1:j
    )
    x
},
pbmc.list,
c(0.25, 0.19, 0.03, 0.23, 0.19),
c(0.05, 0.05, 0.005, 0.025, 0.025),
c(20, 20, 18, 21,22),
SIMPLIFY = FALSE)


## plot changes
null <- mapply(function(.x, .y) {
    # name of the DF prediction can change, so extract the correct column name.
    DF.name = colnames(.x@meta.data)[grepl("DF.classification",
                                           colnames(.x@meta.data))]
    pdf(file.path(out, paste0("Fig.8.", .y, "doublets.pdf")),
        height = 5,
        width = 10)
    print(cowplot::plot_grid(
        ncol = 2,
        DimPlot(.x, group.by = "orig.ident") + NoAxes(),
        DimPlot(.x, group.by = DF.name) + NoAxes()
    ))
    dev.off()
    pdf(file.path(out, paste0(
        "Fig.7.", .y,
        "doublets.nFeatures_RNA.pdf"
    )),
    height = 5,
    width = 10)
    print(VlnPlot(
        .x,
        features = "nFeature_RNA",
        group.by = DF.name,
        pt.size = 0.1
    ))
    dev.off()
}, pbmc.list, names(pbmc.list), SIMPLIFY = FALSE)

## remove doublets
pbmc.list <- lapply(pbmc.list, function(.x) {
    DF.name = colnames(.x@meta.data)[grepl("DF.classification",
                                           colnames(.x@meta.data))]
    
    .x <- .x[, .x@meta.data[, DF.name] == "Singlet"]
    .x
})

saveRDS(pbmc.list,
        file = file.path(out,
                         "02.Doublets.removed.Seurat.8.objs.RDS"))
## only version 3 
pbmc.list <- readRDS(file = file.path(out,
                                      "02.Doublets.removed.Seurat.8.objs.RDS"))
## re-do integration
pbmc.list <- lapply(
    X = pbmc.list,
    FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x,
                                  selection.method = "vst",
                                  nfeatures = 2000)
    }
)

# select all features across datasets for integration
features <- SelectIntegrationFeatures(object.list = pbmc.list)

# RPCA-based integration runs significantly faster, and also represents a more
# conservative approach where cells in different biological states are less
# likely to 'align' after integration.
pbmc.anchors <- FindIntegrationAnchors(
    object.list = pbmc.list,
    anchor.features = features,
    reduction = "cca"
)
# this command creates an 'integrated' data assay of all genes
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors,
                               features.to.integrate = features)

# specify that we will perform downstream analysis on the corrected
# data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(pbmc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.combined <- ScaleData(
    pbmc.combined,
    features = features,
    vars.to.regress = c("nFeature_RNA",
                        "percent_mt"),
    verbose = FALSE
)

pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

##  18 PC for mesenchyme
pdf(
    file.path(out, 
              "Fig 8.Elbow plot showing PC significance post filtering.pdf"),
    height = 10,
    width = 20
)
ElbowPlot(pbmc.combined, ndims = 50, reduction = "pca")
dev.off()
npc <- 20
future.seed <- 12334
pbmc.combined <- RunUMAP(pbmc.combined,
                         reduction = "pca", dims = 1:npc)
pbmc.combined <- RunTSNE(pbmc.combined,
                         reduction = "pca", dims = 1:npc)
pbmc.combined <- FindNeighbors(pbmc.combined,
                               reduction = "pca", dims = 1:npc)

## reanalysis
resolution <- 0.7
pbmc.combined <- FindClusters(pbmc.combined,
                              random.seed = 0,
                              resolution = resolution)

# Visualization

saveRDS(pbmc.combined, 
        file.path(out, 
                  "03.Doublets.removed.Seurat.integrated.5.objs.RDS"))


pdf(file.path(
    out,
    paste0(
        "Fig.11.UMAP.showing.clusters.by.sample.",
        resolution,
        ".pdf"
    )
),
width = 25, height = 5)
DimPlot(pbmc.combined, reduction = "umap",
        split.by = "orig.ident")

dev.off()
pdf(file.path(
    out,
    paste0(
        "Fig.11.UMAP.showing.clusters.",
        resolution,
        ".pdf"
    )
),
width = 12, height = 5)
p1 <- DimPlot(pbmc.combined, reduction = "umap",
              group.by = "ident")
p2 <- DimPlot(pbmc.combined,
              reduction = "umap",
              label = TRUE,
              label.box = TRUE,
              repel = TRUE)
p1 + p2
dev.off()


# Visualization
pdf(file.path(out, paste0(
    "Fig 12.tSNE showing clusters.",
    resolution, ".pdf"
)), width = 12, height = 5)
p1 <- DimPlot(pbmc.combined, reduction = "tsne",
              group.by = "ident")
p2 <- DimPlot(pbmc.combined,
              reduction = "tsne",
              label = TRUE,
              repel = TRUE)
p1 + p2
dev.off()

pdf(file.path(
    out,
    paste0("Fig 13.tSNE showing clusters by sample.",
           resolution, ".pdf")
),
width = 25, height = 5)
DimPlot(pbmc.combined, reduction = "tsne", split.by = "orig.ident")
dev.off()


## 
ref_se <- MonacoImmuneData(ensembl = TRUE, cell.ont = "all")

counts <- as.matrix(pbmc.combined@assays$RNA@counts)
pbmc_se <- SummarizedExperiment(assays = list(counts = counts),
                                    colData = pbmc.combined@meta.data)
saveRDS(pbmc_se, file = file.path(out, "pbmc_scrnaseq_se.RDS"))

pbmc_se@NAMES <- gsub(".+?-(ENSG\\d+)", "\\1", pbmc_se@NAMES)
row_names <- rownames(pbmc_se@assays@data@listData$counts)

## cell-level annotation
anno_pbmc <- SingleR(test = pbmc_se, ref = ref_se,
                     assay.type.test = "counts",
                     assay.type.ref = "logcounts",
                     clusters = NULL,
                     labels = ref_se$label.fine,
                     BPPARAM = MulticoreParam())


short_name <- c("cMN", "Naive B","h+fhT", "NK", 
                "h+fhT", "Naive CD4+", "Naive CD8+",
                "h+fhT", "h+fhT","cmCD8+", "Treg",
                "h+fhT", "MAIT", "mB", "GDT","mB",
                "mB", "intMN", "teCD4+","te+emCD8+",
                "Prog","GDT", "Pb", "te+emCD8+","mDC",
                "pDC", "ncMN")
names(short_name) <- unique(anno_pbmc$labels)

pbmc.combined[["SingleR.labels"]] <- as.factor(short_name[anno_pbmc$labels])

table(short_name[anno_pbmc$labels])
# cmCD8+     cMN         GDT       h+fhT      intMN       MAIT        mB 
#  189       6346        657       3073        319        661        1121 
# mDC      Naive B   Naive CD4+  Naive CD8+    ncMN        NK         Pb 
#  52        1796       4401        1912        16        1112        14 
# pDC       Prog     te+emCD8+     teCD4+       Treg 
#  9         16         228          51          526 


pdf(file.path(
    out,
    paste0(
        "Fig.9.",
        npc,
        "PCs.UMAP.showing.clusters.by.labels.",
        resolution,
        ".pdf"
    )
),
width = 18, height = 7)
p1 <- DimPlot(
    pbmc.combined,
    reduction = "umap",
    label = TRUE,
    label.box = TRUE,
    repel = TRUE,
    group.by = "SingleR.labels"
)
p2 <- DimPlot(
    pbmc.combined,
    reduction = "umap",
    label.box = TRUE,
    label = TRUE,
    repel = TRUE
)
p1 + p2

dev.off()

pdf(file.path(
    out,
    paste0(
        "Fig.10.",
        npc,
        "PCs.tSNE.showing.clusters.by.labels.",
        resolution,
        ".pdf"
    )
),
width = 18, height = 7)
p1 <- DimPlot(
    pbmc.combined,
    reduction = "tsne",
    label = TRUE,
    label.box = TRUE,
    repel = TRUE,
    group.by = "SingleR.labels"
)
p2 <- DimPlot(
    pbmc.combined,
    reduction = "tsne",
    label = TRUE,
    label.box = TRUE,
    repel = TRUE
)
p1 + p2

dev.off()

counts <- as.matrix(pbmc.combined@assays$RNA@counts)
pbmc_se <- SummarizedExperiment(assays = list(counts = counts),
                                colData = pbmc.combined@meta.data)

saveRDS(pbmc_se, file = file.path(out, "04.annotated.pbmc_scrnaseq_se.RDS"))


# For performing differential expression after integration,
## we switch back to the original data
DefaultAssay(pbmc.combined) <- "RNA"
pbmc.markers <- FindAllMarkers(
    pbmc.combined,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

saveRDS(pbmc.combined,
        file = file.path(out, "03.Doublets.removed.Seurat.integrated.RDS"))

DefaultAssay(pbmc.combined) <- "RNA"


feature2plot_0 <-
    c(  "CD33", "CD86", "FCGR1A",
        "S100A8",
        "S100A9",
        "S100A12",
        "FCGR1B",
        "LYZ",
        "FCGR2A",
        "CST3",
        "LGALS3",
        "CD14",
        "CD68",
        "TLR1",
        "TLR2",
        "TLR4",
        "TLR6",
        "TLR7",
        "TLR8",
        "CCR2",
        "MS4A7",
        "ITGAM",
        "FCER1G",
        "S100A4",
        "CD44",
        "SELL",
        "LTB",
        "CCR7",
        "CD27",
        "CD3D",
        "CD3E",
        "CD3G",
        "IL7R",
        "CD28",
        "CD4",
        "CD8A",
        "CD8B",
        "PF4",
        "PPBP",
        "PTCRA",
        "FOXP3",
        "IL2RA",
        "KLRB1",
        "GZMK",
        "ALOX5AP",
        "GNLY",
        "NKG7",
        "FCGR3A",
        "TRDC",
        "GZMH",
        "GZMB",
        "KLRC2",
        "KLRC1",
        "CD79A",
        "TCL1A",
        "CD19",
        "TLR10",
        "IGHA1",
        "IGHG1",
        "TSPAN13",
        "HLA-DRA",
        "HLA-DOA",
        "HLA-DQA1",
        "MS4A1",
        "CD1C",
        "FCER1A",
        "IL1B",
        "CD38",
        "TYMS"
    )

feature2plot <- rownames(pbmc.combined)[gsub("-ENS.+", "",
                                             rownames(pbmc.combined)) %in%  feature2plot_0]
names(feature2plot) <-
    gsub("-ENS.+", "", feature2plot, perl = TRUE)
feature2plot <- feature2plot[feature2plot_0]

feature2plot <- unname(feature2plot)
pdf(file.path(
    out,
    paste0("Fig 15.Feature.plots.on.TSNE-3.",
           resolution, ".pdf")
), height = 85, width = 19)
FeaturePlot(
    pbmc.combined,
    features = feature2plot,
    order = TRUE,
    reduction = "tsne",
    ncol = 4,
    min.cutoff = "q5"
)
dev.off()

## violin plots
pdf(file.path(out, paste0(
    "Fig 16.violin plots by groups.",
    resolution, ".pdf"
)),
height = 60, width = 7)
VlnPlot(pbmc.combined, features = feature2plot, ncol = 2)
dev.off()

## Dot Plots
pdf(file.path(
    out,
    paste0("Fig 17. DopPlot showing feature expression.",
           resolution, ".pdf")
),
width = 10, height = 4)
DotPlot(pbmc.combined, features = feature2plot,
        dot.scale = 3) + theme(
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 8
            ),
            axis.text.y = element_text(size = 8),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)
        ) +
    scale_x_discrete(breaks = feature2plot,
                     labels = gsub("-ENS.+", "", feature2plot, perl = TRUE)) +
    scale_y_discrete(limits = c(
        "2",
        "3",
        "7",
        "13",
        "4",
        "0",
        "1",
        "12",
        "10",
        "11",
        "8",
        "6",
        "9",
        "5"
    ))
dev.off()


Idents(pbmc.combined) <- "SingleR.labels"
feature2plot_0 <-
    c(  "CD33", "CD86", "FCGR1A",
        "S100A8",
        "S100A9",
        "S100A12",
        "FCGR1B",
        "LYZ",
        "FCGR2A",
        "CD14",
        "TLR1",
        "TLR2",
        "TLR4",
        "TLR6",
        "TLR7",
        "TLR8",
        "CCR2",
        "CD1C",
        "CST3",
        "CD68",
        "PTCRA",
        "MS4A7",
        "FCGR3A",
        "CD4",
        "ITGAM",
        "LGALS3",
        "FCER1G",
        "S100A4",
        "CD44",
        "SELL",
        "LTB",
        "CCR7",
        "CD27",
        "CD3D",
        "CD3E",
        "CD3G",
        "IL7R",
        "CD28",
        "CD8A",
        "CD8B",
        "FOXP3",
        "IL2RA",
        "KLRB1",
        "GZMK",
        "ALOX5AP",
        "GNLY",
        "NKG7",
        "TRDC",
        "GZMH",
        "GZMB",
        "KLRC2",
        "KLRC1",
        "CD79A",
        "TCL1A",
        "CD19",
        "TLR10",
        "IGHA1",
        "IGHG1",
        "TSPAN13",
        "HLA-DRA",
        "HLA-DOA",
        "HLA-DQA1",
        "MS4A1",
        "FCER1A",
        "IL1B",
        "CD38",
        "TYMS",
        "PF4",
        "PPBP"
    )

feature2plot <- rownames(pbmc.combined)[gsub("-ENS.+", "",
                                             rownames(pbmc.combined)) %in%  feature2plot_0]
names(feature2plot) <-
    gsub("-ENS.+", "", feature2plot, perl = TRUE)
feature2plot <- feature2plot[feature2plot_0]

feature2plot <- unname(feature2plot)
pdf(file.path(
    out,
    paste0("Fig 18. DopPlot showing feature expression.singleR.",
           resolution, ".pdf")
),
width = 10, height = 4)
DotPlot(pbmc.combined, features = feature2plot,
        #idents  = "SingleR.labels",
        dot.scale = 3) + theme(
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 8
            ),
            axis.text.y = element_text(size = 8),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)
        ) +
    scale_x_discrete(breaks = feature2plot,
                     labels = gsub("-ENS.+", "", feature2plot, perl = TRUE)) +
    scale_y_discrete(limits = c(
        "cMN", "intMN", "ncMN","mDC", "pDC",    
         "Naive CD4+", "h+fhT", "Treg", "teCD4+",
         "Naive CD8+","cmCD8+", "te+emCD8+", "MAIT","GDT", "NK",
         "Naive B", "mB", "Pb", "Prog" 
    ))
dev.off()


## cell proportion per condition
metadata <- pbmc.combined@meta.data
cell_group <- split(metadata$seurat_clusters, metadata$orig.ident)

prop <- lapply(cell_group, function(.x) {
    freq <- table(.x)
    prop <- round(freq / sum(freq) * 100, digits = 2)
})

prop <- as.data.frame(do.call(cbind, prop))
prop$Cluster <- factor(rownames(prop),
                       levels = sort(as.numeric(as.character(rownames(
                           prop
                       )))))
cell_num <- table(metadata$seurat_clusters)
prop <- cbind(prop, cell_num)

prop <- prop[, c(6, 1:4)]
colnames(prop)[1] <- "Cells"

write.table(
    prop,
    file = file.path(out, "Cell proportion.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

save.image(file.path(out, "05162022.RData"))
