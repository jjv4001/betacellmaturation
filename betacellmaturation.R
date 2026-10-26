library(Seurat)

#Create Seurat Objects
AZD<-Read10X(data.dir = "/home/fs01/jjv4001/betacellchemperturbscRNA/outs/per_sample_outs/AZD/sample_filtered_feature_bc_matrix")
AZD<- CreateSeuratObject(counts = AZD, project = "AZD", min.cells = 3, min.features = 200)
AZD$treatment<-rep("AZD")
Mycophemofetil<-Read10X(data.dir = "/home/fs01/jjv4001/betacellchemperturbscRNA/outs/per_sample_outs/Mycophemofetil/sample_filtered_feature_bc_matrix")
Mycophemofetil<- CreateSeuratObject(counts = Mycophemofetil, project = "Mycophemofetil", min.cells = 3, min.features = 200)
Mycophemofetil$treatment<-rep("Mycophemofetil")
Mycophenoacid<-Read10X(data.dir = "/home/fs01/jjv4001/betacellchemperturbscRNA/outs/per_sample_outs/Mycophenoacid/sample_filtered_feature_bc_matrix")
Mycophenoacid<- CreateSeuratObject(counts = Mycophenoacid, project = "Mycophenoacid", min.cells = 3, min.features = 200)
Mycophenoacid$treatment<-rep("Mycophenoacid")
WT<-Read10X(data.dir = "/home/fs01/jjv4001/betacellchemperturbscRNA/outs/per_sample_outs/WT/sample_filtered_feature_bc_matrix")
WT<- CreateSeuratObject(counts = WT, project = "WT", min.cells = 3, min.features = 200)
WT$treatment<-rep("WT")

#merge objexts
merged<-merge(WT, y=c(AZD, Mycophemofetil, Mycophenoacid))
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
merged<-JoinLayers(merged)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
ElbowPlot(merged)
merged <- FindNeighbors(merged, dims = 1:20)
merged <- FindClusters(merged, resolution = 6.0)
merged <- RunUMAP(merged, dims = 1:20)

##remove doublets 
# get raw counts
counts <- GetAssayData(merged, assay = "RNA", layer = "counts")

# library size per cell
lib.size <- Matrix::colSums(counts)

# CPM
cpm <- t(t(counts) / lib.size) * 1e6

# log2CPM
log2cpm <- log2(cpm + 1)

marker_genes <- c("INS", "GCG", "PPY", "SST", "PRSS1")
marker_genes_present <- marker_genes[marker_genes %in% rownames(log2cpm)]

marker_expr <- log2cpm[marker_genes_present, , drop = FALSE]

# count how many marker genes are > 13 in each cell
n_high_markers <- Matrix::colSums(marker_expr > 13)

# potential doublets = more than one marker above threshold
marker_doublet <- n_high_markers > 1
table(marker_doublet)

merged$marker_based_doublet <- marker_doublet[colnames(merged)]

#remove doublets in scRNA data
merged1<-subset(merged,marker_based_doublet=="FALSE")
merged1 <- NormalizeData(merged1)
merged1 <- FindVariableFeatures(merged1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged1)
merged1 <- ScaleData(merged1, features = all.genes)
merged1 <- RunPCA(merged1, features = VariableFeatures(object = merged1))
ElbowPlot(merged1)
merged1 <- FindNeighbors(merged1, dims = 1:20)
merged1 <- FindClusters(merged1, resolution = 1.0)
merged1 <- RunUMAP(merged1, dims = 1:20)

##Assign celltypes
#Run function to assign cell types
assign_celltypes_by_score <- function(seurat_obj, marker_lists, score_name_prefix = "CellTypeScore", min_score = 0) {
  # 1. Run AddModuleScore on the gene sets
  seurat_obj <- AddModuleScore(seurat_obj, features = marker_lists, name = score_name_prefix)
  
  # 2. Extract score columns
  score_cols <- paste0(score_name_prefix, 1:length(marker_lists))
  scores <- seurat_obj@meta.data[, score_cols]
  
  # 3. Cell type names (in the same order as marker_lists)
  cell_types <- names(marker_lists)
  
  # 4. Assign each cell to cell type with max score
  max_idx <- max.col(scores, ties.method = "first")
  max_scores <- apply(scores, 1, max)
  
  predicted <- ifelse(max_scores >= min_score, cell_types[max_idx], "Unknown")
  
  # 5. Add predicted cell type to metadata
  seurat_obj$predicted_celltype <- predicted
  
  return(seurat_obj)
}

#Define different markers for the respective cell types

markers <- list(
  Beta = c("INS","IAPP","PDX1"),
  Alpha = c("GCG","ARX","IRX2"),
  Delta = c("SST","HHEX","RBP4"),
  Enterochromaffin = c("TPH1","DDC","SLC18A1","CHGA"),
  Ductal = c("KRT19","KRT8","KRT18","SOX9","CFTR"),
  Acinar = c("PRSS1","CPA1","CTRB1"),
  EndocrineProgenitors = c("NEUROG3","NEUROD1","NKX2-2")
)

merged1 <- assign_celltypes_by_score(merged1, markers, min_score = 0.5)
Idents(merged1) <- merged1$predicted_celltype

#Reorder identities
new.cluster.ids <- c("Delta cells","Beta cells","Unknown","Alpha cells","Enterochromaffin cells","Endocrine Progenitors","Ductal cells","Acinar cells")
names(new.cluster.ids) <- levels(beta)
beta <- RenameIdents(beta, new.cluster.ids)

new_levels <- c("Beta cells","Alpha cells","Delta cells","Endocrine Progenitors","Enterochromaffin cells","Ductal cells","Acinar cells","Unknown")
Idents(beta) <- factor(Idents(beta), levels = new_levels)

new_levels <- c("Unknown","Acinar cells","Ductal cells","Enterochromaffin cells","Endocrine Progenitors","Delta cells","Alpha cells","Beta cells")
Idents(beta) <- factor(Idents(beta), levels = new_levels)

features= c("INS","IAPP","PDX1",
            "GCG","ARX","IRX2",
            "SST","HHEX","RBP4",
            "NEUROG3","NEUROD1","NKX2-2",
            "TPH1","DDC","SLC18A1","CHGA",
            "KRT19","KRT8","KRT18","SOX9","CFTR",
            "PRSS1","CPA1","CTRB1")

DotPlot(beta, features=features, dot.scale=5)+RotatedAxis()+scale_color_viridis_c(option = "A", direction=-1) 

#save assigned rds file
saveRDS(beta,"/athena/chenlab/scratch/jjv4001/betacellassigned.rds")

#subset only beta cells, normalize, scale
betaonly<-subset(beta, idents="Beta cells")
betaonly <- NormalizeData(betaonly)
betaonly <- FindVariableFeatures(betaonly, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(betaonly)
betaonly <- ScaleData(betaonly, features = all.genes)
betaonly <- RunPCA(betaonly, features = VariableFeatures(object = betaonly))
ElbowPlot(betaonly)
betaonly <- FindNeighbors(betaonly, dims = 1:20)
betaonly <- FindClusters(betaonly, resolution = 0.1)
betaonly <- RunUMAP(betaonly, dims = 1:20)
Idents(betaonly)<-betaonly$treatment

#obtain differentially expressed gene lists
mycophenoacid.markers <- FindMarkers(betaonly, ident.1 = "Mycophenoacid", ident.2 ="WT")
mycophemofetil.markers <- FindMarkers(betaonly, ident.1 = "Mycophemofetil", ident.2 ="WT")
AZD.markers <- FindMarkers(betaonly, ident.1 = "AZD", ident.2 ="WT")

#obtain differentially expressed gene lists with statistically significant values
mycophenoacid.markers<-subset(mycophenoacid.markers, p_val_adj<0.05)
mycophemofetil.markers<-subset(mycophemofetil.markers, p_val_adj<0.05)
AZD.markers<-subset(AZD.markers, p_val_adj<0.05)
saveRDS(betaonly, "/athena/chenlab/scratch/jjv4001/betaonlychemperturb.rds")
write.csv(mycophenoacid.markers,"/athena/chenlab/scratch/jjv4001/mycophenoacid.markers.csv")
write.csv(mycophemofetil.markers,"/athena/chenlab/scratch/jjv4001/mycophemofetil.markers.csv")
write.csv(AZD.markers,"/athena/chenlab/scratch/jjv4001/AZDmarkers.csv")

#load previous studies of S6 sc-derived beta cells and human islets
otherbeta<-readRDS("/athena/chenlab/scratch/jjv4001/S6cellmerged.rds")
otherbeta$study<-rep("Otherstudies")
humanislets<-readRDS("/athena/chenlab/scratch/jjv4001/ctrl_ml_clustered_final.rds")
humanislets$study<-rep("humanislets")
humanislets$Identity<-humanislets$Cell_Type
beta$study<-rep("thisstudy")

#merge previous studies of S6 sc-derived beta cells and human islets with current study data
allmerged<-merge(beta, y=c(otherbeta, humanislets))
allmerged<-JoinLayers(allmerged)
allmerged[["RNA"]] <- split(allmerged[["RNA"]], f = allmerged$study)
allmerged <- NormalizeData(allmerged)
allmerged <- FindVariableFeatures(allmerged)
allmerged <- ScaleData(allmerged)
allmerged <- RunPCA(allmerged)
allmerged <- FindNeighbors(allmerged , dims = 1:30, reduction = "pca")
allmerged <- FindClusters(allmerged , resolution = 2, cluster.name = "unintegrated_clusters")
allmerged  <- RunUMAP(allmerged , dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
#integrate previous studies of S6 sc-derived beta cells and human islets with current study data
allmerged[["RNA"]] <- split(allmerged[["RNA"]], f = allmerged$study)
allmerged <- IntegrateLayers(object = allmerged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
allmerged[["RNA"]] <- JoinLayers(allmerged[["RNA"]])

#run UMAP and rename data
allmerged <- FindNeighbors(allmerged, reduction = "integrated.cca", dims = 1:30)
allmerged<- FindClusters(allmerged, resolution = 1)
allmerged <- RunUMAP(allmerged, dims = 1:30, reduction = "integrated.cca", reduction.name="umap.integrated")
Idents(allmerged)<-allmerged$Identity
allmerged$Identity[allmerged$Identity == "Acinar"] <- "Acinar cells"
allmerged$Identity[allmerged$Identity == "Alpha"] <- "Alpha cells"
allmerged$Identity[allmerged$Identity == "Beta"] <- "Beta cells"
allmerged$Identity[allmerged$Identity == "Delta"] <- "Delta cells"
allmerged$Identity[allmerged$Identity == "Ductal"] <- "Ductal cells"
allmerged$Identity[allmerged$Identity == "Endocrine progenitors"] <- "Endocrine Progenitors"
allmerged$Identity[allmerged$Identity == "Polyhormonal"] <- "Polyhormonal cells"
saveRDS(allmerged, "/athena/chenlab/scratch/jjv4001/allmergedbetacellchempeturb.rds")

#subset only beta cells from all studies 
betacellmerged<-subset(allmerged, Identity=="Beta cells")
betacellmerged <- NormalizeData(betacellmerged)
betacellmerged <- FindVariableFeatures(betacellmerged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(betacellmerged)
betacellmerged <- ScaleData(betacellmerged, features = all.genes)
betacellmerged <- RunPCA(betacellmerged, features = VariableFeatures(object = betacellmerged))
ElbowPlot(betacellmerged)
betacellmerged <- FindNeighbors(betacellmerged, dims = 1:20)
betacellmerged <- FindClusters(betacellmerged, resolution = 0.1)
betacellmerged <- RunUMAP(betacellmerged, dims = 1:20)
DimPlot(betacellmerged, reduction="umap.integrated", group.by="study")
saveRDS(betacellmerged,"/athena/chenlab/scratch/jjv4001/betacellallstudiesmerged.rds")

#plot bulk PCA plot profiles
library(dplyr)
library(Seurat)

# Create grouping variable
betacell$group <- paste0(betacell$study,
                              betacell$treatment
)

pb <- AggregateExpression(
  betacell,
  assays = "RNA",
  group.by = "group",
  slot = "counts",
  return.seurat = FALSE
)

pb_mat <- pb$RNA

# log transform
pb_log <- log1p(pb_mat)

# remove genes with zero variance
pb_log <- pb_log[apply(pb_log, 1, var) > 0, ]

# run PCA
pca <- prcomp(t(pb_log), scale. = TRUE)

library(ggplot2)
library(ggrepel)

pca_df <- data.frame(
  group = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

ggplot(pca_df, aes(PC1, PC2, color = group, label = group)) +
  geom_point(size = 4) +
  geom_text_repel() +
  theme_classic() +
  labs(title = "Pseudobulk PCA by group")

  
  library(matrixStats)

# convert to matrix
pb_log <- as.matrix(pb_log)

# compute gene variances
vars <- rowVars(pb_log)

# select top variable genes
top_genes <- names(sort(vars, decreasing = TRUE))[1:100]

# PCA
pca <- prcomp(t(pb_log[top_genes, ]), scale. = TRUE)

library(ggplot2)
library(ggrepel)

pca_df <- data.frame(
  group = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

ggplot(pca_df, aes(PC1, PC2, color = group, label = group)) +
  geom_point(size = 4) +
  geom_text_repel() +
  theme_classic() +
  labs(title = "Pseudobulk PCA by group")

##plot module scores based on beta cell maturity genes-10 genes
beta_maturity_genes <- c(
  "MAFA",
  "UCN3",
  "IAPP",
  "G6PC2",
  "SLC30A8",
  "PCSK1",
  "ABCC8",
  "KCNJ11",
  "GCK",
  "GLP1R"
)

# keep only genes present in object
beta_maturity_genes_present <- intersect(beta_maturity_genes, rownames(beta))

beta <- AddModuleScore(
  object = beta,
  features = list(beta_maturity_genes_present),
  name = "BetaMaturity",
  assay = DefaultAssay(beta)
)

##plot module scores based on beta cell maturity genes-11 genes
library(dplyr)
library(ggplot2)

avg_scores <- beta@meta.data %>%
  group_by(treatmentstudy) %>%
  summarise(avg_score = mean(BetaMaturity1, na.rm = TRUE))

ggplot(avg_scores, aes(x = reorder(treatmentstudy, avg_score), y = avg_score, fill = avg_score)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  labs(x = "Cell type", y = "Average beta maturity score")


beta_maturity_genes <- c(
  "MAFA",
  "UCN3",
  "IAPP",
  "G6PC2",
  "SLC30A8",
  "PCSK1",
  "ABCC8",
  "KCNJ11",
  "GCK",
  "GLP1R","PDX1"
)

# keep only genes present in object
beta_maturity_genes_present <- intersect(beta_maturity_genes, rownames(beta))

beta <- AddModuleScore(
  object = beta,
  features = list(beta_maturity_genes_present),
  name = "BetaMaturity",
  assay = DefaultAssay(beta)
)

library(dplyr)
library(ggplot2)

avg_scores <- beta@meta.data %>%
  group_by(treatmentstudy) %>%
  summarise(avg_score = mean(BetaMaturity1, na.rm = TRUE))

ggplot(avg_scores, aes(x = reorder(treatmentstudy, avg_score), y = avg_score, fill = avg_score)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  labs(x = "Cell type", y = "Average beta maturity score")



