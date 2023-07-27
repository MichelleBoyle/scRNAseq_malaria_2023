# scRNAseq_malaria_2023 - https://www.biorxiv.org/content/10.1101/2022.11.16.516822v1 
3' scRNAseq of PBMCs from malaria

##### Session/Software info ####
# Typical install time on a "normal" desktop computer ~ 30 minutes
# Typical run time on a "normal" desktop computer ~ 45 minutes. Exception:"PBMC Data Integration and Clustering"
##### Set working directory ######
# Set working directory as parent folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#### Set Lib Path ####
myPaths <- .libPaths()
myPaths <- c(myPaths, "./renviron")
.libPaths(myPaths)  # add new path
# Downloaded packages are in zip folder, Change R library path to access them or alternatively download packages using code below (<20 minutes)
##### Install packages#####
install.packages("Seurat", "./renviron")
install.packages("SeuratObject", "~/renviron")
install.packages("dplyr", "~/renviron")
install.packages("stringr", "~/renviron")
install.packages("ggplot2", "~/renviron")
install.packages("ggpubr", "~/renviron")
install.packages("tibble", "~/renviron")
install.packages("UpSetR", "~/renviron")
install.packages("ggpubr", "~/renviron")
install.packages("BiocManager", "~/renviron")
BiocManager::install("MAST", lib = "~/renviron", force = TRUE)
##### Load packages #####
library(Seurat)
library(SeuratObject)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(UpSetR)
library(BiocManager)
library(MAST)
# Total runtime expected for all code after package loaded ~45 minutes
##### Import single-cell RNA-seq count data #####
# Ensure there is a folder within working directory "Cell_ranger_outputs" that contains the below folders with these files in each folder "barcodes.tsv.gz", "features.tsv.gz","matrix.tsv.gz". Optionally data can be downloaded from this https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217930
# We use the for loop to execute two commands for each sample - (1) read in the count data (Read10X()) and (2) create the Seurat objects from the read in data (CreateSeuratObject()):
for (file in c("A_filtered_feature_bc_matrix" ,"B_filtered_feature_bc_matrix","C_filtered_feature_bc_matrix","D_filtered_feature_bc_matrix","E_filtered_feature_bc_matrix","F_filtered_feature_bc_matrix","G_filtered_feature_bc_matrix","H_filtered_feature_bc_matrix","I_filtered_feature_bc_matrix","J_filtered_feature_bc_matrix","K_filtered_feature_bc_matrix","L_filtered_feature_bc_matrix","M_filtered_feature_bc_matrix","N_filtered_feature_bc_matrix","O_filtered_feature_bc_matrix","P_filtered_feature_bc_matrix","Q_filtered_feature_bc_matrix","R_filtered_feature_bc_matrix","S_filtered_feature_bc_matrix","T_filtered_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("Cell_ranger_outputs/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3,
                                   min.features = 200,
                                   project = file)
  assign(file, seurat_obj)
}
# Take a quick look at the metadata to see how it looks
head(A_filtered_feature_bc_matrix@meta.data)
head(B_filtered_feature_bc_matrix@meta.data)
head(C_filtered_feature_bc_matrix@meta.data)
head(D_filtered_feature_bc_matrix@meta.data)
head(E_filtered_feature_bc_matrix@meta.data)
head(F_filtered_feature_bc_matrix@meta.data)
head(G_filtered_feature_bc_matrix@meta.data)
head(H_filtered_feature_bc_matrix@meta.data)
head(I_filtered_feature_bc_matrix@meta.data)
head(J_filtered_feature_bc_matrix@meta.data)
head(K_filtered_feature_bc_matrix@meta.data)
head(L_filtered_feature_bc_matrix@meta.data)
head(M_filtered_feature_bc_matrix@meta.data)
head(N_filtered_feature_bc_matrix@meta.data)
head(O_filtered_feature_bc_matrix@meta.data)
head(P_filtered_feature_bc_matrix@meta.data)
head(Q_filtered_feature_bc_matrix@meta.data)
head(R_filtered_feature_bc_matrix@meta.data)
head(S_filtered_feature_bc_matrix@meta.data)
head(T_filtered_feature_bc_matrix@meta.data)
# We will merge these objects together into a single Seurat object. This will make it easier to run the QC steps for all sample groups together and enable us to easily compare the data quality for all the samples.
merged_seurat <- merge(x = A_filtered_feature_bc_matrix,
                       y = c(B_filtered_feature_bc_matrix,C_filtered_feature_bc_matrix,D_filtered_feature_bc_matrix,E_filtered_feature_bc_matrix,F_filtered_feature_bc_matrix,G_filtered_feature_bc_matrix,H_filtered_feature_bc_matrix,I_filtered_feature_bc_matrix,J_filtered_feature_bc_matrix,K_filtered_feature_bc_matrix,L_filtered_feature_bc_matrix,M_filtered_feature_bc_matrix,N_filtered_feature_bc_matrix,O_filtered_feature_bc_matrix,P_filtered_feature_bc_matrix,Q_filtered_feature_bc_matrix,R_filtered_feature_bc_matrix,S_filtered_feature_bc_matrix,T_filtered_feature_bc_matrix),
                       add.cell.id = c("adult1day0", "adult1day7", "adult1day28","adultc1","child1day0","child1day7","child1day28","adultc2","adult2day0", "adult2day28", "adult3day0", "adult3day28", "child2day0", "child2day28", "child3day0", "child3day28", "adult2day7", "adult3day7", "child2day7", "child3day7"))
#Look at the first and last few rows of the counts matrix within the merged_seurat object to check if it has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
# Explore merged metadata
View(merged_seurat@meta.data)
# Aside from the number of UMIs per cell (nCount_RNA) and number of genes detected per cell (nFeature_RNA), We need to calculate some additional metrics (1) the number of genes detected per UMI: this metric gives us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data) (2) mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes
# The number of genes per UMI for each cell is quite easy to calculate, and we will log10 transform the result for better comparison between samples.
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
# Create the metadata dataframe by extracting the meta.data slot from the Seurat object
metadata <- merged_seurat@meta.data
head(metadata)
metadata$cells <- rownames(metadata)
# Add cell IDs to metadata
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample columns: Get sample names for each of the cells based on the cell prefix
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^adult1day0_"))] <- "adult1day0"
metadata$sample[which(str_detect(metadata$cells, "^adult1day7_"))] <- "adult1day7"
metadata$sample[which(str_detect(metadata$cells, "^adult1day28_"))] <- "adult1day28"
metadata$sample[which(str_detect(metadata$cells, "^adultc1_"))] <- "adultc1"
metadata$sample[which(str_detect(metadata$cells, "^child1day0_"))] <- "child1day0"
metadata$sample[which(str_detect(metadata$cells, "^child1day7_"))] <- "child1day7"
metadata$sample[which(str_detect(metadata$cells, "^child1day28_"))] <- "child1day28"
metadata$sample[which(str_detect(metadata$cells, "^adultc2_"))] <- "adultc2"
metadata$sample[which(str_detect(metadata$cells, "^adult2day0_"))] <- "adult2day0"
metadata$sample[which(str_detect(metadata$cells, "^adult2day28_"))] <- "adult2day28"
metadata$sample[which(str_detect(metadata$cells, "^adult3day0_"))] <- "adult3day0"
metadata$sample[which(str_detect(metadata$cells, "^adult3day28_"))] <- "adult3day28"
metadata$sample[which(str_detect(metadata$cells, "^child2day0_"))] <- "child2day0"
metadata$sample[which(str_detect(metadata$cells, "^child2day28_"))] <- "child2day28"
metadata$sample[which(str_detect(metadata$cells, "^child3day0_"))] <- "child3day0"
metadata$sample[which(str_detect(metadata$cells, "^child3day28_"))] <- "child3day28"
metadata$sample[which(str_detect(metadata$cells, "^adult2day7_"))] <- "adult2day7"
metadata$sample[which(str_detect(metadata$cells, "^adult3day7_"))] <- "adult3day7"
metadata$sample[which(str_detect(metadata$cells, "^child2day7_"))] <- "child2day7"
metadata$sample[which(str_detect(metadata$cells, "^child3day7_"))] <- "child3day7"
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata
# Save R object as RDS to load at any time
saveRDS(merged_seurat, "merged_seurat.RDS")
##### Assess quality metrics #####
# Visualise the number of cells per sample
metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells per sample")
ggsave("raw_number of cells per sample.pdf")
# Visualize the number UMIs/transcripts per cell
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMI counts(transcripts per cell)")
ggsave("raw_number of UMIs or transcripts per sample.pdf")
# Visualize the distribution of genes detected per cell via histogram
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("nGenes detected per cell")
ggsave("raw_distribution of genes detected per cell-histogram.pdf")
# Visualize the distribution of genes detected per cell via boxplot
metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("raw_distribution of genes detected per cell-boxplot.pdf")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>%
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) +
  ggtitle("Correlation between genes and UMIs and presence of mitoRNA")
ggsave("raw_correlation between genes and UMIs and presence of mitoRNA.pdf")
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>%
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.5)+
  ggtitle("Mitochondrial counts ratio")
ggsave("raw_mitochondrial counts ratio.pdf")
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  ggtitle("Complexity (genes detected per UMI)")
ggsave("raw_complexity, genes detected per UMI.pdf")
# Visualise QC by ViolinPlot
VlnPlot(merged_seurat, features = c("nUMI", "nGene","log10GenesPerUMI", "mitoRatio"), group.by = "sample", pt.size = 0, ncol = 4) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("raw_Vlnplot_sample.pdf", h=4.01, w=15.79)
# Create meta.data for timpeoints and controls
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult1day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult2day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult3day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child1day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child2day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child3day0"] <- "Day0"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult1day7"] <- "Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult2day7"] <-"Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult3day7"] <-"Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child1day7"] <- "Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child2day7"] <-"Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child3day7"] <- "Day7"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult1day28"] <- "Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult2day28"] <-"Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adult3day28"] <-"Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child1day28"] <- "Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child2day28"] <-"Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "child3day28"] <- "Day28"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adultc1"] <- "Control"
merged_seurat@meta.data$day[merged_seurat@meta.data$sample== "adultc2"] <- "Control"
# Define an order of cluster identities
day_levels <- c("Day0", "Day7", "Day28", "Control")
# Re-level object@ident
merged_seurat@meta.data$day <- factor(x = merged_seurat@meta.data$day, levels = day_levels)
# QC violin plot according to Day
VlnPlot(merged_seurat, features = c("nUMI", "nGene","log10GenesPerUMI", "mitoRatio"), group.by = "day", pt.size = 0, ncol = 4) +theme(axis.text.x = element_text(angle = 90))
ggsave("raw_Vlnplot_day.pdf", h=4.01, w=5.41)
##### Cell Filtering#####
# Filter out low quality reads using manually/supervised selected thresholds to clean data and create consistency- these will change with experiment
filtered_seurat <- subset(x = merged_seurat,
                          subset= (nUMI >= 500) &
                            (nGene >= 250) &
                            (log10GenesPerUMI > 0.80) &
                            (mitoRatio < 0.20))
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# Create a clean (post-filtered) metadata dataframe by extracting the meta.data slot from the flitered Seurat object
metadata_clean <- filtered_seurat@meta.data
##### Assess quality metrics post-filter #####
# Visualise the number of cells per sample
metadata_clean %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells per sample")
ggsave("filtered_number of cells per sample.pdf")
# Visualize the number UMIs/transcripts per cell
metadata_clean %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMI counts(transcripts per cell)")
ggsave("filtered_number of UMIs or transcripts per sample.pdf")
# Visualize the distribution of genes detected per cell via histogram
metadata_clean %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)+
  ggtitle("nGenes detected per cell")
ggsave("filtered_distribution of genes detected per cell-histogramt.pdf")
# Visualize the distribution of genes detected per cell via boxplot
metadata_clean %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
ggsave("filtered_distribution of genes detected per cell-boxplot.pdf")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>%
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample) +
  ggtitle("Correlation between genes and UMIs and presence of mitoRNA")
ggsave("filtered_correlation between genes and UMIs and presence of mitoRNA.pdf")
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_clean %>%
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.5)+
  ggtitle("Mitochondrial counts ratio")
ggsave("filtered_mitochondrial counts ratio.pdf")
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  ggtitle("Complexity (genes detected per UMI)")
ggsave("filtered_complexity, genes detected per UMI.pdf")
# QC violin plot according to Sample
VlnPlot(filtered_seurat, features = c("nUMI", "nGene","log10GenesPerUMI", "mitoRatio"), group.by = "sample", pt.size = 0, ncol = 4) +theme(axis.text.x = element_text(angle = 90))
ggsave("filtered_Vlnplot_sample.pdf", h=4.01, w=15.79)
# QC violin plot according to Day
VlnPlot(filtered_seurat, features = c("nUMI", "nGene","log10GenesPerUMI", "mitoRatio"), group.by = "day", pt.size = 0, ncol = 4) +theme(axis.text.x = element_text(angle = 90))
ggsave("filtered_Vlnplot_day.pdf", h=4.01, w=5.41)
# Save R object as RDS to load at any time
saveRDS(filtered_seurat, "filtered_seurat.RDS")



# NOT SURE IF THIS IS NECESSARY MAYBE JUST A TEST
##### Cell Cycle scoring ####



#probably not needed
# Define G2M and S phase cell cycle based on cannonical markers, (?????reference??????)


g2m_genes <- c("NCAPD2", "ANLN", "TACC3", "HMMR", "GTSE1", "NDC80", "AURKA", "TPX2", "BIRC5", "G2E3", "CBX5", "RANGAP1", "CTCF", "CDCA3", "TTK", "SMC4", "ECT2", "CENPA", "CDC20", "NEK2", "CENPF", "TMPO", "HJURP", "CKS2", "DLGAP5", "PIMREG", "TOP2A", "PSRC1", "CDCA8", "CKAP2", "NUSAP1", "KIF23", "KIF11", "KIF20B", "CENPE", "GAS2L3", "KIF2C", "NUF2", "ANP32E", "LBR", "MKI67", "CCNB2", "CDC25C", "HMGB2", "CKAP2L", "BUB1", "CDK1", "CKS1B", "UBE2C", "CKAP5", "AURKB", "CDCA2", "TUBB4B", "JPT1")
s_genes <- c("UBR7", "RFC2", "RAD51", "MCM2", "TIPIN", "MCM6", "UNG", "POLD3", "WDR76", "CLSPN", "CDC45", "CDC6", "MSH2", "MCM5", "POLA1", "MCM4", "RAD51AP1", "GMNN", "RPA2", "CASP8AP2", "HELLS", "E2F8", "GINS2", "PCNA", "NASP", "BRIP1", "DSCC1", "DTL", "CDCA7", "CENPU", "ATAD2", "CHAF1B", "USP1", "SLBP", "RRM1", "FEN1", "RRM2", "EXO1", "CCNE2", "TYMS", "BLM", "PRIM1", "UHRF1")
# Score cells for cell cycle 
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE)
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)








##### Remove residual RBC #####
#Filter out residual haemoglobin-associated genes
counts <- GetAssayData(filtered_seurat, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('HBB','HBA1','HBA2','HBG2','HBG1','HBD'))),]
filtered_seurat <- subset(filtered_seurat, features = rownames(counts))
options(future.globals.maxSize = 4000 * 1024^4)
# Save R object as RDS to load at any time
saveRDS(filtered_seurat, "filtered_seurat.RDS")
##### PBMC Data Integration and Clustering #####
### WARNING: This task requires significant processing power: ~104gB RAM, ~6CPUs, 5hr run-time
### If reviewers do not have access to a device capable, we recommend to skip to "PBMC cluster visualisation" and import object - "seurat_integrated.rds" in submission folder
# Split seurat object by sample to perform cell cycle scoring and normalization on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")
split_seurat <- split_seurat[c("adult1day0", "adult1day7", "adult1day28","adultc1","child1day7","child1day28","adultc2","adult2day0", "adult2day28", "adult3day0", "adult3day28", "child2day0", "child2day28", "child3day0", "child3day28", "adult2day7", "adult3day7", "child2day7", "child3day7")]
# Defined G2M and S phase cell cycle based on cannonical markers, (?????reference??????)
g2m_genes <- c("NCAPD2", "ANLN", "TACC3", "HMMR", "GTSE1", "NDC80", "AURKA", "TPX2", "BIRC5", "G2E3", "CBX5", "RANGAP1", "CTCF", "CDCA3", "TTK", "SMC4", "ECT2", "CENPA", "CDC20", "NEK2", "CENPF", "TMPO", "HJURP", "CKS2", "DLGAP5", "PIMREG", "TOP2A", "PSRC1", "CDCA8", "CKAP2", "NUSAP1", "KIF23", "KIF11", "KIF20B", "CENPE", "GAS2L3", "KIF2C", "NUF2", "ANP32E", "LBR", "MKI67", "CCNB2", "CDC25C", "HMGB2", "CKAP2L", "BUB1", "CDK1", "CKS1B", "UBE2C", "CKAP5", "AURKB", "CDCA2", "TUBB4B", "JPT1")
s_genes <- c("UBR7", "RFC2", "RAD51", "MCM2", "TIPIN", "MCM6", "UNG", "POLD3", "WDR76", "CLSPN", "CDC45", "CDC6", "MSH2", "MCM5", "POLA1", "MCM4", "RAD51AP1", "GMNN", "RPA2", "CASP8AP2", "HELLS", "E2F8", "GINS2", "PCNA", "NASP", "BRIP1", "DSCC1", "DTL", "CDCA7", "CENPU", "ATAD2", "CHAF1B", "USP1", "SLBP", "RRM1", "FEN1", "RRM2", "EXO1", "CCNE2", "TYMS", "BLM", "PRIM1", "UHRF1")
for (i in 1:length(x = split_seurat)) {
  split_seurat[[i]] <- NormalizeData(object = split_seurat[[i]], verbose = FALSE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- FindVariableFeatures(object = split_seurat[[i]], selection.method = "vst",
                                            nfeatures = 2000, verbose = FALSE)
}
# Find Integration Anchors and integrate data 
integ.anchors <- FindIntegrationAnchors(object.list =  split_seurat, dims = 1:20)
seurat_integrated <- IntegrateData(anchorset = integ.anchors, dims = 1:20)
DefaultAssay(seurat_integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seurat_integrated <- ScaleData(seurat_integrated, vars.to.regress = c("mitoRatio"))
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30)
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30, reduction = "pca")
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:30)
seurat_integrated <- FindClusters(object = seurat_integrated, resolution=0.6)
DimPlot(seurat_integrated, reduction = "umap", label = TRUE, label.size = 6)
saveRDS(seurat_integrated, "seurat_integrated.RDS")
##### PBMC Cluster visualisation ####
seurat_integrated <-readRDS("seurat_integrated.rds")
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult1day0"] <- "Adult1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult2day0"] <- "Adult2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult3day0"] <- "Adult3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child2day0"] <- "Child2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child3day0"] <- "Child3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult1day7"] <- "Adult1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult2day7"] <- "Adult2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult3day7"] <- "donor7"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child1day7"] <- "Child1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child2day7"] <- "Child2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child3day7"] <- "Child3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult1day28"] <- "Adult1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult2day28"] <- "Adult2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adult3day28"] <- "Adult3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child1day28"] <- "Child1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child2day28"] <- "Child2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="child3day28"] <- "Child3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adultc1"] <- "Control1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adultc2"] <- "Control2"
#Set identity at integrated_snn_res.0.6 to compare clusters at chosen resolution
Idents(seurat_integrated) <- "integrated_snn_res.0.6"
# Visualise canonical PBMC gene expression in cluster 
DotPlot(seurat_integrated, features=c("CD14","LYZ", "S100A8", "S100A9", "IL1A", "IL1B", "TNF", "FCER1G", "CTSS", "FCGR3A", "CX3CR1",
                                      "CDKN1C", "IFITM3", "FCER1A", "CD1C", "HLA-DRA", "HLA-DRB1", "CD74", "CLEC4C", "LILRA4", "PLD4",
                                      "SERPINF1", "IL3RA", "CCR7", "GIMAP7", "IL7R", "MAL", "LTB", "NELL2", "PASK", "CD8A", "CD8B",
                                      "KLRB1", "KLRG1", "TRGC1", "TRGC2", "TRDC", "NCAM1", "GNLY", "NKG7", "GZMA", "GZMB", "KLRC1",
                                      "KLRD1", "CD19", "CD79A", "BANK1", "IGHM","IGHD", "MZB1", "CD38", "JCHAIN", "IGHA1", "IGHG1",
                                      "IGHG3", "MKI67", "PCNA", "STMN1", "PCLAF", "CD34", "SPINK2", "SMIM24", "AVP", "PPBP", "ITGB3",
                                      "PF4", "NRGN"),cols=c("white", "red3")) +coord_flip()+ theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
ggsave("Clustering_compare-Dotplot.pdf", h =16, w=14)
# Based on canonical genes rename cluster identitites 
seurat_integrated <- RenameIdents(object = seurat_integrated,   
                                  "0"= "CD4 T cells",
                                  "1"="CD4 T cells",
                                  "2"= "CD8 T cells",
                                  "3" = "NK cells",
                                  "4" = "CD8 T cells",
                                  "5" = "CD4 T cells",
                                  "6" = "γδ T cells",
                                  "7" = "CD14 Monocytes",
                                  "8" = "CD4 T cells",
                                  "9" = "NK cells",
                                  "10" = "B cells",
                                  "11" = "B cells",
                                  "12" = "CD8 T cells",
                                  "13" = "CD4 T cells",
                                  "14" = "CD4 T cells",
                                  "15" = "Uncharacterised",
                                  "16" = "CD16 Monocytes",
                                  "17" = "NKT cells",
                                  "18" = "Plasma cells",
                                  "19" = "cDCs",
                                  "20" = "pDCs",
                                  "21" = "HSPCs",
                                  "22" = "Platelets")
# Visualise mKi67 expression to find proliferating cells
FeaturePlot(seurat_integrated, reduction = "umap", features = "MKI67", raster = TRUE)
ggsave("MIK67-FeaturePlot.pdf", h =6, w = 6)
# Using cell selector rename distinct MKI67 cluster as 'Proliferating cells'
plot <- DimPlot(seurat_integrated)
cells.proliferating <- CellSelector(plot = plot)
#Click done in GUI, then select cells again below
seurat_integrated <- CellSelector(plot = plot, object = seurat_integrated, ident = "Proliferating cells")
# Check cell selector has worked
DimPlot(seurat_integrated, reduction = "umap", label = FALSE, raster = TRUE)
# Add cluster annotation in active identity to meta.data and order identities for dotplot visualisation
seurat_integrated[["annotate_cluster"]] <-Idents(object = seurat_integrated)
seurat_integrated@meta.data$annotate_cluster <- factor(seurat_integrated@meta.data$annotate_cluster, levels = c("CD14 Monocytes","CD16 Monocytes","cDCs","pDCs","CD4 T cells","CD8 T cells","γδ T cells","NKT cells","NK cells","B cells","Plasma cells","Proliferating cells","HSPCs","Platelets","Uncharacterised"))
#Check annotation
unique(seurat_integrated@meta.data$annotate_cluster)
#Set annotation as active identity
Idents(seurat_integrated) <- "annotate_cluster"
#Check annotated clusters expression of canonical genes
DotPlot(seurat_integrated, features=c("CD14","LYZ", "S100A8", "S100A9", "IL1A", "IL1B", "TNF", "FCER1G", "CTSS", "FCGR3A", "CX3CR1",
                                      "CDKN1C", "IFITM3", "FCER1A", "CD1C", "HLA-DRA", "HLA-DRB1", "CD74", "CLEC4C", "LILRA4", "PLD4",
                                      "SERPINF1", "IL3RA", "CCR7", "GIMAP7", "IL7R", "MAL", "LTB", "NELL2", "PASK", "CD8A", "CD8B",
                                      "KLRB1", "KLRG1", "TRGC1", "TRGC2", "TRDC", "NCAM1", "GNLY", "NKG7", "GZMA", "GZMB", "KLRC1",
                                      "KLRD1", "CD19", "CD79A", "BANK1", "IGHM","IGHD", "MZB1", "CD38", "JCHAIN", "IGHA1", "IGHG1",
                                      "IGHG3", "MKI67", "PCNA", "STMN1", "PCLAF", "CD34", "SPINK2", "SMIM24", "AVP", "PPBP", "ITGB3",
                                      "PF4", "NRGN"),cols=c("white", "red3")) +coord_flip()+ theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
ggsave("AnnotatedClustering_compare-Dotplot.pdf", h =16, w=14)
##### Identify DEGs between days####
# This is an example using CD14 monocytes of how we analysed our clusters and identified DEGs between days
# Create Donor meta.data
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("adult1day0","adult1day7","adult1day28" )] <- "Adult1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("adult2day0","adult2day7","adult2day28" )] <- "Adult2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("adult3day0","adult3day7","adult3day28" )] <- "Adult3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("child1day7","child1day28" )] <- "Child1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("child2day0","child2day7","child2day28" )] <- "Child2"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident %in% c("child3day0","child3day7","child3day28" )] <- "Child3"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adultc1"] <- "Control1"
seurat_integrated@meta.data$donor[seurat_integrated@meta.data$orig.ident=="adultc2"] <- "Control2"
# Create Day meta.data
seurat_integrated@meta.data$Day[seurat_integrated@meta.data$orig.ident %in% c("adult1day0","adult2day0","adult3day0","child2day0","child3day0")] <- "Day0"
seurat_integrated@meta.data$Day[seurat_integrated@meta.data$orig.ident %in% c("adult1day7","adult2day7","adult3day7","child1day7","child2day7","child3day7")] <- "Day7"
seurat_integrated@meta.data$Day[seurat_integrated@meta.data$orig.ident %in% c("adult1day28","adult2day28","adult3day28","child1day28","child2day28","child3day28")] <- "Day28"
seurat_integrated@meta.data$Day[seurat_integrated@meta.data$orig.ident %in% c("adultc1","adultc2")] <- "Control"
# Add meta.data for the group and day
seurat_integrated$celltype_day <- paste(seurat_integrated$annotate_cluster, seurat_integrated$Day,  sep = "_")
##use the celltype.day to find DEGs for each subset and day combination i.e. find DEGS for CD14 Monocytes comparing Day 0 with Day 7, Day 0 with Day 28, Day 7 with Day 28.
Idents(seurat_integrated) <- "celltype_day"
DefaultAssay(seurat_integrated) <- "RNA"
CD14Monos_d0v7 <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day0", ident.2 = "CD14 Monocytes_Day7", test.use = "MAST")
CD14Monos_d0v7['Comparison'] <- "0v7"
write.csv(CD14Monos_d0v7, "CD14Monos_d0v7.csv")
with(CD14Monos_d0v7, sum(p_val_adj < 0.05))
#1556 DEGs with p_val_adj < 0.05
CD14Monos_d0v28 <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day0", ident.2 = "CD14 Monocytes_Day28", test.use = "MAST")
CD14Monos_d0v28['Comparison'] <- "0v28"
write.csv(CD14Monos_d0v28, "CD14Monos_d0v28.csv")
with(CD14Monos_d0v28, sum(p_val_adj < 0.05))
#1677 DEGs with p_val_adj < 0.05
CD14Monos_d7v28 <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day7", ident.2 = "CD14 Monocytes_Day28", test.use = "MAST")
CD14Monos_d7v28['Comparison'] <- "7v28"
write.csv(CD14Monos_d7v28, "CD14Monos_d7v28.csv")
with(CD14Monos_d7v28, sum(p_val_adj < 0.05))
#286 DEGs with p_val_adj < 0.05
CD14Monos_d0vC <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day0", ident.2 = "CD14 Monocytes_Control", test.use = "MAST")
CD14Monos_d0vC['Comparison'] <- "0vC"
write.csv(CD14Monos_d0vC, "CD14Monos_d0vC.csv")
with(CD14Monos_d0vC, sum(p_val_adj < 0.05))
#1698 DEGs with p_val_adj < 0.05
CD14Monos_d7vC <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day7", ident.2 = "CD14 Monocytes_Control", test.use = "MAST")
CD14Monos_d7vC['Comparison'] <- "7vC"
write.csv(CD14Monos_d7vC, "CD14Monos_d7vC.csv")
with(CD14Monos_d7vC, sum(p_val_adj < 0.05))
#506 DEGs with p_val_adj < 0.05
CD14Monos_d28vC <- FindMarkers(seurat_integrated, ident.1 = "CD14 Monocytes_Day28", ident.2 = "CD14 Monocytes_Control", test.use = "MAST")
CD14Monos_d28vC['Comparison'] <- "28vC"
write.csv(CD14Monos_d28vC, "CD14Monos_d28vC.csv")
with(CD14Monos_d28vC, sum(p_val_adj < 0.05))
#764 DEGs with p_val_adj < 0.05
saveRDS(seurat_integrated, "seurat_integrated_annot.rds")
##### Sub-cluster analysis example gamma-delta T cells#### 
# Here is an example of how we re-import our PBMC subclusters as their own seurat object and perform further clustering,classification and DEG analysis. This code will take less than 20 minutes to run
seurat_integrated_annot <-readRDS("seurat_integrated_annot.rds")
Idents(seurat_integrated_annot) <- "annotate_cluster"
annotated_gd <-subset(seurat_integrated_annot, idents = "γδ T cells")
#QC 
pdf("QC_gd.pdf", h = 8, w = 12)
VlnPlot(annotated_gd, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), pt.size = 0.001)
FeatureScatter(annotated_gd, "nFeature_RNA", "nCount_RNA")
dev.off()
##Based on these graphs we'll remove nFeature_RNA > 500 & nFeature_RNA < 3500 & mitoRatio > 0.015
table(annotated_gd@meta.data$Day)
# Control    Day0   Day28    Day7 
# 1151       867    2877     2263 
annotated_gdF <- subset(annotated_gd, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & mitoRatio > 0.015)  
table(annotated_gdF@meta.data$Day)
# Control    Day0   Day28    Day7 
# 1110       819    2794     2202 
##Saving filtered seurat object
#Filtered data QC
pdf("QC_postfilter.pdf", h = 8, w = 12)
VlnPlot(annotated_gdF, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),  pt.size = 0.001)
FeatureScatter(annotated_gdF, "nFeature_RNA", "nCount_RNA")
dev.off()
#Findvariable genes
annotated_gdF <- FindVariableFeatures(annotated_gdF, selection.method = 'vst', nfeatures = 10000)
# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(annotated_gdF), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(annotated_gdF)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#> When using repel, set xnudge and ynudge to 0 for optimal results
pdf("Top20_variable genes.pdf", h = 8, w= 15)
plot1 + plot2
dev.off()
# Scale Data and Run PCA 
DefaultAssay(annotated_gdF) <- "integrated"
annotated_gdF <- ScaleData(annotated_gdF, vars.to.regress = c("mitoRatio"))
annotated_gdF <- RunPCA(annotated_gdF, npcs = 30)
# Examine and visualize PCA results a few different ways
print(annotated_gdF[["pca"]], dims = 1:5, nfeatures = 10)
VizDimLoadings(annotated_gdF, dims = 1:2, reduction = "pca")
VizDimLoadings(annotated_gdF, dims = 3:4, reduction = "pca")
VizDimLoadings(annotated_gdF, dims = 5, reduction = "pca")
DimHeatmap(annotated_gdF, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(annotated_gdF)
# Run UMAP
annotated_gdF <- RunUMAP(annotated_gdF, dims = 1:30, reduction = "pca")
# Determine the K-nearest neighbor graph
annotated_gdF <- FindNeighbors(object = annotated_gdF, dims = 1:30)
annotated_gdF <- FindClusters(object = annotated_gdF, resolution = 0.6, verbose = TRUE)
# Visualise Clustering
DimPlot(annotated_gdF, reduction = "umap", label = TRUE, label.size = 3)
ggsave("gdClustering_umap.pdf", h =6, w=6)
# Using cell selector investigate cells @UMAP-1 ~15
plot <- DimPlot(annotated_gdF)
cells.proliferating <- CellSelector(plot = plot)
# Click done in GUI, then select cells again below
annotated_gdF <- CellSelector(plot = plot, object = annotated_gdF, ident = "8")
annotated_gdF[["integrated_snn_res.0.6_V2"]] <-Idents(object = annotated_gdF)
Idents(annotated_gdF) <- "integrated_snn_res.0.6_V2"
DimPlot(annotated_gdF, reduction = "umap", label = FALSE)
# Check annotated clusters expression of canonical genes
DotPlot(annotated_gdF, features=c("CD14","LYZ", "S100A8", "S100A9", "IL1A", "IL1B", "TNF", "FCER1G", "CTSS", "FCGR3A", "CX3CR1",
                                      "CDKN1C", "IFITM3", "FCER1A", "CD1C", "HLA-DRA", "HLA-DRB1", "CD74", "CLEC4C", "LILRA4", "PLD4",
                                      "SERPINF1", "IL3RA", "CCR7", "GIMAP7", "IL7R", "MAL", "LTB", "NELL2", "PASK", "CD8A", "CD8B",
                                      "KLRB1", "KLRG1", "TRGC1", "TRGC2", "TRDC", "NCAM1", "GNLY", "NKG7", "GZMA", "GZMB", "KLRC1",
                                      "KLRD1", "CD19", "CD79A", "BANK1", "IGHM","IGHD", "MZB1", "CD38", "JCHAIN", "IGHA1", "IGHG1",
                                      "IGHG3", "MKI67", "PCNA", "STMN1", "PCLAF", "CD34", "SPINK2", "SMIM24", "AVP", "PPBP", "ITGB3",
                                      "PF4", "NRGN"),cols=c("white", "red3")) +coord_flip()+ theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
ggsave("gdClustering_Filtered-investigate.pdf", h =16, w=8)
# Based on canonical gene expression c8 &c7 look like artifact of clustering likely pDCs/CD14 monocytes & CD16 monocytes respectively- remove from analysis
keep_idents <- c("0", "1", "2", "3", "4", "5", "6")
annotated_gdF_sub <- subset(annotated_gdF, integrated_snn_res.0.6_V2 %in% keep_idents)
Idents(annotated_gdF_sub) <- "integrated_snn_res.0.6_V3"
DimPlot(annotated_gdF_sub, reduction = "umap", label = FALSE)
# Re-run UMAP
annotated_gdF_sub <- RunUMAP(annotated_gdF_sub, dims = 1:30, reduction = "pca")
# Re-run Determine the K-nearest neighbor graph
annotated_gdF_sub <- FindNeighbors(object = annotated_gdF_sub, dims = 1:30)
annotated_gdF_sub <- FindClusters(object = annotated_gdF_sub, resolution = 0.6, verbose = TRUE)
DimPlot(annotated_gdF_sub, reduction = "umap", label = TRUE, label.size = 3)
ggsave("gdClustering_Filtered-umap.pdf", h =6, w=6)
Idents(annotated_gdF_sub) <- "integrated_snn_res.0.6"
DimPlot(annotated_gdF_sub,reduction = "umap", split.by = "Day")
ggsave("gdClustering_Filtered-umap_day.pdf", h =4, w=8)
DimPlot(annotated_gdF_sub,reduction = "umap", split.by = "Phase")
ggsave("gdClustering_Filtered-umap_phase.pdf", h =4, w=8)
VlnPlot(annotated_gdF_sub, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), group.by = "integrated_snn_res.0.6", pt.size = 0.001)
ggsave("gdClustering_Filtered-QC.pdf", h =4, w=8)
# This code here sets all RNA rows (genes) as all.genes then scales the gene expression across cells/barcodes
DefaultAssay(annotated_gdF_sub) <- "RNA"
all.genes <- rownames(annotated_gdF_sub)
annotated_gdF_sub <- ScaleData(annotated_gdF_sub, features = all.genes)
# find markers for every cluster compared to all remaining cells, report
# only the positive ones
annotated_gdF_sub.markers <- FindAllMarkers(annotated_gdF_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
annotated_gdF_sub.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(annotated_gdF_sub.markers, "gdClustering_TOPDEGbwClusters.csv") 
# Visualisation of Top20 DEG
pdf("gdClustering_Top20_DEGplots.pdf", h=14, w =14)
VlnPlot(annotated_gdF_sub, features = c("YBX3", "RGCC", "CCL4L2","CCL4","KLRB1","LTB","CEBPD","CD8A",
                                    "GZMB","GNLY","AQP3","FGFBP2","GNLY","HLA-DRA","CD74"))
FeaturePlot(annotated_gdF_sub,reduction = "umap",  keep.scale = "feature", features = c("YBX3", "RGCC", "CCL4L2","CCL4","KLRB1","LTB","CEBPD","CD8A",
                                                                                    "GZMB","GNLY","AQP3","FGFBP2","GNLY","HLA-DRA","CD74"))
dev.off()
annotated_gdF_sub.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(annotated_gdF_sub, features = top20$gene) 
# find all markers of cluster 0
pdf("gdClustering_cluster_vlnandFeatureplots.pdf", h = 12, w = 16)
cluster0.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 0, ident.2 = c(1:6),  min.pct = 0.25, only.pos = TRUE)
cluster0.markers <- cluster0.markers[order(cluster0.markers$avg_log2FC, rev(cluster0.markers$p_val_adj), decreasing = TRUE), ]
head(cluster0.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("FGFBP2", "GNLY", "GZMH","GZMB", "NKG7","FCGR3A","KLRD1","LGALS1","KLRD1","CD247","CST7","CCL5", "PRF1", "KLRF1")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("FGFBP2", "GNLY", "GZMH","GZMB", "NKG7","FCGR3A","KLRD1","LGALS1","KLRD1","CD247","CST7","CCL5", "PRF1", "KLRF1"))
# find all markers of cluster 1
cluster1.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 1, ident.2 = c(0,2:6), min.pct = 0.25, only.pos = TRUE)
cluster1.markers <- cluster1.markers[order(cluster1.markers$avg_log2FC, rev(cluster1.markers$p_val_adj), decreasing = TRUE), ]
head(cluster1.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("LTB", "AQP3", "TNFRSF4","ICOS", "CD4","CCR6","CD2","CCR7","IL7R","FLT3LG"))
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("LTB", "AQP3", "TNFRSF4","ICOS", "CD4","CCR6","CD2","CCR7","IL7R","FLT3LG"))
# find all markers of cluster 2
cluster2.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 2, ident.2 = c(0:1,3:6), min.pct = 0.25, only.pos = TRUE)
cluster2.markers <- cluster2.markers[order(cluster2.markers$avg_log2FC, rev(cluster2.markers$p_val_adj), decreasing = TRUE), ]
head(cluster2.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("CEBPD", "CD8A", "KLRB1","FAM43A", "NCR3","GZMK","AQP3","IL4I1","IL7R","FURIN")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("CEBPD", "CD8A", "KLRB1","FAM43A", "NCR3","GZMK","AQP3","IL4I1","IL7R","FURIN"))
# find all markers of cluster 3
cluster3.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 3, ident.2 = c(0:2,4:6), min.pct = 0.25, only.pos = TRUE)
cluster3.markers <- cluster3.markers[order(cluster3.markers$avg_log2FC, rev(cluster3.markers$p_val_adj), decreasing = TRUE), ]
head(cluster3.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("CCL4L2", "CCL4", "CCL3","IFNG", "TNF","EGR1","XCL2","CD69","NFKBID", "BTG2", "IER2","PMAIP1")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("CCL4L2", "CCL4", "CCL3","IFNG", "TNF","EGR1","XCL2","CD69","NFKBID", "BTG2", "IER2","PMAIP1"))
# find all markers of cluster 4
cluster4.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 4, ident.2 = c(0:3,5:6), min.pct = 0.25, only.pos = TRUE)
cluster4.markers <- cluster4.markers[order(cluster4.markers$avg_log2FC, rev(cluster4.markers$p_val_adj), decreasing = TRUE), ]
head(cluster4.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("YBX3", "RGCC", "ID1","VIM", "IER5L","CREN","CD7","IFRD1","CXCR3")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("YBX3", "RGCC", "ID1","VIM", "IER5L","CREN","CD7","IFRD1","CXCR3"))
# find all markers of cluster 5
cluster5.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 5, ident.2 = c(0:4,6), min.pct = 0.25, only.pos = TRUE)
cluster5.markers <- cluster5.markers[order(cluster5.markers$avg_log2FC, rev(cluster5.markers$p_val_adj), decreasing = TRUE), ]
head(cluster5.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("KLRB1", "CCR6","LTB", "CD8A","TNF", "IER3","IL7R","NCR3","CD69","EGR1", "GZMK", "AQP3")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("KLRB1", "CCR6","LTB", "CD8A","TNF", "IER3","IL7R","NCR3","CD69","EGR1", "GZMK", "AQP3"))
# find all markers of cluster 6
cluster6.markers <- FindMarkers(annotated_gdF_sub, ident.1 = 6, ident.2 = c(0:5), min.pct = 0.25, only.pos = TRUE)
cluster6.markers <- cluster6.markers[order(cluster6.markers$avg_log2FC, rev(cluster6.markers$p_val_adj), decreasing = TRUE), ]
head(cluster6.markers, n = 20)
VlnPlot(annotated_gdF_sub, features = c("HLA-DRA", "CD74","COTL1", "CD8B","TRAC", "HLA-DRB1","HLA-DQB1","GZMK","CXCR3","HLA-DPA1", "HLA-DRB5", "HLA-DPB1","IL32")) 
FeaturePlot(annotated_gdF_sub, reduction = "umap", features = c("HLA-DRA", "CD74","COTL1", "CD8B","TRAC", "HLA-DRB1","HLA-DQB1","GZMK","CXCR3","HLA-DPA1", "HLA-DRB5", "HLA-DPB1","IL32"))
dev.off()
# c2 and c5 look similar
GD_cluster_markers2V5_res0.6 <- FindMarkers(annotated_gdF, ident.1 = "2", ident.2 = "5", grouping.var = "sample")
GD_cluster_markers2V5_res0.6 <- GD_cluster_markers2V5_res0.6[order(GD_cluster_markers2V5_res0.6$avg_log2FC, rev(GD_cluster_markers2V5_res0.6$p_val_adj), decreasing = TRUE), ]
head(GD_cluster_markers2V5_res0.6, n = 20)
tail(GD_cluster_markers2V5_res0.6, n = 20)
# Difference is mostly inflammatory gene expression but difference is not great plan to combine
# Write CSV for records
write.csv(cluster0.markers, "cluster0_markers.csv")
write.csv(cluster1.markers, "cluster1_markers.csv")
write.csv(cluster2.markers, "cluster2_markers.csv")
write.csv(cluster3.markers, "cluster3_markers.csv")
write.csv(cluster4.markers, "cluster4_markers.csv")
write.csv(cluster5.markers, "cluster5_markers.csv")
write.csv(cluster6.markers, "cluster6_markers.csv")
write.csv(GD_cluster_markers2V5_res0.6, "GD_cluster_markers2V5_res0.6.csv")
# based upon above cluster investigation above and gene list from paper https://pubmed.ncbi.nlm.nih.gov/33893173/ ; in particular Naive and Type3
#Cytotoxic GNLY	GZMB	GZMH	NKG7	FCGR3A FGFBP2
# Inflammatory "CCL4l2", CCL4 "CCL3","IFNG", "TNF"
#antigen-presenting "HLA-DRA", "CD74", "HLA-DQB1", "HLA-DRB1",""HLA-DPA1", "HLA-DPB1","HLA-DRB5"
#transitional YBX3, ER5L ypel5  ifrd1 akna(cd40 related)
# type 3 - il7r klrb1 ncr3 RORA based on https://www.science.org/doi/epdf/10.1126/sciimmunol.abf0125
# naive - ltb ccr7 aqp3 cxcr3 based on https://www.science.org/doi/epdf/10.1126/sciimmunol.abf0125
Idents(annotated_gdF_sub)<-"integrated_snn_res.0.6"
DimPlot(annotated_gdF_sub, reduction = "umap")
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 == 0] <- "Cytotoxic"
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 == 2] <- "Naive"
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 %in% c(1,5)] <- "Type3"
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 == 3] <- "Inflammatory"
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 == 4] <- "Transitional"
annotated_gdF_sub@meta.data$res0.6_collapsed[annotated_gdF_sub@meta.data$integrated_snn_res.0.6 == 6] <- "Antigen_presenting"
# Select cluster annotation
Idents(annotated_gdF_sub) <- "res0.6_collapsed"
# Visualise cluster annotation with UMAP
pdf("annotatedGD_umap.pdf", h = 8, w = 8)
DimPlot(annotated_gdF_sub, reduction = "umap")
DimPlot(annotated_gdF_sub, reduction = "umap", split.by = "Day", ncol = 4)
DimPlot(annotated_gdF_sub, reduction = "umap", split.by = "sample", ncol = 4)
dev.off()
# Visualize TOP DEG gene expressions with violin plots comparing clusters
DefaultAssay(annotated_gdF_sub) <- "RNA"
GD_cluster_markers_res0.6_collapsed <- FindAllMarkers(annotated_gdF_sub, grouping.var = "res0.6_collapsed")
GD_cluster_markers_res0.6_collapsed <- GD_cluster_markers_res0.6_collapsed %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(GD_cluster_markers_res0.6_collapsed, "GD_cluster_markers_res0.6_collapsed.csv")
pGD_cluster_markers_res0.6_collapsed <- subset(GD_cluster_markers_res0.6_collapsed, p_val_adj < 0.0500)
write.csv(pGD_cluster_markers_res0.6_collapsed, "pGD_cluster_markers_res0.6_collapsed.csv")
# Order clusters in meta.data
annotated_gdF_sub@meta.data$res0.6_collapsed <- factor(annotated_gdF_sub@meta.data$res0.6_collapsed, levels = c("Cytotoxic","Inflammatory","Antigen_presenting","Transitional","Type3", "Naive"))
# Dotplot of Avg Expression
DotPlot(annotated_gdF_sub, features=c("GNLY","GZMB","GZMH", "NKG7", "FCGR3A","FGFBP2",
                                  "CCL4L2", "CCL4", "CCL3", "IFNG", "TNF", 
                                  "HLA-DRA", "CD74", "HLA-DQB1", "HLA-DPA1", 
                                  "IER5L", "YPEL5", "IFRD1", "IL21R",
                                  "KLRB1", "NCR3", "RORA","IL7R",
                                  "LTB", "AQP3", "CCR7"),
        cols=c("white", "red3"), col.min = 0, col.max = 2.5, group.by = "res0.6_collapsed") +
  theme_bw(base_size=14)+
  theme(axis.text.x = element_text(angle = 90,  vjust = 1, hjust=1))+
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
ggsave("annotatedGD_dotplotID.pdf", h = 2.71, w = 9.02)
# Add age meta.data
annotated_gdF_sub@meta.data$age[annotated_gdF_sub@meta.data$donor %in% c("Adult1","Adult2","Adult3", "Control1", "Control2")] <- "Adult"
annotated_gdF_sub@meta.data$age[annotated_gdF_sub@meta.data$donor %in% c("Child1","Child2","Child3")] <- "Child"
# Subcluster frequencies across days line plots
countcluster <- table(annotated_gdF_sub@meta.data$sample,annotated_gdF_sub@meta.data$res0.6_collapsed)
write.csv(countcluster, "countcluster.csv")
countcluster <- read.csv("countcluster.csv")
colnames(countcluster)[which(colnames(countcluster) == "X")] <- "sample"
# Assuming the Seurat object is named 'seurat_obj'
metadata_df <- annotated_gdF_sub@meta.data %>%
  select(Day, donor, age, sample) %>% # Select the columns you want to add
  distinct() %>% # Remove duplicate rows
  rename_all(~tolower(.)) # Rename the columns to lowercase
# Add the columns to the original table
countcluster_with_metadata <- countcluster %>%
  left_join(metadata_df, by = "sample")
countcluster_with_metadata$Total <- rowSums(countcluster_with_metadata[, 2:7])
# Define levels and comparisons
countcluster_with_metadata$donor <- factor(countcluster_with_metadata$donor)
countcluster_with_metadata$day <- factor(countcluster_with_metadata$day, levels = c("Day0","Day7","Day28","Control"))
ST <- list(c("Day0","Day7"), c("Day0", "Day28"),c("Day0","Control"))
# create a vector of column names to loop over
cols <- c("Naive", "Type3", "Transitional", "Antigen_presenting", "Inflammatory", "Cytotoxic")
# loop over each column and create a graph
for (col in cols) {
  # convert columns to factor if necessary
  if (col %in% c("donor", "day")) {
    countcluster_with_metadata[[col]] <- factor(countcluster_with_metadata[[col]])
  }
  # create graph
  graph <- ggplot(countcluster_with_metadata, aes(x = day, y = ((.data[[col]] / Total) * 100))) +
    geom_line(aes(group = donor)) +
    geom_point(aes(shape = age), size = 2) +
    theme_bw() +
    scale_y_continuous(trans = "log10", limits = c(0.5, 100), expand = expansion(mult = c(0.1, 0.1))) +
    stat_compare_means(comparisons = ST, paired = FALSE, method = "wilcox.test", tip.length = 0, size = 2) +
    labs(y = paste(col, "Log (Proportions of cell type)"), x = "") +
    theme(axis.text.x = element_text(size = 6, face = "bold"),
          axis.text.y = element_text(size = 6, face = "bold"),
          axis.title.y = element_text(size = 6, face = "bold"),
          axis.title.x = element_text(size = 6, face = "bold"),
          strip.text = element_text(size = 6, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 4))

  # save graph as a pdf file
  ggsave(paste0(col, "_lineplot_cluster.pdf"), graph, height = 2, width = 2.5)
}
# Begin DEG investigation of Clusters across day0 and day28
annotated_gdF_sub@meta.data$cluster_day <- paste(annotated_gdF_sub@meta.data$res0.6_collapsed, annotated_gdF_sub@meta.data$Day, sep="_")
Idents(annotated_gdF_sub) <- "cluster_day"
DefaultAssay(annotated_gdF_sub) <- "RNA"
Naive_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Naive_Day0", ident.2 = "Naive_Day28", only.pos = TRUE)
Naive_0v28 <- Naive_0v28[order(Naive_0v28$avg_log2FC, rev(Naive_0v28$p_val_adj), decreasing = TRUE), ]
Naive_0v28 <- Naive_0v28 %>% rownames_to_column(var = "gene")
Naive_0v28['cluster'] <- "Naive"
Type3_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Type3_Day0", ident.2 = "Type3_Day28", only.pos = TRUE)
Type3_0v28 <- Type3_0v28[order(Type3_0v28$avg_log2FC, rev(Type3_0v28$p_val_adj), decreasing = TRUE), ]
Type3_0v28 <- Type3_0v28 %>% rownames_to_column(var = "gene")
Type3_0v28['cluster'] <- "Type3"
Transitional_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Transitional_Day0", ident.2 = "Transitional_Day28", only.pos = TRUE)
Transitional_0v28 <- Transitional_0v28[order(Transitional_0v28$avg_log2FC, rev(Transitional_0v28$p_val_adj), decreasing = TRUE), ]
Transitional_0v28 <- Transitional_0v28 %>% rownames_to_column(var = "gene")
Transitional_0v28['cluster'] <- "Transitional"
Inflammatory_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Inflammatory_Day0", ident.2 = "Inflammatory_Day28", only.pos = TRUE)
Inflammatory_0v28 <- Inflammatory_0v28[order(Inflammatory_0v28$avg_log2FC, rev(Inflammatory_0v28$p_val_adj), decreasing = TRUE), ]
Inflammatory_0v28 <- Inflammatory_0v28 %>% rownames_to_column(var = "gene")
Inflammatory_0v28['cluster'] <- "Inflammatory"
Cytotoxic_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Cytotoxic_Day0", ident.2 = "Cytotoxic_Day28", only.pos = TRUE)
Cytotoxic_0v28 <- Cytotoxic_0v28[order(Cytotoxic_0v28$avg_log2FC, rev(Cytotoxic_0v28$p_val_adj), decreasing = TRUE), ]
Cytotoxic_0v28 <- Cytotoxic_0v28 %>% rownames_to_column(var = "gene")
Cytotoxic_0v28['cluster'] <- "Cytotoxic"
Antigen_presenting_0v28 <- FindMarkers(annotated_gdF_sub, ident.1 = "Antigen_presenting_Day0", ident.2 = "Antigen_presenting_Day28", only.pos = TRUE)
Antigen_presenting_0v28 <- Antigen_presenting_0v28[order(Antigen_presenting_0v28$avg_log2FC, rev(Antigen_presenting_0v28$p_val_adj), decreasing = TRUE), ]
Antigen_presenting_0v28 <- Antigen_presenting_0v28 %>% rownames_to_column(var = "gene")
Antigen_presenting_0v28['cluster'] <- "Antigen_presenting"
# remove degs with p values greater than 0.05
pNaive_0v28 <- subset(Naive_0v28, p_val_adj < 0.0500)
pType3_0v28 <- subset(Type3_0v28, p_val_adj < 0.0500)
pTransitional_0v28 <- subset(Transitional_0v28, p_val_adj < 0.0500)
pInflammatory_0v28 <- subset(Inflammatory_0v28, p_val_adj < 0.0500)
pCytotoxic_0v28 <- subset(Cytotoxic_0v28, p_val_adj < 0.0500)
pAntigen_presenting_0v28 <- subset(Antigen_presenting_0v28, p_val_adj < 0.0500)
write.csv(Naive_0v28, "Naive_0v28_DEG.csv")
write.csv(Type3_0v28, "Type3_0v28_DEG.csv")
write.csv(Transitional_0v28, "Transitional_0v28_DEG.csv")
write.csv(Inflammatory_0v28, "Inflammatory_0v28_DEG.csv")
write.csv(Cytotoxic_0v28, "Cytotoxic_0v28_DEG.csv")
write.csv(Antigen_presenting_0v28, "Antigen_presenting_0v28_DEG.csv")
write.csv(pNaive_0v28, "pNaive_0v28_DEG.csv")
write.csv(pType3_0v28, "pType3_0v28_DEG.csv")
write.csv(pTransitional_0v28, "pTransitional_0v28_DEG.csv")
write.csv(pInflammatory_0v28, "pInflammatory_0v28_DEG.csv")
write.csv(pCytotoxic_0v28, "pCytotoxic_0v28_DEG.csv")
write.csv(pAntigen_presenting_0v28, "pAntigen_presenting_0v28_DEG.csv")
# Create data frame of DEGs with adjusted P values less than 0.0500 combined
pCombined_0v28 <- bind_rows(pNaive_0v28, pType3_0v28, pTransitional_0v28, pInflammatory_0v28, pCytotoxic_0v28, pAntigen_presenting_0v28)
# Remove genes that start with "RP" or "MT" from the "gene" column
pCombined_0v28_noMTRP <- pCombined_0v28[!grepl("^RP|^MT", pCombined_0v28$gene), ]
# ID top 20 genes from each column with increased expression at day 0
top20.pCombined_0v28_noMTRP <- pCombined_0v28_noMTRP %>%
  group_by(cluster)%>%
  top_n(n=20, wt = avg_log2FC)
head(pCombined_0v28_noMTRP)
# Define Paramaters for dotplot
DEGdf <- pCombined_0v28
DEGdf$log_p <- -log10(DEGdf$p_val_adj)
DEGdf$log_p_cat[DEGdf$log_p <2] <- "<2"
DEGdf$log_p_cat[DEGdf$log_p >2 & DEGdf$log_p<=10] <- "2-10"
DEGdf$log_p_cat[DEGdf$log_p >10 ] <- ">10"
DEGdf$log_p_cat <- factor(DEGdf$log_p_cat, levels = c("<2","2-10", ">10"))
DEGdf$cluster <- factor(DEGdf$cluster, levels = c("Cytotoxic","Inflammatory","Antigen_presenting","Transitional","Type3", "Naive"))
DEGdf2 <- DEGdf
DEGdf2 <- DEGdf2 %>%
  subset(gene %in% c ("VIM",	"CCL4L2",	"BCL2A1",	"CCL4",	"GZMB",	"SMAP2",	"NFKBIZ",	"HLA-DQB1",	"IL7R",	"LAG3",	"GRASP",	"GPR183",	"ISG20",	"RNF19A",	"ICOS", "STAT5A", "IRF4", "SELL","FASLG","NFKB1","PIM3","RELB","DUSP10",	"FAM177A1","OTULIN","CYTIP","ZNF267","FURIN","RAB8B","MBNL1","IFITM2","CD82","SNX9","PLAAT4","TANK","NDUFA13","SYNJ2","TNIP2","MARCHF2","GNG2","EFHD2","CD74","NPDC1","CRTAM","CPT1A","CIB1","SRGN","TMSB10","TXNIP","H1-10","FKBP11","CD2","CLEC2B","PDE4D","LUZP1","FTH1","HIF1A","RNF145","IFI16","JAK3","TNFRSF4","EVA1B","PPP1R14B","BATF","PDCD1","HAVCR2"))
DEGdf2$gene <- factor(DEGdf2$gene, levels = c ("VIM",	"CCL4L2",	"BCL2A1",	"CCL4",	"GZMB",	"SMAP2",	"NFKBIZ",	"HLA-DQB1",	"IL7R",	"LAG3",	"GRASP",	"GPR183",	"ISG20",	"RNF19A",	"ICOS", "STAT5A", "IRF4", "SELL","FASLG","NFKB1","PIM3","RELB","DUSP10",	"FAM177A1","OTULIN","CYTIP","ZNF267","FURIN","RAB8B","MBNL1","IFITM2","CD82","SNX9","PLAAT4","TANK","NDUFA13","SYNJ2","TNIP2","MARCHF2","GNG2","EFHD2","CD74","NPDC1","CRTAM","CPT1A","CIB1","SRGN","TMSB10","TXNIP","H1-10","FKBP11","CD2","CLEC2B","PDE4D","LUZP1","FTH1","HIF1A","RNF145","IFI16","JAK3","TNFRSF4","EVA1B","PPP1R14B","BATF","PDCD1","HAVCR2"))
# Create dotplot from selected genes
DEGdf_dotplot <- ggplot(DEGdf2, aes(reorder(gene, cluster), cluster))+
  geom_point(aes(colour=avg_log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low = "white", high = "red")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))
DEGdf_dotplot
# Use Venn diagram webtool to identify shared DEGS for 0 v 28 between clusters, export text file and input here in code
# from GD_DEG_upset <- read_csv("venn_DEG_0v28_INPUT.csv")
GD_DEG_upset<-c("Antigen_presenting&Cytotoxic&Inflammatory&Naive&Transitional&Type3" =3,
                "Antigen_presenting&Cytotoxic&Inflammatory&Naive&Type3" =1,
                "Antigen_presenting&Cytotoxic&Naive&Transitional&Type3" =1,
                "Cytotoxic&Inflammatory&Naive&Transitional&Type3" =1,
                "Antigen_presenting&Cytotoxic&Inflammatory&Naive" =2,
                "Antigen_presenting&Cytotoxic&Transitional&Type3" =1,
                "Antigen_presenting&Cytotoxic&Naive&Type3" =2,
                "Cytotoxic&Inflammatory&Transitional&Type3" =1,
                "Cytotoxic&Naive&Transitional&Type3" =1,
                "Antigen_presenting&Cytotoxic&Naive" =1,
                "Antigen_presenting&Naive&Type3" =3,
                "Cytotoxic&Inflammatory&Transitional" =1,
                "Cytotoxic&Inflammatory&Naive" =1,
                "Cytotoxic&Inflammatory&Type3" =2,
                "Cytotoxic&Transitional&Type3" =1,
                "Cytotoxic&Naive&Type3" =3,
                "Inflammatory&Transitional&Type3" =1,
                "Inflammatory&Naive&Type3" =1,
                "Antigen_presenting&Cytotoxic" =3,
                "Cytotoxic&Inflammatory" =3,
                "Cytotoxic&Transitional" =5,
                "Cytotoxic&Naive" =6,
                "Cytotoxic&Type3" =6,
                "Inflammatory&Type3" =3,
                "Transitional&Type3" =1,
                "Naive&Type3" =7,
                "Naive" =22,
                "Type3" =51,
                "Transitional" =6,
                "Antigen_presenting" =3,
                "Inflammatory" =20,
                "Cytotoxic" =34
)

pdf("UpsetR_DEG_0v28.pdf", h=8, w =14)
upset(fromExpression(GD_DEG_upset),
      nintersects = 40,
      nsets = 6,
      order.by = c("freq"),
      decreasing = T,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1,
)
dev.off()
#Save Final version of seurat object
saveRDS(annotated_gdF_sub, "gdFINAL.rds")
