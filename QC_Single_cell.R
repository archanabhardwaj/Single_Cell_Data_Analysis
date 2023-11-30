## 
#   This code highlights the basics for Seurat based single cell data analysis pipeline, including 
#    -> loading data, 
#    -> quality control, 
#    -> visualization of QC metrics, 
#    -> data filtration, and 
#   -> downstream analysis steps such as integration, dimensionality reduction, clustering, and visualization using UMAP. 


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

## Read cell ranger datafiles #### 
files_all <- list.dirs(path = "SINGLE_CELL/Rh_PBMCs/Cellranger", full.names = TRUE, recursive = FALSE)
files <- grep(pattern = "-L1", files_all, value = TRUE)
object_sample_id <- list()

for( i in 1:8){
sample_id <- files[i]
read <- paste(sample_id,"outs/filtered_feature_bc_matrix",sep="/")
expression_matrix <- Read10X(data.dir = read)
expression_matrix_n <- expression_matrix[- grep("mm10__", rownames(expression_matrix)),]
name <- gsub("/SINGLE_CELL/Rh_PBMCs/Cellranger/", "", files[i])
object_sample_number <- name
object_sample_id[i] <- CreateSeuratObject(counts = expression_matrix_n, project = object_sample_number, min.cells = 3, min.features = 200)
}

saveRDS(object_sample_id, "Rh_seurat_proceeded.RData")

###### perform QC of the data ####
qc <- readRDS(file="Rh_seurat_proceeded.RData")
object_sample_id <- qc 

alldata <- merge(object_sample_id[[1]],c(object_sample_id[[2]],object_sample_id[[3]],object_sample_id[[4]],object_sample_id[[5]],
object_sample_id[[6]],object_sample_id[[7]],object_sample_id[[8]]),add.cell.ids = c("532","533","534","535","536", "537","538", "539"))
gene_Listt <- rownames(alldata)
gene_List <- data.frame(gene_Listt)
gene_List_f <- gsub("GRCh38-", "", gene_List$gene_Listt) 
alldataf <- alldata
alldataf@assays$RNA@counts@Dimnames[[1]] <- gene_List_f 
alldataf@assays$RNA@data@Dimnames[[1]] <- gene_List_f 

## Get mitochrondial gene percentage per sample # 
alldataf <- PercentageFeatureSet(alldataf, "^MT-", col.name = "percent_mito")
total_counts_per_cell <- colSums(alldataf@assays$RNA@counts)
mito_genes <- rownames(alldataf)[grep("MT-", rownames(alldataf))]
alldataf$percent_mito <- colSums(alldataf@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
metadata$cells <- rownames(metadata)
metadata$seq_folder =  metadata$orig.ident 
metadata$nCount_RNA =  metadata$nUMI 
metadata$nFeature_RNA =  metadata$nGene        
alldataf$log10GenesPerUMI <- log10(alldataf$nFeature_RNA) / log10(alldataf$nCount_RNA)
alldataf$mitoRatio <- PercentageFeatureSet(object = alldataf, pattern = "^MT-")
alldataf$mitoRatio <- alldataf@meta.data$mitoRatio / 100
metadata <- alldataf@meta.data
mycolors = c(brewer.pal(name="Dark2", n = 8))

tiff("cell_count_per_sample.tiff", units="in", width=9, height=9, res=300)
metadata %>% 
  	ggplot(aes(x=orig.ident, fill=orig.ident)) + 
   geom_bar() +
  scale_color_manual(values = mycolors[1:8])+
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells") 
dev.off()

metadata %>% 
  	ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
tiff("genes_per_cell.tiff", units="in", width=9, height=9, res=300)
metadata %>% 
  	ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
tiff("genes_per_cell_boxplot.tiff", units="in", width=9, height=9, res=300)
metadata %>% 
  	ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")
  	dev.off()
  	
# Visualize the distribution of mitochondrial gene expression detected per cell
tiff("mitochondrial_gene_per_cell.tiff", units="in", width=9, height=9, res=300)
metadata %>% 
  	ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.25)

dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
tiff("correlation_gene_nUMIs.tiff", units="in", width=9, height=9, res=300)

metadata %>% 
  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~orig.ident)
dev.off()

## Apply filteration ### 



################### data filteration ######################################
filtered_seurat <- subset(x = alldataf, 
                         subset= (nCount_RNA >= 500) & 
                           (nFeature_RNA >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.25))
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seuratf <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data

tiff("mitochondrial_gene_per_cell_filtered.tiff", units="in", width=9, height=9, res=300)
# Visualize the distribution of mitochondrial gene expression detected per cell
 metadata_clean %>% 
  	ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() ## +
  	#geom_vline(xintercept = 0.25)
dev.off()

# Visualize the distribution of genes detected per cell via boxplot
tiff("genes_per_cell_boxplot_filtered.tiff", units="in", width=9, height=9, res=300)
metadata_clean %>% 
  	ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")
  	dev.off()


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
tiff("correlation_gene_nUMIs_filtered.tiff", units="in", width=9, height=9, res=300)
metadata_clean %>% 
  	ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~orig.ident)
dev.off()

write.table(file="metadata_clean",metadata_clean)
save(alldataf, file="alldataf_noqc.RData")
save(filtered_seurat, file="seurat_filtered.RData")



