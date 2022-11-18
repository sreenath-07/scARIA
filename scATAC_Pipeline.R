# library(Signac)
# library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
# library(patchwork)
library(hdf5r)

install.packages("Seurat")
library(Seurat)
install.packages(c('dplyr','patchwork'))
library(dplyr)
library(patchwork)
install.packages("Signac")
library(Signac)

require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)

##### READING 10x Matrix File Format

# peak-bc matrix
mex_dir_path <- 'C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC'

mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '\\')
feature_path <- paste(mex_dir_path, "peaks.bed", sep = '\\')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '\\')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

# # tf-bc matrix
# mex_dir_path <- "/opt/sample345/outs/filtered_tf_bc_matrix"
# 
# mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
# feature_path <- paste(mex_dir_path, "motifs.tsv", sep = '/')
# barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')
# 
# features <- readr::read_tsv(feature_path, col_names = c('feature', 'common_name'))
# barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
# 
# mtx <- Matrix::readMM(mtx_path) %>%
#   magrittr::set_rownames(features$feature) %>%
#   magrittr::set_colnames(barcodes$barcode)


expression_matrix <- Read10X(data.dir = mex_dir_path)

##### Creating Chromatin Assay and Seurat Object 10x Matrix File Format

chrom_assay <- CreateChromatinAssay(
  counts = mtx,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = features,
  min.cells = 10,
  min.features = 200
)

bmmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)


# seurat_object = CreateSeuratObject(counts = mtx)

bmmc[['peaks']]
granges(bmmc)

bmmc <- RunTFIDF(bmmc)
bmmc <- FindTopFeatures(bmmc, min.cutoff = 'q0')
bmmc <- RunSVD(bmmc)
DepthCor(bmmc)

bmmc <- RunUMAP(object = bmmc, reduction = 'lsi', dims = 2:30)
bmmc <- FindNeighbors(object = bmmc, reduction = 'lsi', dims = 2:30)
bmmc <- FindClusters(object = bmmc, verbose = FALSE, algorithm = 3)
DimPlot(object = bmmc, label = TRUE) + NoLegend()


#-------------

set.seed(1234)

data_dir <- 'C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)

counts <- Read10X_h5(filename = "C:\\Users\\ssrikrishnan6\\scATAC\\scATAC_GEO_Dataset\\CLL_BMMC_ATAC\\GSM4829412_cll_atac_filtered_matrix.mtx.gz")
metadata <- read.csv(
  file = "../vignette_data/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)