library(Seurat)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]

sample_dirs <- list.dirs(input_path, recursive = FALSE)

for (sample_dir in sample_dirs) {
  seurat_obj <- Read10X(data.dir = file.path(sample_dir, "outs", "filtered_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = basename(sample_dir), min.cells = 3, min.features = 200)

  # QC
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

  saveRDS(seurat_obj, file = file.path(output_path, paste0(basename(sample_dir), "_seurat.rds")))
}
