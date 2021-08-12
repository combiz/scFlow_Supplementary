# Preparation of Zhou_et_al_2020_mouse data for scflow
# Nurun 16/02/2021

# Dataset was downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788


library(dplyr)
library(tidyverse)
library(DropletUtils)
library(Matrix)

setwd("/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/ximerakis_et_al/")
#source("/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/scflow_resources/scFlow/R/read_sparse_matrix.R")

source_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/ximerakis_et_al/ximerakis2019_sc_data/"
dest_dir <- "filtered_matrix"

mat <- read.delim("ximerakis2019_sc_data/ximerakis2019.raw.counts.txt")

# qs::qsave(mat, file = "ximerakis2019.raw.counts.qs")
ensemble_mapping <- read.delim("/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/scflow_resources/ensembl_mappings_mouse.tsv")

# The gene names were converted to Ensembl gene id

# mat <- qs::qread("ximerakis2019.raw.counts.qs")
mat$ensembl_id <- ensemble_mapping[match(mat$X, ensemble_mapping$external_gene_name), "ensembl_gene_id"]
mat <- mat %>%
  dplyr::filter(!is.na(ensembl_id))
rownames(mat) <- mat$ensembl_id

# qs::qsave(mat, file = "ximerakis2019.raw.counts.qs")
metadata <- read.delim("ximerakis2019_sc_data/ximerakis2019.metadata.txt")

# write the files into respective sample directory in cellranger format

file_list <- list.files("raw_data", full.names = TRUE)
for(file in file_list){
  tmp_mat <- read.delim(file)
  cell_barcode <- colnames(tmp_mat)
  col_idx <- match(cell_barcode, colnames(mat))
  col_idx <- col_idx[!is.na(col_idx)]
  tmp_mat <- mat[ , col_idx]
  tmp_mat <- as.matrix(tmp_mat)
  tmp_mat <- as(tmp_mat, "dgCMatrix")
  sample_name <- purrr::map_chr(basename(file), ~ strsplit(., "_")[[1]][1])
  dest_path <- file.path(dest_dir, sample_name)
  dir.create(dest_path, recursive = TRUE)
  DropletUtils::write10xCounts(path = dest_path,
                 x = tmp_mat,
                 barcodes = colnames(tmp_mat),
                 gene.id = rownames(tmp_mat),
                 gene.symbol = rownames(tmp_mat),
                 gene.type = "Gene Expression",
                 overwrite = TRUE,
                 type = "auto",
                 genome = "unknown",
                 version = "3")
}


# SampleSheet.tsv and Manifest.txt generation

data_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/ximerakis_et_al/filtered_matrix"
manifest <- data.frame(key = ids::proquint(n = length(dir(data_dir)), 
                                           n_words = 1L, use_cache = TRUE, use_openssl = FALSE),
                       filepath = dir(data_dir, full.names = TRUE))


samplesheet <- read.delim("GSE129788_samplesheet.txt")
samplesheet <- samplesheet %>%
  dplyr::rename(individual = Sample_title,
                sample_id = Sample_geo_accession,
                age = Sample_characteristics_ch1) %>%
  dplyr::mutate(individual = purrr::map_chr(as.character(individual), ~ strsplit(., ":", fixed = TRUE)[[1]][1]),
                age = purrr::map_chr(as.character(age), ~ strsplit(., " ", fixed = TRUE)[[1]][2])) %>%
  dplyr::mutate(age = paste0(age, "m")) %>%
  dplyr::arrange(sample_id)

idx <- match(samplesheet$sample_id, basename(as.character(manifest$filepath)))
manifest <- manifest[idx, ]

samplesheet <- samplesheet %>%
  dplyr::mutate(manifest = manifest$key) %>%
  dplyr::select(c(manifest, dplyr::everything()))

samplesheet <- samplesheet %>%
  dplyr::mutate(group = case_when(age == "2-3m" ~ "Young",
                                  age == "21-22m" ~ "Old"))


dir.create("refs")
write.table(manifest, file = "refs/Manifest.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(samplesheet, file = "refs/SampleSheet.tsv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# save(manifest, metadata, samplesheet, idx, file = "ximerakis_data_processing.rda")
