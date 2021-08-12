# Preparation of Zhou_et_al_2020_mouse data for scflow
# Nurun 16/02/2021

# Dataset was downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140511

library(dplyr)
setwd("/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_mouse")
source_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_mouse/raw_data/"
dest_dir <- "filtered_matrix"

sample_id <- list.files(path = source_dir)
sample_id <- purrr::map_chr(sample_id, ~ strsplit(., "_")[[1]][1]) %>%
  unique()

# move the files into respective sample directory

for (id in sample_id) {
  file_list <- list.files(path = source_dir,
                          pattern = id,
                          full.names = TRUE)
  
  dest_path <- file.path(dest_dir, id)
  dir.create(path = dest_path, recursive = TRUE)
  
  file.copy(from = file_list, 
            to = dest_path)
}


comsub <- function(x) {
  # sort the vector
  x <- sort(x)
  # split the first and last element by character
  d_x <- strsplit(x[c(1,length(x))],"")
  # search for the first not common element and so, get the last matching one
  der_com <- match(FALSE,do.call("==",d_x))-1
  # if there is no matching element, return an empty vector, else return the common part
  ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))
}

# Renaming the file names to feature/barcode/matrix structure

for (id in sample_id) {
  file_path <- file.path("filtered_matrix", id)
  file_list <- list.files(path = file_path, full.names = TRUE)
  prefix <- basename(file_list) %>% comsub()
  for(file in file_list) {
    new_file <- gsub(prefix, "", file)
    file.rename(from = file,
                to = new_file)
  }
}

for(id in sample_id[1]){
  file <- list.files(path = file.path("filtered_matrix", id), pattern = "genes", full.names = TRUE)
  new_file <- gsub("genes", "features", file)
  file.rename(from = file,
              to = new_file)
}


# SampleSheet.tsv and Manifest.txt generation

data_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_mouse/filtered_matrix"
manifest <- data.frame(key = ids::proquint(n = length(dir(data_dir)), 
                                           n_words = 1L, use_cache = TRUE, use_openssl = FALSE),
                       filepath = dir(data_dir, full.names = TRUE))

metadata <- read.delim("GPL24247_samplesheet.txt")
metadata <- metadata %>%
  dplyr::rename(sample_title = Sample_title,
                sample_id = Sample_geo_accession,
                genotype = Sample_characteristics_ch1,
                age = Sample_characteristics_ch1.1) %>%
  dplyr::mutate(genotype = purrr::map_chr(as.character(genotype), ~ strsplit(., " ", fixed = TRUE)[[1]][3]),
                age = purrr::map_chr(as.character(age), ~ strsplit(., " ", fixed = TRUE)[[1]][2])) %>%
  dplyr::mutate(genotype = gsub("[0-9]$", "", SampleSheet$individual),
                age = paste0(age, "m") )
  # dplyr::filter(genotype %in% c("wt", "5XFAD"))

metadata <- metadata[order(metadata$sample_id), ]

idx <- match(metadata$sample_id, basename(as.character(manifest$filepath)))
manifest <- manifest[idx, ]
samplesheet <- data.frame(manifest = manifest$key,
                         sample_id = metadata$sample_id,
                         individual = metadata$sample_title,
                         group = metadata$genotype,
                         age = metadata$age)


dir.create("refs")
write.table(manifest, file = "refs/Manifest.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(samplesheet, file = "refs/SampleSheet.tsv", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# save.image("zhou_mouse_data_processing.rda")

