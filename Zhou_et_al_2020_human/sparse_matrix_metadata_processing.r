# Preparation of Zhou_et_al_2020_human data for scflow
# Nurun 16/02/2021

# Dataset was downloaded from https://www.synapse.org/#!Synapse:syn22264146
# Supplementary files were downloaded from https://www.nature.com/articles/s41591-019-0695-9#Abs1

setwd("/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_human")
##### file formatting script
source_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_human/expression_matrix"

# move the files into respective sample directory
sample_id <- list.files(source_dir)
sample_id <- unique(sample_id)

for (id in sample_id) {
  file_list <- list.files(path = source_dir,
                          pattern = id,
                          full.names = TRUE)
  
  sample_path <- file.path(source_dir, id)
  dir.create(path = sample_path, recursive = TRUE)
  
  file.copy(from = file_list, 
            to = sample_path)
}

# Renaming the file names to feature/barcode/matrix structure

sample_id <- list.files(source_dir, full.names = TRUE)

for (id in sample_id) {
  file_list <- list.files(id, full.names = TRUE)
  id <- basename(id)
  for(file in file_list) {
    new_file <- gsub(paste0(id, "_"), "", file)
    file.rename(from = file,
                to = new_file)
  }
}



# SampleSheet.tsv and Manifest.txt generation


data_dir <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/zhou_et_al_human"


manifest <- data.frame(key = ids::proquint(n = length(dir(file.path(data_dir, "expression_matrix"))), 
                                           n_words = 1L, use_cache = TRUE, use_openssl = FALSE),
                       filepath = dir(file.path(data_dir, "expression_matrix"), full.names = TRUE))


assay_metadata <- read.csv(file.path(data_dir, "snRNAseqAD_TREM2_assay_scRNAseq_metadata.csv"))

#biospeciment_metadata <- read.csv("snRNAseqAD_TREM2_biospecimen_metadata.csv")
#metadata <- merge(biospeciment_metadata, assay_metadata, by = "specimenID")
#rosmap_metadata <- read.csv("ROSMAP_Clinical_2019-05_v3.csv")
#metadata <- merge(metadata, rosmap_metadata, by = "individualID")


metadata <- readxl::read_xlsx(file.path(data_dir,"Copy of 41591_2019_695_MOESM2_ESM.xlsx"), sheet = "Rush samples", skip = 4)

metadata <- merge(metadata, assay_metadata, by.x = "Sample ID in snRNA-seq", by.y = "specimenID")


id_map <- match(purrr::map_chr(as.character(manifest$filepath), 
                               ~ strsplit(., "/", fixed = TRUE)[[1]][11]), 
                metadata$`Sample ID in snRNA-seq`)

metadata <- metadata[id_map, ]


SampleSheet <- data.frame(manifest = manifest$key,
                          sample_id = metadata$`Sample ID in snRNA-seq`,
                          individual = metadata$`Sample Identifier`,
                          pathological_diagnosis = ifelse(as.numeric(as.character(metadata$braaksc)) < 3, 
                                                          "Control", "AD"),
                          braak = as.numeric(as.character(metadata$braaksc)),
                          cerad = metadata$ceradsc,
                          amyloid = metadata$amyloid,
                          tangles = metadata$tangles,
                          cogdx = metadata$`cogdx / CDR`, #clinical diagnosis from ROSMAP
                          apoe = metadata$`APOE Genotype`,
                          trem2 = metadata$`TREM2 Genotype`,
                          age = round(metadata$Age),
                          sex = c("F", "M")[as.factor(metadata$Sex)],
                          PMI = round(metadata$pmi, 1),
                          brain_region = "PFC_BA10",
                          capdate = metadata$rnaBatch,
                          prepdate = metadata$libraryBatch,
                          seqdate = metadata$sequencingBatch)

dir.create(file.path(data_dir, "refs"))

write.table(manifest, file = file.path(data_dir, "refs/Manifest.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(SampleSheet, file = file.path(data_dir, "refs/SampleSheet.tsv"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# save.image("zhou_human_data_processing.RData")
