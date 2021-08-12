# Preparation of Mathys_et_al_2019 SampleSheet.tsv and Manifest.txt
# Nurun 16/02/2021

# Supplementary files were downloaded from https://www.nature.com/articles/s41586-019-1195-2#Abs1

basePath <- "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/tsai_et_al_alzheimers"

id_mapping <- read.csv(file.path(basePath, "SampleMetadata/snRNAseqPFC_BA10_id_mapping.csv"))
id_mapping <- id_mapping[, c(2:3)]
id_mapping <- id_mapping[!duplicated(id_mapping$projid), ]
id_mapping <- id_mapping[order(id_mapping$projid), ]
rownames(id_mapping) <- NULL


supplimentary3 <- read_xlsx(file.path(basePath, "SampleMetadata/41586_2019_1195_MOESM3_ESM.xlsx"), sheet = 2)
supplimentary5 <- read_xlsx(file.path(basePath, "SampleMetadata/41586_2019_1195_MOESM5_ESM.xlsx"))

clinical_info <- read.csv(file.path(basePath, "SampleMetadata/ROSMAP_Clinical_2019-05_v3.csv"))

clinical_info <- clinical_info[ , c("projid" ,"apoe_genotype", "pmi")]

sampleinfo <- merge(id_mapping, supplimentary3, by = "Subject")

idx <- which(colnames(sampleinfo) %in% colnames(supplimentary5))
idx <- idx[-1]

sampleinfo <- sampleinfo[-idx]
sampleinfo <- merge(sampleinfo, supplimentary5, by = "Subject")
sampleinfo <- merge(sampleinfo, clinical_info, by = "projid")

sampleinfo$age_death <- gsub("+", "", sampleinfo$age_death, fixed = TRUE)
sampleinfo$age_death <- as.integer(sampleinfo$age_death)

sampleinfo <- sampleinfo[order(sampleinfo$projid), ]

data_dir <- file.path(basePath, "ExpressionData/projid")

manifest <- data.frame(key = ids::proquint(n = length(dir(data_dir)), 
                                           n_words = 1L, use_cache = TRUE, use_openssl = FALSE),
                       filepath = dir(data_dir, full.names = TRUE))

all(sampleinfo$projid == purrr::map_chr(as.character(manifest$filepath), 
                                                     ~ strsplit(., "/", fixed = TRUE)[[1]][12]))


SampleSheet <- data.frame(manifest = manifest$key, 
                          sample_id = sampleinfo$projid,
                          individual = sampleinfo$Subject,
                          pathological_diagnosis = ifelse(sampleinfo$braaksc < 3, 
                                                          "Control", "AD"),
                          braak = sampleinfo$braaksc,
                          cerad = sampleinfo$ceradsc,
                          amyloid = sampleinfo$amyloid,
                          tangles = sampleinfo$tangles,
                          cogdx = sampleinfo$cogdx,
                          pathological_diagnosis_original = c("Control", "AD")[ as.factor(sampleinfo$`pathologic diagnosis of AD`)], #pathological diagnosis used in the original study
                          group = c("Early-pathology", "Late-pathology", "No-pathology")[as.factor(sampleinfo$pathology.group)],
                          apoe = sampleinfo$apoe_genotype,
                          age = sampleinfo$age_death, 
                          sex = c("F", "M")[as.factor(sampleinfo$msex)],
                          PMI = round(sampleinfo$pmi, 1),
                          brain_region = "PFC_BA10")


write.table(manifest, file = file.path(basePath, "refs/Manifest.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(SampleSheet, file = file.path(basePath, "refs/SampleSheet.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
