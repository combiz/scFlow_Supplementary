# Preparation of Mathys_et_al_2019 data for scflow
# Combiz 21/11/2019

# Dataset was downloaded from https://www.synapse.org/#!Synapse:syn18485175
# note that the cellranger output was first compressed with gzip and renamed
# to the standard features.tsv.gz, barcodes.tsv.gz, and matrix.mtx.gz names

library(scFlow)
library(purrr)
library(dplyr)

# the cellranger output
mat <- read_sparse_matrix("/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/tsai_et_al_alzheimers/ExpressionData/CellRangerOutput")

# we load this to obtain the numeric id to project id mappings
notfiltered_column_metadata <- read.delim(
  "/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/tsai_et_al_alzheimers/ExpressionData/notfiltered/notfiltered_column_metadata.txt")

# the id (numeric 1-48 is affixed to the barcodes)
ids <- purrr::map_chr(
  as.character(notfiltered_column_metadata$TAG),
  ~ strsplit(., "\\.")[[1]][2]
  )

# we obtain a mapping of the numeric 1-48 id to projid
mappings <- unique(data.frame(
  id = as.factor(ids),
  projid = as.factor(notfiltered_column_metadata$projid))
  )

# the ids for the large matrix
mat_ids <- purrr::map_chr(
  as.character(colnames(mat)),
  ~ strsplit(., "\\-")[[1]][2]
)

basepath <- file.path("/rds/general/project/ukdrmultiomicsproject/live/public_datasets/single_nuclei/tsai_et_al_alzheimers/ExpressionData/projid")
for(id in mappings$id) {
  projid <- as.character(mappings[mappings$id == id, ]$projid)
  print(sprintf("Project ID: %s (%s)", projid, id))
  outdir <- file.path(basepath, projid)
  dir.create(outdir, showWarnings = FALSE)
  write_sparse_matrix(
    mat[, mat_ids == id],
    outdir
  )
}

