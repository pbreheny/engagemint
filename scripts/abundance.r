#' Read in miRNA abundance
devtools::load_all(quiet=TRUE)
library(vroom)
mir_count <- vroom(
  data_loc + "miR quant/master_miRcount.csv",
  col_select = 1:26) %>%
  na.omit() %>%
  select(-c(Species_ID,
            grep("fam_sum", colnames(.)),
            grep("geomean", colnames(.)))) %>%
  set_colnames(c(
    "miR Family 6", "miR Family 8", "Seed+m8",
    "MiRBase ID", "Mature Sequence Family",
    "Conservation", "MiRBase Accession", "Left Ventricle<br>Rank",
    "Aorta<br>Rank", "Gastrointestinal<br>Rank", "Hippocampus<br>Rank",
    "Cerebellum<br>Rank", "Frontal Cortex<br>Rank"))
mir_count[,c("MiRBase ID", "Aorta<br>Rank", "Cerebellum<br>Rank", 
             "Frontal Cortex<br>Rank", "Gastrointestinal<br>Rank",
             "Hippocampus<br>Rank", "Left Ventricle<br>Rank",
             "Conservation", "Seed+m8", "miR Family 6", "miR Family 8",
             "Mature Sequence Family")] %>%
  saveRDS('rds/abundance.rds')
