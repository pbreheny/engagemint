#' Read in eCLIP data
devtools::load_all(quiet=TRUE)

eclipfile <- read.table(data_loc + "format_All_targetsite_final.bed", header = TRUE)

process_ec <- function(x) {
  mutate(x, Tissue = as.factor(Tissue)) %>%
    mutate(strand = as.factor(strand)) %>%
    mutate(`Min -Log10 P Value` = ifelse(`MinPval.Log10.` > 100, 
                                         100.00, `MinPval.Log10.`)) %>%
    relocate(`Min -Log10 P Value`, .after = Cluster_Coord) %>%
    select(-`MinPval.Log10.`)
}

if (!dir.exists('app/rds')) dir.create('app/rds')
process_ec(eclipfile) %>%
  saveRDS('app/rds/eclip.rds')
