#' Read in eCLIP data
devtools::load_all(quiet=TRUE)

bedpaths <- list.files(data_loc + "bed12", full.names = TRUE, pattern = '*.bed')
bedfiles <- list()
for(i in 1:length(bedpaths)) bedfiles[[i]] <- fread(bedpaths[i])
names(bedfiles) <- gsub(".*\\/", "", bedpaths)

# Process eclip
process_ec <- function(x) {
  select(x, c(Chromosome = V1, Strand = V6, `Gene Symbol` = V4, 
              `-Log10 P Value` = V5, `Window Start` = V2, `Window End` = V3)) %>%
    mutate(`-Log10 P Value` = ifelse(`-Log10 P Value` > 100, 
                                     100.00, round(`-Log10 P Value`, 2)))
}
lapply(bedfiles, process_ec) %>%
  saveRDS('rds/eclip.rds')
