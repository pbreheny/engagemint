#' Read in original HITS-CLIP peak data
devtools::load_all(quiet=TRUE)
suppressPackageStartupMessages(library(dplyr))
library(readxl)
path <- data_loc + "miR-CLIP database build example table_Boudreau lab.xlsx"
mir_data <- lapply(excel_sheets(path), read_excel, 
                  path = path, na = "NA")
names(mir_data) <- excel_sheets(path)

# Process each hits-clip sample
process_hc <- function(x) {
  x %>%
    mutate(`P-value` = ifelse(`P-value` == 0 | `P-value` < 1e-100,
                              1e-100,
                              `P-value`)) %>%
    mutate(`Avg Peak Height` = round(`Avg Peak Height`, 2),
           `-Log10 P Value` = round(-log10(`P-value`), 2)) %>%
    relocate(`-Log10 P Value`, .after = `Avg Peak Height`) %>%
    select(-`P-value`) %>%
    mutate(`Genomic Coordinates Link (Hg19)` =
             hyperlink(`UCSC Link (can just be a button)`,
                       `Genomic Coordinates (Hg19)`)) %>%
    relocate(`Genomic Coordinates Link (Hg19)`) %>%
    select(-`UCSC Link (can just be a button)`)
}

lapply(mir_data, process_hc) %>%
  saveRDS('rds/hits-clip.rds')
