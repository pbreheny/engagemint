#' Read in original HITS-CLIP peak data
devtools::load_all(quiet=TRUE)
suppressPackageStartupMessages(library(dplyr))
library(readxl)
path <- data_loc + "miR-CLIP database build example table_Boudreau lab.xlsx"
hitsfile <- lapply(excel_sheets(path), read_excel, 
                  path = path, na = "NA")
names(hitsfile) <- excel_sheets(path)

# Process each hits-clip sample
process_hc <- function(x) {
  temp <- x
  for(i in names(temp)) {
    temp[[i]]$Tissue <- as.factor(i)
    
    temp[[i]] <- temp[[i]] %>%
      mutate(`P-value` = ifelse(`P-value` == 0 | `P-value` < 1e-100,
                                1e-100,
                                `P-value`)) %>%
      mutate(`Avg Peak Height` = round(`Avg Peak Height`, 2),
             `-Log10 P Value` = round(-log10(`P-value`), 2)) %>%
      mutate(`Genomic Coordinates Link (Hg19)` =
               hyperlink(`UCSC Link (can just be a button)`,
                         `Genomic Coordinates (Hg19)`)) %>%
      relocate(`-Log10 P Value`, .after = `Avg Peak Height`) %>%
      relocate(`Genomic Coordinates Link (Hg19)`) %>%
      relocate(Tissue) %>%
      select(-`P-value`) %>%
      select(-`UCSC Link (can just be a button)`)
  }
  
  do.call(rbind, temp)
}

if (!dir.exists('app/rds')) dir.create('app/rds')
process_hc(hitsfile) %>%
  saveRDS('app/rds/hits-clip.rds')
