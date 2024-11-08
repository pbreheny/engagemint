#' Ryan sent a list of SNPs with rs IDs; this scripts converts them to
#' genomic coordinates for downstream analysis
#' e.g.: rs12345 -> chr1_1012429_G_C_b38
devtools::load_all(quiet=TRUE)
library(GenomicRanges) |> suppressPackageStartupMessages()
library(SNPlocs.Hsapiens.dbSNP155.GRCh38) |> suppressPackageStartupMessages()
library(XtraSNPlocs.Hsapiens.dbSNP144.GRCh38) |> suppressPackageStartupMessages()

# Read SNPs
rsid <- readxl::read_excel(data_loc + 'New eCLIP SNPs Oct 2024.xlsx', col_names = 'v1')$v1

# Map SNPs
gr_snp <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, rsid, ifnotfound='drop')
out_snp <- data.table(
  gr_snp$RefSNP_id,
  seqnames(gr_snp) |> as.character(),
  start(gr_snp),
  IUPAC_CODE_MAP[gr_snp$alleles_as_ambig])

# Map Indels
gr_ind <- snpsById(XtraSNPlocs.Hsapiens.dbSNP144.GRCh38, rsid, ifnotfound='drop')
out_ind <- data.table(
  gr_ind$RefSNP_id,
  seqnames(gr_ind) %>% as.character() %>% str_replace('ch', ''),
  start(gr_ind),
  NA)

# Merge and save
out <- rbind(out_snp, out_ind)
fwrite(out, glue('{data_loc}/new-snps.txt'), sep='\t', col.names=FALSE)

# Missing SNPs
miss <- setdiff(rsid, out$V1) 
pct <- 100*length(miss)/nrow(out)
glue('{length(miss)} missing / {nrow(out)} total ({breheny::.f(pct, 1)} %) ')
head(miss)
out
