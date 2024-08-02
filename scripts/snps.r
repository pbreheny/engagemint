#' Ryan sent a list of SNPs; are any of them eQTLs in our tissues of interest?
#' First run scripts/rs-to-gc.r 
devtools::load_all(quiet=TRUE)

# Read in snps
snp <- fread(data_loc + 'snps.txt')
snp$id <- glue_data(snp, 'chr{V2}_{V3}')

# Loop over eQTLs
f <- list.files('~/lss/Fisher/gtex/v8/qtl/GTEx_Analysis_v8_eQTL/', '*.gene_pairs.txt.gz',)
tissue <- str_split_i(f, '\\.', 1)
dat_list <- vector('list', length(f))
pb <- txtProgressBar(0, length(f), style=3)
for (i in seq_along(f)) {
  # Read in eQTL
  eqtl <- fread(str_glue('~/lss/Fisher/gtex/v8/qtl/GTEx_Analysis_v8_eQTL/{f[i]}'))
  eqtl[, id := str_replace(variant_id, '^([^_]+)_([^_]+).*', '\\1_\\2')]
  
  # Merge
  dat_list[[i]] <- merge(snp, eqtl, by='id') %>%
    .[, .(rsid=V1, variant_id, tissue=tissue[i], gene_id = str_split_i(gene_id, '\\.', 1), tss_distance, slope=slope, pval=pval_nominal)]
  setTxtProgressBar(pb, i)
}
close(pb)
dat <- rbindlist(dat_list)

library("org.Hs.eg.db") %>% suppressPackageStartupMessages()
orgdb <- org.Hs.eg.db
# Annotate genes
anno <- select(orgdb, keys=unique(dat$gene_id), columns=c("SYMBOL", "GENENAME"), keytype="ENSEMBL") %>%
  data.table() %>%
  .[!duplicated(ENSEMBL)] %>%
  .[, .(gene_id=ENSEMBL, symbol=SYMBOL, name=GENENAME)]
out <- merge(dat, anno, by='gene_id') %>%
  .[, .(rsid, variant_id, tissue, gene_id, tss_distance, symbol, name, slope, pval)]
setkey(out, 'variant_id')
fwrite(out, '~/tmp/snp-eqtl.txt', sep='\t')

# Upload
system('cloud-share ~/tmp/snp-eqtl.txt')
