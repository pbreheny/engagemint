#' Generates images for each gene
#' This is "prototype" code, as the real tracks should include multiple samples

library(data.table)
library(magrittr)
library(GenomicRanges)
library(readxl)
library(Homo.sapiens) # human annotation package

# read in summarized HITS-CLIP data
path <- "R:/rdss_pbreheny/grade/CLIP_bg/miR-CLIP database build example table_Boudreau lab.xlsx"
miRdata <- lapply(excel_sheets(path), read_excel,
                  path = path, na = "NA")
names(miRdata) <- excel_sheets(path)

# retrieve list of HITS-CLIP genes
hits_genes <- lapply(miRdata, "[", "Gene Symbol") %>%
  unlist(use.names = FALSE) %>%
  na.omit() %>%
  unique()

# match genes to respective GRanges
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg19.knownGene
allcoords <- transcripts(Homo.sapiens, columns = "SYMBOL")
hits_coords <- subset(allcoords, SYMBOL %in% hits_genes)

# subset to standard chromosomes - will revisit more options later, this was for the simplest integration
hits_coords <- keepSeqlevels(hits_coords, standardChromosomes(hits_coords)[1:24], 
                             pruning.mode = "coarse")

# create union of all gene coordinates from same gene
genenames <- unique(unlist(hits_coords$SYMBOL))
genedat <- GRanges()
for(i in 1:length(genenames)) {
  genedat[i] <- range(subset(hits_coords, SYMBOL %in% genenames[i]))
} 
genedat$SYMBOL <- genenames

# read in raw peak data - this is for one individual I believe
bgfile <- fread("R:/rdss_pbreheny/grade/CLIP_bg/bedgraphs/4_1715_mapped_hg19_PLUS.bg", skip = 1)

peakranges <- GRanges(seqnames = bgfile$V1,
                  ranges = IRanges(start = bgfile$V2, end = bgfile$V3),
                  count = bgfile$V5)

# generate images for each gene, save as png
for(i in 1:length(genenames)) {
  
  temp <- subsetByOverlaps(peakranges, genedat[genedat$SYMBOL == genenames[i]])
  
  if(!(length(temp$count) == 0)) {
    
    chr <- temp@seqnames@values
    cov <- na.omit(coverage(temp, weight = "count")[chr][[1]])
    
    png(filename = paste0("Images/", genenames[i], ".png"), width = 1150, height = 300)
    plot(start(cov)[2]:length(cov), cov[start(cov)[2]:length(cov)], type = "s",
         xlab = "Coordinate", ylab = "Read Coverage", las = 1, frame = FALSE)
    box(bty = "l")
    dev.off()
    
  }
}

#### image for presentation - gene QKI ####

# chr6 <- bgfile[bgfile$V1 == "chr6"]
# 
# rangeschr6 <- GRanges(seqnames = chr6$V1,
#                       ranges = IRanges(start = chr6$V2, end = chr6$V3),
#                       count = chr6$V5)
# 
# cov <- coverage(rangeschr6[133778:133915], weight = "count")[[1]]
# 
# plot(start(cov)[2]:length(cov), cov[start(cov)[2]:length(cov)], type = "s", 
#      xlab = "Chr6 Coordinate", ylab = "Read Coverage", las = 1,
#      main = "QKI HITS-CLIP Peak")
# 
# qki <- subsetByOverlaps(rangeschr6, genedat[genedat$SYMBOL == "QKI"])
# cov <- na.omit(coverage(qki, weight = "count")[[1]])
# 
# plot(start(cov)[2]:length(cov), cov[start(cov)[2]:length(cov)], type = "s", 
#      xlab = "Chr6 Coordinate", ylab = "Peak Height", las = 1,
#      main = "QKI HITS-CLIP")
