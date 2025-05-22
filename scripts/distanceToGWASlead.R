## How far are emVars from the lead SNP? What about how far the lead SNPs are from all other variants, not just emVars? 
## Generates figure s1d

## Load libraries
library(data.table)
library(readxl)
library(GenomicRanges)

## GWAS: 
gwas <- fread("/work/users/m/a/marielle/gwas/jansen/AD_sumstats_Jansenetal_2019sept_filtered_0.05.txt.gz")
gwas <- gwas[,c(2,3,8,6)]
colnames(gwas) <- c("chrom", "pos", "p", "snp")
gwas$chrom <- paste0("chr", gwas$chrom)

## MPRA: 
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
## Add on risk and protective columns, and if corrected lfc = negative, swap a1 and a2
## In the end, a2 should always be alt and a1 should always be ref, where a1=ref="protective" and a2=alt="risk"
#Loop to adjust when corrected LFC < 0 but keep if corrected LFC > 0
for (i in 1:nrow(adjusted)){
  if (adjusted$corrected_logFC[i]>0){
    adjusted$risk[i] <- adjusted$a2[i]
    adjusted$protective[i] <- adjusted$a1[i]
  } else if (adjusted$corrected_logFC[i]<0){
    adjusted$risk[i] <- adjusted$a1[i] #SWAP 
    adjusted$protective[i] <- adjusted$a2[i]
  }
}

allelic <- adjusted[adjusted$sig==TRUE,]

### MPRA active
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
activeRSID <- adjusted[adjusted$id %in% active,]$rsid

## Add activity: 
adjusted$activity <- ""
adjusted[adjusted$id %in% active,]$activity <- TRUE
adjusted[adjusted$activity=="",]$activity <- FALSE
adjusted <- adjusted[!adjusted$seqnames=="chr3",]

## Select putative causal variants to look at 
variants <- adjusted

## Overlap with GWAS loci 
loci <- read_xlsx("/work/users/m/a/marielle/work/MPRA/align/jansenLonelySNPs.xlsx") |> 
  as.data.frame()
loci <- loci[!is.na(loci$pos),] #this is the 29 loci from Jansen et al 2019
regular <- loci[!is.na(loci$`Start bp`),]
lonely <- loci[is.na(loci$`Start bp`),]
regular <- GRanges(seqnames = Rle(regular$Chr), 
                   ranges = IRanges(start = regular$`Start bp`, end = regular$`End bp`), 
                   rsid = regular$LeadSNPs)
lonely <- GRanges(seqnames = Rle(lonely$Chr), 
                  ranges = IRanges(start = lonely$pos, end = lonely$pos+1), 
                  rsid = lonely$LeadSNPs)
#combine
loci <- c(regular, lonely)
seqlevelsStyle(loci) <- "UCSC"
loci <- sort(loci)

## Convert to GRanges: 
buffer <- 0.2E6
putative <- GRanges(seqnames = Rle(variants$seqnames), 
                    ranges = IRanges(start = variants$pos-buffer, end = variants$pos+buffer), 
                    id = variants$id, 
                    rsid = variants$rsid, 
                    pos = variants$pos, 
                    active = variants$activity, 
                    allelic = variants$sig)
levels <- paste0("chr", 1:22)
levels <- levels[levels %in% seqnames(putative)]
seqlevels(putative) <- levels
putative <- sort(putative)

## Determine which GWAS locus each of the variants belongs to 
putative$locus <- nearest(putative, loci)

## Figure out the classes for each of these variants: 
emvar <- adjusted[adjusted$sig==TRUE & adjusted$activity==TRUE,]$rsid
allelic <- adjusted[adjusted$sig==TRUE & !adjusted$rsid %in% emvar,]$rsid
active <- adjusted[adjusted$activity==TRUE & !adjusted$rsid %in% emvar,]$rsid
inactive <- adjusted[!adjusted$rsid %in% c(emvar, allelic, active),]$rsid

putative$class <- ""
putative[putative$rsid %in% emvar,]$class <- "emvar"
putative[putative$rsid %in% allelic,]$class <- "allelic"
putative[putative$rsid %in% active,]$class <- "active"
putative[putative$rsid %in% inactive,]$class <- "inactive"

## For each locus, determine how far the variants are 
distances <- list()
for (i in 1:length(loci)){
  leadPos <- start(loci[i])
  sub <- putative[putative$locus==i]
  sub$dist <- abs(sub$pos - leadPos)
  distances[[i]] <- sub
}
data <- do.call(`c`, distances) |> 
  as.data.frame()
data$dist <- data$dist/100e5

distanceDensity <- 
  ggplot(data, aes(x = dist, fill = class)) + 
  geom_density(alpha = 0.7, color = NA) + 
  xlim(0,0.5E6) + 
  scale_fill_manual(values = c("#9FD141", "#66ABF2", "#158C7E", "#D3DAE0")) + 
  theme_minimal() + 
  xlab("Distance between variant\nand GWAS lead SNP (100kb)") + ylab("Density") + 
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8))
save(distanceDensity, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/distanceDensity_lmer.rda")


