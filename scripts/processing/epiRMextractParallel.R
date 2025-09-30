## Script to extract DNase-seq scores from epiRM cell types from all MPRA-elements 
## Note I wrote this script to be run within a series of sbatch scripts so it can be run in parallel, it takes ~10-20 min to run for each cell type in my experience 

## Load libraries 
library(plotgardener)
library(dplyr)
library(ggplot2)

## Read in data: 
## MPRA activity at variant level
load("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9_mpraset.rda")
## For RNA and DNA, get the median counts across both alleles and conditions
rnamat$median <- apply(rnamat, 1, median)
dnamat$median <- apply(dnamat, 1, median)
mpraCounts <- data.frame(id = rownames(rnamat), 
                         rna = rnamat$median, 
                         dna = dnamat$median)
mpraCounts$mpraRatio <- mpraCounts$rna/mpraCounts$dna
mpraCounts$logRatio <- log(mpraCounts$mpraRatio)
## Get column for coords: 
mpraCounts$chr <- unlist(strsplit(mpraCounts$id, "_"))[seq(1,(nrow(mpraCounts)*2),by=2)]
mpraCounts$pos <- unlist(strsplit(mpraCounts$id, "_"))[seq(2,(nrow(mpraCounts)*2),by=2)] |> as.integer()

## Remove NAs
mpraCounts <- na.omit(mpraCounts)

## Add  quantile information
mpraCounts <- mpraCounts[order(mpraCounts$logRatio, decreasing = TRUE),]
n <- 10
mpraCounts <- 
  mpraCounts %>%
  mutate(quantile = ntile(mpraCounts$logRatio, n))
mpraCounts$quantile <- factor(mpraCounts$quantile, levels = c(1:n))
## Scale to center at 0
mpraCounts$logRatio <- scale(mpraCounts$logRatio, center = TRUE, scale = FALSE)
mpraCounts$atac <- ""
#parameters
buffer <- 200

## Pass SLURM arguments into script: 
args <- commandArgs(trailingOnly = TRUE)
atacBW <- args[1]
name <- unlist(strsplit(atacBW, "/work/users/m/a/marielle/external/epiRM/atac/signal/"))[2]
name <- unlist(strsplit(name, "-"))[1]

## Loop through each variant
for (i in 1:nrow(mpraCounts)){
  print(i)
  atac <- readBigwig(file = atacBW,
                     chrom = mpraCounts$chr[i],
                     chromstart = mpraCounts$pos[i]-buffer,
                     chromend = mpraCounts$pos[i]+buffer)
  atacScore <- max(atac$score)
  if (length(atacScore)==0){
    atacScore <- 0
  }
  mpraCounts$atac[i] <- as.numeric(atacScore)
  
}
mpraCounts$celltype <- name
mpraCounts$atac <- as.numeric(mpraCounts$atac)
write.table(mpraCounts, file = paste0("/work/users/m/a/marielle/work/AD3D/data/epiRM/atac", name, "_DNase_data.txt"), row.names = FALSE, sep = "\t", quote = F)
