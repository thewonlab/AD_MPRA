## Script to run motifBreakR on all MPRA variants 
## This is to be run via various sbatch jobs to be run in parallel 
## Each variant is passed in as an argument for this script 

## Motif BreakR in parallel 
library(tidyverse)
library(magrittr)
library(data.table)
library(motifbreakR)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(nullranges)

## Get motif data
motif_list <- MotifDb[which(mcols(MotifDb)$organism == "Hsapiens")]
## remove stamlab Motifs, because they are not annotated
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "stamlab")]
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "cisbp_1.02")]
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "hPDI")]
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "jolma2013")]
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "SwissRegulon")]
motif_list <- motif_list[which(mcols(motif_list)$dataSource != "UniPROBE")]

## From the N9N6 pooling: 
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

mpraGR <- 
  GRanges(adjusted$seqnames,
          IRanges(adjusted$pos, adjusted$pos),
          SNP_id = adjusted$rsid,
          id = adjusted$id, 
          sig = adjusted$sig)
mpraGR$active <- FALSE
mpraGR[mpraGR$id %in% active,]$active <- TRUE
mpraGR$putative <- FALSE
mpraGR[mpraGR$sig==TRUE & mpraGR$active==TRUE]$putative <- TRUE

## Set whether we are using the matched set or the putative causal variants: 
## Emvars (aka putative)
## Note, I commented out which variants were not being actively run 
vars <- mpraGR[mpraGR$putative==TRUE]$SNP_id #putative
#vars <- mpraGR[mpraGR$putative==FALSE]$SNP_id # background all not putative 

adjusted <- adjusted[adjusted$rsid %in% vars,]

## Convert into the appropriate format for motifBreakR
vars <- 
  GRanges(adjusted$seqnames,
          IRanges(adjusted$pos, adjusted$pos),
          SNP_id = adjusted$rsid,
          REF = DNAStringSet(adjusted$a2),
          ALT = DNAStringSet(adjusted$a1))

## Format
seqlevels(vars) <- paste0("chr", c(1:22))
names(vars) <- vars$SNP_id
vars <- sort(vars)
seqlengths(vars) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:22]
isCircular(vars) <- isCircular(BSgenome.Hsapiens.UCSC.hg19)[1:22]
genome(vars) <- "hg19"
strand(vars) <- "*"
attributes(vars)$genome.package <- "BSgenome.Hsapiens.UCSC.hg19"


args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args)

## For certain variants, swap the ref and the alt allele to match the database
## For all variants being tested, figure out if it matches or if it doesn't. If it doesn't match, then swap
pos <- GPos(seqnames = seqnames(vars),
            pos = start(ranges(vars)), 
            RefSNP_id = vars$SNP_id, 
            alleles_as_ambig = as.character(vars$REF))

## Determine what the actual ref and alt alleles for this variant are 
actual <- inferRefAndAltAlleles(pos[i], "hg19")
if(!actual$ref_allele==vars$REF[i]){
  ref <- vars[i]$REF
  alt <- vars[i]$ALT
  vars[i]$REF <- alt
  vars[i]$ALT <- ref
}


## Run
print("Starting Run")
results <- motifbreakR(
  snpList = vars[i], filterp = TRUE,
  pwmList = motif_list,
  method = "ic",
  threshold = 1e-3,
  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
)
print("calculating p-value")
results <- calculatePvalue(results, granularity = 1e-5)
results <- results[which(mcols(results)$effect == "strong")]
results <- data.frame(results)
results <- results[order(results$alleleDiff),]
results$geneSymbol <- factor(results$geneSymbol, levels = unique(results$geneSymbol))
results$motifPos <- sapply(results$motifPos, paste, collapse = ",")

## Save the output to a dataframe
print("writing output")
write.table(as.data.frame(results), file = paste0("/work/users/m/a/marielle/work/AD3D/data/motifbreakr/allBackground/", i, "_matched_results.txt"), quote = F, sep = "\t", row.names = FALSE)
