## Script to extend the target gene analysis to the macrophage response eQTLs from Alasoo 2018 and the microglia ABC from Kosoy 2022.
## Generates the data shown in Figures S6d-e

## Load libraries
library(data.table)
library(dplyr)
library(GenomicRanges)
library(readxl)
library(mariner)
library(liftOver)

### MACROPHAGE eQTL OVERLAP ###

## Read in the Alasoo QTL data
load("/work/users/m/a/marielle/qtl/alasoo/alasoo_filtered.rds")
## Read in the MPRA data
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T, sep = "\t")
## Active variants
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
## Add activity: 
mpra$activity <- ""
mpra[mpra$id %in% active,]$activity <- TRUE
mpra[mpra$activity=="",]$activity <- FALSE
emvar <- mpra[mpra$sig==TRUE & mpra$activity==TRUE,]$rsid

## For each of the alasoo datasets, figure out which have an eQTL (and what gene it connects to) for our emvars: 
overlapQTL <- list()
for (i in 1:length(qtlData)){
  sub <- qtlData[[i]]
  overlapQTL[[i]] <- as.data.frame(sub[sub$anchor1.rsid %in% emvar,])
}

overlap <- do.call(rbind, overlapQTL) #this is a subset for signifcant alasoo eQTLs for emvars only
table(overlap$dataset)

## add a category for stimulated vs just naive 
overlap$class <- "proinflammatory"
overlap[overlap$dataset=="naive",]$class <- "resting"

## how many unique variants: 
unique(overlap$anchor1.rsid) ##18
## how many unique genes: 
unique(overlap$anchor2.hgnc) ## 18

## which rsids are connected to each group? 
restingVar <- overlap[overlap$class=="resting",]$anchor1.rsid |> unique()
inflamVar <- overlap[overlap$class=="proinflammatory",]$anchor1.rsid |> unique()
ggvenn(data = list(resting = restingVar, inflammatory = inflamVar), auto_scale = TRUE) ## 3 unique to resting, 8 shared, 7 unique to inflammatory

## get the context specific variants
restOnly <- setdiff(restingVar, inflamVar)
inflamOnly <- setdiff(inflamVar, restingVar)
contextSpecific <- c(restOnly, inflamOnly)

## emvars that are significant interaction terms
interaction <- fread("/work/users/m/a/marielle/work/AD3D/tables/MPRA_allelic_interaction.txt") |> 
  filter(conditionalInteraction==TRUE) |> 
  pull(rsid)
contextSpecific[contextSpecific %in% interaction] ## 7 of these emVars that overlap a context-specific QTL show an interaction effect from MPRA 

## eGenes: 
egenes <- overlap$anchor2.hgnc |> unique()

### MICROGLIA ABC OVERLAP ###
## Read in the ABC pairs from Kosoy 2022
micABC <- fread("/work/users/m/a/marielle/external/kosoy2022/ABC_results_valid.txt")
## Rearrange the columns so we can easily convert to GInteractions: 
micABC <- 
  micABC[,c(11:13,1:3,4:9,14:24)] |> 
  as_ginteractions()

## Convet variants to GRanges object
var <- 
  GRanges(seqnames = Rle(mpra$seqnames), 
          ranges = IRanges(start = mpra$pos, end = mpra$pos), 
          rsid = mpra$rsid, 
          active = mpra$activity, 
          allelic = mpra$sig)

## Lift over to hg38 
ch <- "/work/users/m/a/marielle/ref/liftOver/hg19ToHg38.over.chain"
var38 <- liftOver(var, chain = import.chain(ch)) |> 
  unlist()
emvar <- var38[var38$active==TRUE & var38$allelic==TRUE]

## How many of the emvars overlap a microglia ABC pair? 
subsetByOverlaps(emvar, micABC) ##6
subsetByOverlaps(micABC, emvar) ## 19 pairs 
subsetByOverlaps(micABC, emvar)$TargetGene_name |> unique() ## 19 genes
abcGenes <- subsetByOverlaps(micABC, emvar)$TargetGene_name |> unique()
abcGenes <- abcGenes[!abcGenes==""]

