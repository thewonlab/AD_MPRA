## Script to find ane example region to show in Figure 4 that fits the following criteria: 
## (1) Significant eQTL
## (2) Significant ABC pair
## (3) gene associated with AD fr5om target gene enrichment: 

## Load libraries
library(plotgardener)
library(data.table)
library(readxl)
library(grid)
library(ggplot2)
library(tidyr)
library(ggh4x)
library(cowplot)
library(mariner)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(RColorBrewer)
library(purrr)

## Read in mpra allelic and active elements
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
adjusted$active <- FALSE
adjusted[adjusted$id %in% active,]$active <- TRUE
emvar <- adjusted[adjusted$sig==TRUE & adjusted$active == TRUE,]$rsid

## Read in genomic data
hicFiles <- c("/proj/phanstiel_lab/Data/processed/LIMA/hic/LIMA_THP1_WT_LPIF_CMB_S_0.0.0_megaMap/aligned/LIMA_THP1_WT_LPIF_CMB_S_0.0.0_omegaMap_inter_withNorms.hic")
rnaBW <- "/work/users/m/a/marielle/proj/LIMA/signal/rna/LIMA_RNA_MERGE_0_1440.bw"
atacBW <- "/work/users/m/a/marielle/proj/LIMA/signal/atac/LIMA_ATAC_MERGE_0_1440.bw"
enhBW <- "/work/users/m/a/marielle/proj/LIMA/signal/cnr/LIMA_H3K27ac_MERGE_0_1440.bw"

## GWAS: 
gwas <- fread("/work/users/m/a/marielle/gwas/jansen/AD_sumstats_Jansenetal_2019sept_filtered_0.05.txt.gz")
gwas <- gwas[,c(2,3,8,6)]
colnames(gwas) <- c("chrom", "pos", "p", "snp")
gwas$chrom <- paste0("chr", gwas$chrom)

## reformat mpra data
mpra <- adjusted[,c(10,11,6,13,9)]
colnames(mpra) <- c("chrom", "pos", "p", "snp", "sig")

## QTL Data: 
qtl <- list.files("/work/users/m/a/marielle/work/AD3D/data/qtl/plot/", "eqtl", full.names = TRUE)
ch <- import.chain("/work/users/m/a/marielle/ref/liftOver/hg38ToHg19.over.chain")

## ABC
abc <- read.table("/work/users/m/a/marielle/work/AD3D/data/abc/LIMA_ABC_enhOverlap.txt", header=TRUE)

## subet for variants of interest
variants <- adjusted[adjusted$rsid %in% emvar,]

## Set parameters
colors <- brewer.pal(8, "YlGnBu")
tp <- c("0", "24")

## Convert to GRanges: 
buffer <- 0.01E6
putative <- GRanges(seqnames = Rle(variants$seqnames), 
                    ranges = IRanges(start = variants$pos-buffer, end = variants$pos+buffer), 
                    id = variants$id, 
                    putative = variants$rsid, 
                    pos = variants$pos)
levels <- paste0("chr", 1:22)
levels <- levels[levels %in% seqnames(putative)]
seqlevels(putative) <- levels
putative <- sort(putative)

## For figure 4 we look at two variants: 
#4f: rs72838287 (i=4)
#4g: rs12721109 (i=43)

## PLOT 
pageCreate(width = 8.75, height = 8, showGuides = FALSE, xgrid = 0.0, ygrid = 0.0)

## bin1
params <- 
  pgParams(chrom = "chr2",
           chromstart = 127795786,
           chromend = 127950786,
           assembly = "hg19",
           x = 4,
           width = 2.5,
           height = 0.5)

## apoe
params <- 
  pgParams(chrom = "chr19",
           chromstart = 45287221,
           chromend = 45607221,
           assembly = "hg19",
           x = 4,
           width = 2.5,
           height = 0.5)

pst <- params
pst$height <- 0.3
params$chrom <- unlist(strsplit(params$chrom, "chr"))[2]

hic_ypos <- c(1.725)
atac_ypos <- c((4.35))
h3k_ypos <- c((atac_ypos+pst$height+0.05))
rna_ypos <- c((h3k_ypos+pst$height+0.05))

## Plot HiC
plotHicRectangle(data = hicFiles, params = params, y = hic_ypos, zrange = c(0,1000), height = 0.9, res = 5000)

## GWAS
gwasMax <- gwas[gwas$chrom == pst$chrom & gwas$pos>pst$chromstart & gwas$pos< pst$chromend,]$p
gwasMax <- gwasMax[gwasMax>0]
gwasMax <- log10(gwasMax) * -1
gwasMax <- max(gwasMax) |> round()
#gwasMax <- 100 #toggle for APOE 

gwasMan <- plotManhattan(data = gwas, 
                         params = pst,
                         range = c(0,gwasMax+6),
                         y = 2.7, height = 0.5, 
                         sigLine = TRUE, col = "gray90", fill = "#B2BBBF",
                         baseline = TRUE, baseline.lwd = 2, baseline.color = "black",
                         snpHighlights = data.frame(snp = c(putative[i]$putative),
                                                    pch = c(23),
                                                    cex = c(0.5),
                                                    col = c("black"), 
                                                    fill = c("black")))
annoYaxis(
  plot = gwasMan,
  at = c(seq(0, gwasMax, by = gwasMax %/% 3), gwasMax+5),
  axisLine = TRUE, fontsize = 8
)

plotText(
  label = "-log10(p)", x = 3.6, y = 2.95, rot = 90,
  fontsize = 6, fontface = "bold", just = "center",
  default.units = "inches"
)

# ## MPRA
mpraMax <- mpra[mpra$chrom == pst$chrom & mpra$pos>pst$chromstart & mpra$pos< pst$chromend,]$p |> na.omit()
mpraMax <- log10(mpraMax) * -1
mpraMax <- max(mpraMax) |> ceiling()
mpraMan <- plotManhattan(mpra,
                         params = pst,
                         y = 3.25, height = 0.5,
                         sigLine = TRUE, col = "gray90", sigVal = 0.1, fill = "#B2BBBF",
                         #fill = colorby("ratio", colorRampPalette(c("gray90", "#4F52A3"))),
                         baseline = TRUE, baseline.lwd = 2, baseline.color = "black", range = c(0,mpraMax+0.5),
                         snpHighlights = data.frame(snp = c(putative[i]$putative),
                                                    pch = c(23),
                                                    cex = c(0.5),
                                                    col = c("black"),
                                                    fill = c("black")))

if(mpraMax==1){
  annoYaxis(
    plot = mpraMan,
    at = c(0,2),
    axisLine = TRUE, fontsize = 8
  )
} else {
  annoYaxis(
    plot = mpraMan,
    at = c(seq(0, mpraMax, by = mpraMax %/% 2), mpraMax),
    axisLine = TRUE, fontsize = 8
  )
}

plotText(
  label = "-log10(padj)", x = 3.6, y = 3.5, rot = 90,
  fontsize = 6, fontface = "bold", just = "center",
  default.units = "inches"
)

## QTL: 
## First just read in the specific chromosome: 
qtlDataFull <- read.table((qtl[grep(pst$chrom, qtl)])[1], header=TRUE)
## Get just the variants that are also in the GWAS data
qtlData <- qtlDataFull[qtlDataFull$Variant %in% gwas$snp,]
colnames(qtlData)[4] <- "snp"
qtlData <- merge(gwas, qtlData, by = "snp")
## if you know what gene you're looking at, do so here, if not, comment out
qtlData <- qtlData[qtlData$Gene=="ENSG00000136717",] #BIN1 - 2

## Get ymax
eqtlMax <- qtlData[qtlData$pos>params$chromstart & qtlData$pos< params$chromend,]$pvalue #tried using BH before, too strict? 
eqtlMax <- eqtlMax[eqtlMax>0]
eqtlMax <- log10(eqtlMax) * -1
eqtlMax <- max(eqtlMax) |> round()
#eqtlMax <- 75
#rearrange
qtlData <- qtlData[,c(2,3,15,1,17)] #16 is for BH, 15 for pval
colnames(qtlData) <- c("chrom", "pos", "p", "snp", "sig")
qtlMan <- plotManhattan(qtlData,  
                        params = pst, 
                        range = c(0,eqtlMax+1),
                        y = 3.8, height = 0.5,
                        sigLine = TRUE, sigVal = 0.05, col = "gray90", fill = "#B2BBBF",
                        baseline = TRUE, baseline.lwd = 2, baseline.color = "black", 
                        #fill = colorby("sig", colorRampPalette(c("gray80", "#2D2976"))),
                        snpHighlights = data.frame(snp = c(putative$putative[i]),
                                                   pch = c(23),
                                                   cex = c(0.5),
                                                   col = c("black"), 
                                                   fill = c("black")))

annoYaxis(
  plot = qtlMan,
  at = c(seq(0, eqtlMax, by = eqtlMax %/% 3), eqtlMax),
  axisLine = TRUE, fontsize = 8
)

plotText(
  label = "-log10(p)", x = 3.6, y = 4.05, rot = 90,
  fontsize = 6, fontface = "bold", just = "center",
  default.units = "inches"
)

## ABC Arches
abcSub <- abc[abc$anchor2.TargetGene=="APOE" & abc$rsid %in% putative$putative[i],]
abcSub <- abcSub[,c(1:3,6:8)]


## If just plotting ABC, use y = 3.8. If plotting both, use y = 4.3
archPlot <- plotPairsArches(
  data = abcSub, params = pst,
  fill = "gray10", alpha = 1, flip = TRUE,
  x = params$x, y = 3.8, height = 0.5, lwd = 0.2,
  just = c("left", "top"),
  default.units = "inches", 
)

## ATAC 
atacData <- lapply(atacBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
atacRange <- c(0, lapply(atacData, `[[`, 4) |> unlist() |> max())
#atacRange <- c(0, 100)
## Plot ATAC
atacS <- pmap(list(atacData, atac_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = atacRange, linecolor = colors[4], fill = colors[4]))

### H3K27ac
h3kData <- lapply(enhBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
h3kRange <- c(0, lapply(h3kData, `[[`, 4) |> unlist() |> max())
#h3kRange <- c(0, 100)
## Plot H3K27ac
h3k27S <- pmap(list(h3kData, h3k_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = h3kRange, linecolor = colors[6], fill = colors[6]))

## RNA
rnaData <- lapply(rnaBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
rnaRange <- c(0, lapply(rnaData, `[[`, 4) |> unlist() |> max())
#rnaRange <- c(0, 5000) #APOE 
## Plot RNA
rnaS <- pmap(list(rnaData, rna_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = rnaRange, linecolor = colors[8], fill = colors[8]))

## Plot Genes
geneTrack <- plotGenes(params = pst, height = 0.25, y=rna_ypos+pst$height+0.05, geneHighlights = data.frame(gene = "APOE", color = "#158C7E"))

## Plot coords
genCoord <- annoGenomeLabel(geneTrack, x = pst$x, y=rna_ypos+pst$height+0.05+0.25+0.05, scale="bp", fontsize = 8)

## Add assay text
plotText("GWAS", x=4.05, y=2.75, fontface = "bold", fontsize = 8, just = c("left", "center"))
plotText("MPRA", x=4.05, y=3.35, fontface = "bold", fontsize = 8, just = c("left", "center"))
plotText("eQTL", x=4.05, y=3.95, fontface = "bold", fontsize = 8, just = c("left", "center"))
plotText("ABC", x=4.05, y=3.95, fontface = "bold", fontsize = 8, just = c("left", "center"))

## Plot timepoint labels
plotText(label = "Hi-C", fontsize = 8, x = params$x+0.04, y = hic_ypos + 0.12, just = c('left', 'bottom'))

## Plot assay labels
plotText("ATAC", x=params$x-0.1, y=mean(atac_ypos+0.1), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[4])
plotText("H3K27ac", x=params$x-0.1, y=mean(h3k_ypos+0.1), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[6])
plotText("RNA", x=params$x-0.1, y=mean(rna_ypos+0.1), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[8])

