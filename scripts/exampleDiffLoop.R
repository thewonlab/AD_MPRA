## Script to generate a plotgardener figure for a differential loop with concordant increases in ATAC/K27/RNA
## Generates figure 2b

## Load libraries
library(GenomicRanges)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(purrr)
library(RColorBrewer)
library(data.table)
library(mariner)
library(stringr)

## Read in LIMA data instead (all available from Reed et al 2022 Cell Reports GEO superseries GSE201376)
hicFiles <- c("/proj/phanstiel_lab/Data/processed/LIMA/hic/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap/aligned/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_inter.hic", 
              "/proj/phanstiel_lab/Data/processed/LIMA/hic/LIMA_THP1_WT_LPIF_1440_S_0.0.0_megaMap/aligned/LIMA_THP1_WT_LPIF_1440_S_0.0.0_megaMap_inter.hic")
rnaBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/rna/LIMA_THP1_WT_LPIF_S/proc/signal/", "MERGE", full.names = TRUE)[c(2,16)]
atacBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/atac/LIMA_ATAC_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]
enhBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/chip/H3K27ac/LIMA_h3k27ac_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]

## Read in loops
## Loops: 
loops <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/hic/hicLFC.txt") #list of differential loops 
loops$id <- rownames(loops)

## Create plotgardener page
pageCreate(width = 8.75, height = 5, showGuides = FALSE, xgrid = 0.0, ygrid = 0.0)

## Set parameters for differential loop of interest
params <-
  pgParams(chrom = "6",
           chromstart = 113820000,
           chromend = 114300000,
           assembly = "hg19",
           x = 2.1,
           width = 2,
           height = 0.6)

## different parameters for signal tracks
pst <- params
pst$chrom <- paste0("chr", pst$chrom)
pst$height <- 0.1

#Y positions: 
hic_ypos <- c(0.1, 0.1+params$height+0.05)
atac_ypos <- c((hic_ypos[2]+params$height+0.05),(hic_ypos[2]+params$height+0.05)+pst$height+0.05)
h3k_ypos <- c((atac_ypos[2]+pst$height+0.05),(atac_ypos[2]+pst$height+0.05)+pst$height+0.05)
rna_ypos <- c((h3k_ypos[2]+pst$height+0.05),(h3k_ypos[2]+pst$height+0.05)+pst$height+0.05)

## Plot HiC
hicPlots <- pmap(list(hicFiles, hic_ypos), \(x, y) plotHicRectangle(data = x, params = params, y = y, zrange = c(0,80)))

## Highlight the differential loop: (this is the differential loop shown in figure 2b)
pixels <- annoPixels(
  plot = hicPlots[[1]], data = loops[loops$id=="loop18395",], type = "circle"
)
pixels <- annoPixels(
  plot = hicPlots[[2]], data = loops[loops$id=="loop18395",], type = "circle"
)

## ATAC 
atacData <- lapply(atacBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
atacRange <- c(0, lapply(atacData, `[[`, 4) |> unlist() |> max())
#atacRange <- c(0,1000) #change to custom because one peak is reall reall tall and static 
## Plot ATAC
atacS <- pmap(list(atacData, atac_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = atacRange, linecolor = colors[4], fill = colors[4]))

### H3K27ac
h3kData <- lapply(enhBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
h3kRange <- c(0, lapply(h3kData, `[[`, 4) |> unlist() |> max())
## Plot H3K27ac
h3k27S <- pmap(list(h3kData, h3k_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = h3kRange, linecolor = colors[6], fill = colors[6]))

## RNA
rnaData <- lapply(rnaBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
rnaRange <- c(0, lapply(rnaData, `[[`, 4) |> unlist() |> max())
#rnaRange <- c(0,20000) #change to custom because one peak is reall reall tall and static 
## Plot RNA
rnaS <- pmap(list(rnaData, rna_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = rnaRange, linecolor = colors[8], fill = colors[8]))

## Plot Genes
geneTrack <- plotGenes(params = pst, height = 0.25, y=rna_ypos[2]+pst$height+0.05)

## Plot coords
genCoord <- annoGenomeLabel(geneTrack, x = params$x, y=rna_ypos[2]+pst$height+0.05+0.25+0.05, scale="bp", fontsize = 8)

## Add highlight bars 
annoHighlight(
  plot = atacS[[1]],
  chrom = "chr6",
  chromstart = 113950000, chromend = 113960000,
  y = 1.35, height = 1.25, just = c("left", "top"),
  default.units = "inches"
)

annoHighlight(
  plot = atacS[[1]],
  chrom = "chr6",
  chromstart = 114170000, chromend = 114180000,
  y = 1.35, height = 1.25, just = c("left", "top"),
  default.units = "inches"
)

## Plot timepoint labels
plotText(label = tp, fontsize = 6, x = params$x+0.03, y = hic_ypos + 0.1, just = c('left', 'bottom'))
plotText(label = tp, fontsize = 6, fontcolor = colors[4], x = params$x+0.03, y = atac_ypos + 0.03, just = c('left', 'bottom'))
plotText(label = tp, fontsize = 6, fontcolor = colors[6], x = params$x+0.03, y = h3k_ypos + 0.03, just = c('left', 'bottom'))
plotText(label = tp, fontsize = 6, fontcolor = colors[8], x = params$x+0.03, y = rna_ypos + 0.03, just = c('left', 'bottom'))


## Plot assay labels
plotText("ATAC", x=params$x-0.1, y=mean(atac_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[4])
plotText("H3K27ac", x=params$x-0.1, y=mean(h3k_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[6])
plotText("RNA", x=params$x-0.1, y=mean(rna_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[8])
