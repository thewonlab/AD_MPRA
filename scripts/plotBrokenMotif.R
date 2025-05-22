## Script to generate the motif broken plot from figure 3e

## Load libraries
library(plotgardener)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Biostrings)
library(dplyr)
library(readxl)
library(grid)
library(ggplot2)
library(tidyr)
library(ggh4x)
library(cowplot)

## Read in LIMA data instead 
rnaBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/rna/LIMA_THP1_WT_LPIF_S/proc/signal/", "MERGE", full.names = TRUE)[c(2,16)]
atacBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/atac/LIMA_ATAC_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]
enhBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/chip/H3K27ac/LIMA_h3k27ac_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]


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
  if (adjusted$beta[i]>0){
    adjusted$risk[i] <- adjusted$a2[i]
    adjusted$protective[i] <- adjusted$a1[i]
  } else if (adjusted$beta[i]<0){
    adjusted$risk[i] <- adjusted$a1[i] #SWAP
    adjusted$protective[i] <- adjusted$a2[i]
  }
}

allelic <- adjusted[adjusted$sig==TRUE,]
mpra <- adjusted[,c(10,11,6,13,9)]
colnames(mpra) <- c("chrom", "pos", "p", "snp", "sig")

## Active elements 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
adjusted$active <- FALSE
adjusted[adjusted$id %in% active,]$active <- TRUE

## Read in MPRA data at the individual replicate level: 
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batchvar_N9N6_combined_withControls.txt", header=TRUE, sep = "\t")

## Variants = emvars
variants <- adjusted[adjusted$sig==TRUE & adjusted$active==TRUE,]

## Convert to GRanges: 
buffer <- 0.25E5
putative <- GRanges(seqnames = Rle(variants$seqnames), 
                    ranges = IRanges(start = variants$pos-buffer, end = variants$pos+buffer), 
                    id = variants$id, 
                    rsid = variants$rsid, 
                    pos = variants$pos)
levels <- paste0("chr", 1:22)
levels <- levels[levels %in% seqnames(putative)]
seqlevels(putative) <- levels
putative <- sort(putative)

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
                  ranges = IRanges(start = lonely$pos-1, end = lonely$pos+1), 
                  rsid = lonely$LeadSNPs)
#combine
loci <- c(regular, lonely)
seqlevelsStyle(loci) <- "UCSC"
loci <- resize(loci, 2E6)
loci <- sort(loci)

list <- list()
for (i in 1:length(putative)){
  lead <- loci[nearest(putative[i], loci)]$rsid
  mpra[mpra$snp==lead,]$pos
  putative[i]$pos
  range <- sort(c(mpra[mpra$snp==lead,]$pos, putative[i]$pos))
  if(length(range)==1) {
    range <- c(range, range)
  }
  list[[i]] <-
    GRanges(seqnames = seqnames(putative[i]),
            ranges = IRanges(start = range[1]-buffer, end = range[2]+buffer),
            lead = lead,
            putative = putative[i]$rsid)
}
list <- do.call(`c`, list)
#Remove variants where the lead snp wasn't in our mpra
loci <- list[width(list)>(buffer*2)+1]
putative <- loci

## Get negative controls for boxplots
neg <- batchvar[grep("Scrambled", batchvar$variant),]
neg <- neg[,seq(4,(4*15),by=4)]
negCon <- neg[,1:9]
negTrt <- neg[,10:15]
negCon <- apply(negCon, 2, median) |> 
  as.data.frame()
negTrt <- apply(negTrt, 2, median) |> 
  as.data.frame(names = "activity")
colnames(negCon) <- "activity"
colnames(negTrt) <- "activity"
neg <- rbind(negCon, negTrt)
neg$allele <- "negative"
neg$sample <- rownames(neg)
neg$condition <- c(rep("control", 9), rep("treat", 6))
#neg$label <- neg$sample
neg <- neg[,c(2,4,3,1)]
#neg <- rbind(neg, neg)

## Make boxplots: 
conditionalPlots <- list()
for (i in 1:length(putative)){
  print(i)
  #vars <- putative[i]$putative
  vars <- putative[i]
  ref <- adjusted[adjusted$rsid %in% vars$putative,] |> dplyr::select(id, rsid)
  
  df <- batchvar[batchvar$SNP %in% ref$id,]
  #just get the rna/dna ratios: 
  df <- df[,seq(4,(4*15),by=4)]
  # df$allele <- c("a1", "a2")
  ### get the alleles in the correct order 
  alleles <- 
    strsplit(batchvar[batchvar$SNP %in% ref$id,]$variant, "_") |> 
    lapply(`[[`, 4) |> 
    unlist()
  df$allele <- alleles
  con <- df[,c(1:9, 16)]
  trt <- df[,10:16]
  #add con/trt
  con$condition <- "control"
  trt$condition <- "treat"
  #reshape and combine
  df <- 
    rbind(pivot_longer(con, 
                       cols = !c('condition', 'allele'), 
                       names_to = "sample",
                       values_to = "activity"), 
          pivot_longer(trt, 
                       cols = !c('condition', 'allele'), 
                       names_to = "sample",
                       values_to = "activity"))
  
  ## Determine which is risk and which is protective: 
  beta <- adjusted[adjusted$rsid==ref$rsid,]$beta
  df$alleleClass <- ""
  if (beta>0){
    df[df$allele %in% alleles[1],]$alleleClass <- "risk"
    df[df$allele %in% alleles[2],]$alleleClass <- "protective"
  } else if (beta<0) {
    df[df$allele %in% alleles[2],]$alleleClass <- "risk"
    df[df$allele %in% alleles[1],]$alleleClass <- "protective"
  }
  
  ## Add risk/protective information: 
  # df[df$allele=="a1",]$allele <- "protective"
  # df[df$allele=="a2",]$allele <- "risk"
  
  ## Add allele bp
  # al <- adjusted[adjusted$rsid %in% vars,]
  # al <- al[,c(1,19,20)]
  # #al <- al[,c(1,21,22)]
  # colnames(al)[1] <- "label"
  # al <- pivot_longer(al, -label, names_to = "allele", values_to = "bp")
  
  # df$bp <- ""
  # df[df$allele=="protective",]$bp <- al[al$allele=="protective",]$bp
  # df[df$allele=="risk",]$bp <- al[al$allele=="risk",]$bp
  # df$name <- paste0(df$allele, " (", df$bp, ")")
  
  df$name <- paste0(df$alleleClass, " (", df$allele, ")")
  
  df <- df[,c(1:4,6)]
  neg$name <- "neg"
  
  df <- rbind(df, neg)
  mpraBox <- 
    ggplot(df, aes(x = name, y = activity, fill = name)) + 
    facet_wrap(~condition, strip.position = "bottom") + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(fill = name), size = 2, shape = 21, position = position_jitterdodge(jitter.width = 0.01)) + 
    scale_fill_manual(values = c("#D2D7D9", "#B6A4F2", "#B82CE6")) + 
    theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(vjust = 20, size = 5), 
                            theme(plot.margin = unit(c(-1, 0.1, 0.2, 0.1), "cm")), 
                            strip.text = element_text(size = 8, vjust = -2), 
                            plot.title = element_text(size = 12)) + 
    xlab("") + ylab("RNA/DNA") + ylim(round(range(df$activity))) + 
    scale_x_discrete(limits = rev)
  
  conditionalPlots[[i]] <- mpraBox
}

## Set parameters
colors <- brewer.pal(8, "YlGnBu")
tp <- c("0", "24")

## Initiate plotgardener page 
pageCreate(width = 8.75, height = 5, showGuides = FALSE, xgrid = 0.0, ygrid = 0.0)

## Set (i) to be the variant of interest
## i = 37 is the variant rs1871047 we show in figure 3e
## Set parameters 
params <- 
  pgParams(chrom = as.character(seqnames(putative)[i]),
           chromstart = start(putative)[i]+20000,
           chromend = end(putative)[i]-10000,
           assembly = "hg19",
           x = 2.1,
           width = 2,
           height = 0.5)

pst <- params
pst$height <- 0.15

#Y positions: 
atac_ypos <- c((1.8),(1.8)+pst$height+0.05)
h3k_ypos <- c((atac_ypos[2]+pst$height+0.05),(atac_ypos[2]+pst$height+0.05)+pst$height+0.05)
rna_ypos <- c((h3k_ypos[2]+pst$height+0.05),(h3k_ypos[2]+pst$height+0.05)+pst$height+0.05)

## GWAS
gwasMax <- gwas[gwas$chrom == params$chrom & gwas$pos>params$chromstart & gwas$pos< params$chromend,]$p
gwasMax <- gwasMax[gwasMax>0]
gwasMax <- log10(gwasMax) * -1
gwasMax <- max(gwasMax) |> round()

gwasMan <- plotManhattan(data = gwas, 
                         params = params,
                         range = c(0,gwasMax+5),
                         y = 0.5, height = 0.5, 
                         sigLine = TRUE, col = "gray90", fill = "#B2BBBF",
                         baseline = TRUE, baseline.lwd = 2, baseline.color = "black",
                         snpHighlights = data.frame(snp = c(putative[i]$putative),
                                                    pch = c(23),
                                                    cex = c(0.75),
                                                    col = c("black"), 
                                                    fill = c("black")))
annoYaxis(
  plot = gwasMan,
  at = c(seq(0, gwasMax, by = gwasMax %/% 5), gwasMax+5),
  axisLine = TRUE, fontsize = 8
)

plotText(
  label = "-log10(p)", x = 1.65, y = 0.75, rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

## MPRA
mpraMax <- mpra[mpra$chrom == params$chrom & mpra$pos>params$chromstart & mpra$pos< params$chromend,]$p |> na.omit()
mpraMax <- log10(mpraMax) * -1
mpraMax <- max(mpraMax) |> ceiling()
mpraMan <- plotManhattan(mpra,
                         params = params,
                         y = 1.2, height = 0.5,
                         sigLine = TRUE, col = "gray90", sigVal = 0.05, fill = "#B2BBBF",
                         baseline = TRUE, baseline.lwd = 2, baseline.color = "black", range = c(0,mpraMax+0.5),
                         snpHighlights = data.frame(snp = c(loci[i]$putative),
                                                    pch = c(23),
                                                    cex = c(0.75),
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
  label = "-log10(padj)", x = 1.7, y = 1.45, rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

## ATAC 
atacData <- lapply(atacBW, \(x) readBigwig(file = x, params = pst)[c(1:3, 6)])
atacRange <- c(0, lapply(atacData, `[[`, 4) |> unlist() |> max())
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
## Plot RNA
rnaS <- pmap(list(rnaData, rna_ypos), \(x, y) plotSignal(data = x, params = pst, y = y, range = rnaRange, linecolor = colors[8], fill = colors[8]))

## Plot Genes
geneTrack <- plotGenes(params = params, height = 0.25, y=rna_ypos[2]+pst$height+0.05)

## Plot coords
genCoord <- annoGenomeLabel(geneTrack, x = params$x, y=rna_ypos[2]+pst$height+0.05+0.25+0.05, scale="bp", fontsize = 8)

## Plot timepoint labels
plotText(label = tp, fontsize = 6, fontcolor = colors[4], x = params$x+0.03, y = atac_ypos + 0.03, just = c('left', 'bottom'))
plotText(label = tp, fontsize = 6, fontcolor = colors[6], x = params$x+0.03, y = h3k_ypos + 0.03, just = c('left', 'bottom'))
plotText(label = tp, fontsize = 6, fontcolor = colors[8], x = params$x+0.03, y = rna_ypos + 0.03, just = c('left', 'bottom'))

## Plot assay labels
plotText("GWAS", x=2.15, y=0.5, fontface = "bold", fontsize = 8, just = c("left", "center"))
plotText("MPRA", x=2.15, y=1.2, fontface = "bold", fontsize = 6, just = c("left", "center"))
plotText("ATAC", x=params$x-0.1, y=mean(atac_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[4])
plotText("H3K27ac", x=params$x-0.1, y=mean(h3k_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[6])
plotText("RNA", x=params$x-0.1, y=mean(rna_ypos+0.04), fontface = "bold", fontsize = 6, rot = 90, fontcolor = colors[8])

## Add boxplot 
plotGG(conditionalPlots[[i]], x = 4.125, y = 1.5, height = 2.5, width = 2.5)

## Add a highlight bar: 
pos <- adjusted[adjusted$rsid==putative$putative[i],]$pos
annoHighlight(gwasMan, chrom = params$chrom, 
              chromstart = pos - 125, chromend = pos + 125, y = 0.5, height = 2.8)

