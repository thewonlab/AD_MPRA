## Figure 1 locus zoom script

## Survey putative causal variants with just GWAS and MPRA data: 

## Load libraries
install_github(repo = "PhanstielLab/plotgardener", ref = "plotManhattan_features") #specific branch of plotgardener 
library(plotgardener)
library(data.table)
library(readxl)
library(ggplot2)
library(tidyr)
library(ggh4x)
library(cowplot)
library(GenomicRanges)
library(dplyr)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)

## GWAS: 
gwas <- fread("/work/users/m/a/marielle/gwas/jansen/AD_sumstats_Jansenetal_2019sept_filtered_0.05.txt.gz") ## Jansen 2019 summary statistics
gwas <- gwas[,c(2,3,8,6)]
colnames(gwas) <- c("chrom", "pos", "p", "snp")
gwas$chrom <- paste0("chr", gwas$chrom)

## MPRA: 
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T) ## results from MPRA_allelic.R
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

## allelic
allelic <- adjusted[adjusted$sig==TRUE,]
mpra <- adjusted[,c(10,11,6,13,9)]
colnames(mpra) <- c("chrom", "pos", "p", "snp", "sig")

## activity
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
adjusted$active 

## Add activity: 
adjusted$activity <- ""
adjusted[adjusted$id %in% active,]$activity <- TRUE
adjusted[adjusted$activity=="",]$activity <- FALSE

## Select variants to look at--these are the two that we use in figure 1
variants <- adjusted[adjusted$rsid %in% c("rs1151105", "rs6979218"),]

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

## Convert to GRanges: 
buffer <- 0.2E6
putative <- GRanges(seqnames = Rle(variants$seqnames), 
                    ranges = IRanges(start = variants$pos-buffer, end = variants$pos+buffer), 
                    id = variants$id, 
                    rsid = variants$rsid, 
                    pos = variants$pos)
levels <- paste0("chr", 1:22)
levels <- levels[levels %in% seqnames(putative)]
seqlevels(putative) <- levels
putative <- sort(putative)

## Determine which GWAS locus each of the variants belongs to 
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

## Read in MPRA data at the individual replicate level: 
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batchvar_N9N6_combined_withControls.txt", header=TRUE, sep = "\t") ## available through GEO
## Get negative controls: 
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
neg$condition <- c(rep("resting", 9), rep("activated", 6))
neg$label <- neg$sample
neg <- neg[,c(2,4,3,1,5)]
neg <- rbind(neg, neg)
neg$class <- c(rep("MPRA", 15), rep("GWAS", 15))

## Make the locus plots but in ggplot along with the boxplots 
plots <- list()
for (i in 1:length(loci)){
  print(i)
  vars <- c(loci[i]$lead, loci[i]$putative)
  ref <- adjusted[adjusted$rsid %in% vars,] |> dplyr::select(id, rsid)
  ref$class <- ""
  ref[ref$rsid==loci[i]$lead,]$class <- "GWAS"
  ref[ref$rsid==loci[i]$putative,]$class <- "MPRA"
  
  data <- data.frame()
  for (v in ref$id){
    df <- batchvar[batchvar$SNP %in% v,]
    #just get the rna/dna ratios: 
    df <- df[,seq(4,(4*15),by=4)]
    df$allele <- c("a1", "a2")
    con <- df[,c(1:9, 16)]
    trt <- df[,10:16]
    #add con/trt
    con$condition <- "resting"
    trt$condition <- "activated"
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
    df$label <- v
    df$class <- ""
    if (v==ref[ref$class=="GWAS",]$id){
      df$class <- "GWAS"
    } else if (v==ref[ref$class=="MPRA",]$id){
      df$class <- "MPRA"
    }
    data <- rbind(data, df)
  }
  
  ## Add risk/protective information: 
  data[data$allele=="a1",]$allele <- "protect"
  data[data$allele=="a2",]$allele <- "risk"
  
  ## Add allele bp
  al <- adjusted[adjusted$rsid %in% vars,]
  al <- al[,c(1,19,20)]
  colnames(al)[1] <- "label"
  al <- pivot_longer(al, -label, names_to = "allele", values_to = "bp")
  
  data$bp <- ""
  vars <- unique(data$label)
  data[data$allele=="protect" & data$label==vars[1],]$bp <- al[al$allele=="protective" & al$label==vars[1],]$bp
  data[data$allele=="risk" & data$label==vars[1],]$bp <- al[al$allele=="risk" & al$label==vars[1],]$bp
  data[data$allele=="protect" & data$label==vars[2],]$bp <- al[al$allele=="protective" & al$label==vars[2],]$bp
  data[data$allele=="risk" & data$label==vars[2],]$bp <- al[al$allele=="risk" & al$label==vars[2],]$bp
  data$name <- paste0(data$allele, " (", data$bp, ")")
  
  data <- data[,c(1:6,8)]
  neg$name <- "neg"
  
  data <- rbind(data, neg)
  data$condition <- factor(data$condition, levels = c("resting", "activated"))
  
  gwasBox <- 
    ggplot(data[data$class=="GWAS",], aes(x = name, y = activity, fill = name)) + 
    facet_wrap(~condition, strip.position = "bottom") + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(fill = name), size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.01)) + 
    scale_fill_manual(values = c("#D2D7D9", "#C0ADFF", "#C22EF2")) + 
    theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(vjust = 26, size = 5), 
                            theme(plot.margin = unit(c(0.1, 0.1, -1, 0.1), "cm")), 
                            strip.text = element_text(size = 8, vjust = 0), 
                            plot.title = element_text(size = 12)) + 
    xlab("") + ylab("RNA/DNA") + ylim(round(range(data$activity))) + 
    scale_x_discrete(limits = rev)
  
  mpraBox <- 
    ggplot(data[data$class=="MPRA",], aes(x = name, y = activity, fill = name)) + 
    facet_wrap(~condition, strip.position = "bottom") + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(aes(fill = name), size = 1.5, shape = 21, position = position_jitterdodge(jitter.width = 0.01)) + 
    scale_fill_manual(values = c("#D2D7D9", "#C0ADFF", "#C22EF2")) + 
    theme_minimal() + theme(legend.position = "none", axis.text.x = element_text(vjust = 26, size = 5), 
                            theme(plot.margin = unit(c(-1, 0.1, 0.2, 0.1), "cm")), 
                            strip.text = element_text(size = 8, vjust = 0), 
                            plot.title = element_text(size = 12)) + 
    xlab("") + ylab("RNA/DNA") + ylim(round(range(data$activity))) + 
    scale_x_discrete(limits = rev)
  
  ## LOCUS ZOOM
  params <-
    pgParams(chrom = as.character(seqnames(loci)[i]),
             chromstart = start(loci)[i],
             chromend = end(loci)[i],
             assembly = "hg19",
             x = 0.5,
             width = 2,
             height = 1)
  
  ##GWAS
  gwasData <- gwas[gwas$chrom == params$chrom & gwas$pos>params$chromstart & gwas$pos< params$chromend,]
  gwasMan <- 
    ggplot(gwasData, aes(x = pos, y = -log10(p))) + 
    geom_point(color = "#D2D7D9", fill = "#D2D7D9", cex = 1) + 
    geom_point(data = gwasData[gwasData$snp==loci$lead[i]], color = "black", fill = "black", pch = 24, cex = 2) + 
    geom_point(data = gwasData[gwasData$snp==loci$putative[i]], color = "black", fill = "black", pch = 23, cex = 2) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          plot.background = element_rect(fill = "transparent", color = NA)) + 
    xlab("") + 
    xlim(start(loci)[i], end(loci)[i]) + 
    geom_hline(yintercept = 7, lty = 3)
  
  ### MPRA GGPLOT VERSION WITH PAIRED ALLELES: 
  ## First get all of the variants that are within the locus
  vars <- mpra[mpra$chrom == params$chrom & mpra$pos>params$chromstart & mpra$pos< params$chromend,]$snp
  vars <- adjusted[adjusted$rsid %in% vars,]$id
  df <- batchvar[batchvar$SNP %in% vars,]
  #just get the rna/dna ratios: 
  rna <- df[,seq(2,(4*15),by=4)]
  dna <- df[,seq(3,(4*15),by=4)]
  #just get treated: 
  rnaCon <- rna[,1:9]
  dnaCon <- dna[,1:9]
  rnaTrt <- rna[,10:15]
  dnaTrt <- dna[,10:15]
  ratio <- apply(rna, 1, sum)/apply(dna, 1, sum)
  ratioCon <- apply(rnaCon, 1, sum)/apply(dnaCon, 1, sum)
  ratioTrt <- apply(rnaTrt, 1, sum)/apply(dnaTrt, 1, sum)
  df$ratio <- ratio
  df$ratioCon <- ratioCon
  df$ratioTrt <- ratioTrt
  df$allele <- c("a1", "a2")
  df <- df[,c(1,62,63,64,65,66)]
  df$delta <- df$ratioTrt-df$ratioCon
  df$lfc <- log(df$ratioTrt/df$ratioCon)
  df <- na.omit(df)
  ## separate position 
  df$position <- unlist(strsplit(df$SNP, "_"))[seq(2,(2*nrow(df)), by=2)] |> as.numeric()
  leadID <- adjusted[adjusted$rsid==loci$lead[i],]$id 
  mpraID <- adjusted[adjusted$rsid==loci$putative[i],]$id 
  
  mpraMan <- ggplot(df, aes(x = position, y = ratio, color = allele)) + 
    geom_line(aes(group = SNP), color = "#D2D7D9") + 
    geom_point(cex = 1) + 
    scale_color_manual(values = c("#C0ADFF", "#C22EF2")) + 
    geom_point(data = df[df$SNP==leadID,], fill = "black", pch = 24, cex = 2) + 
    geom_point(data = df[df$SNP==mpraID,], fill = "black", pch = 23, cex = 2) + 
    theme_classic() + 
    theme(legend.position = "none", 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          axis.title.y = element_text(size = 8)) + 
    xlab("") + 
    #ylab("Log Fold Change(Activated \nRNA/DNA to Resting RNA/DNA)") + 
    ylab("RNA/DNA") + 
    #ylim(0.42,2) + 
    xlim(start(loci)[i], end(loci)[i]) + 
    geom_hline(yintercept = 0, lty = 3)
  
  
  ## Gene tracks
  hg19_txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  seqlevels(hg19_txdb) <- params$chrom
  genes <- genes(hg19_txdb, columns = "GENEID")
  genes <- genes[start(genes)>params$chromstart & end(genes)<params$chromend]
  genes <- unique(genes)
  genes$SYMBOL <- mapIds(org.Hs.eg.db, keys = as.character(genes$GENEID), keytype = "ENTREZID", column = "SYMBOL")
  genes <- as.data.frame(genes)
  genes$x <- "+"
  genes[genes$strand=="-",]$x <- "-"
  genes$mid <- (genes$start+genes$end)/2
  if (i == 2){
    goi <- c("MS4A2", "MS4A6A", "MS4A6E")
  } else if (i == 1){
    goi <- c("STAG3", "GPC2", "SPDYE3", "ZCWPW1", "MEPCE")
  }
  genes$label <- ""
  genes[genes$SYMBOL %in% goi,]$label <- goi
  
  geneTrack <- 
    ggplot(genes) + 
    geom_linerange(data = genes[genes$strand=="+",], aes(xmin = start, xmax = end, y = x), color = "#D2D7D9", size = 2) + 
    geom_linerange(data = genes[genes$strand=="-",], aes(xmin = end, xmax = start, y = x), color = "#D2D7D9", size = 2) + 
    xlim(start(loci)[i], end(loci)[i]) + 
    theme_classic() + 
    geom_text(aes(x = mid, y = x, label = label, color = strand), vjust = -1, size = 2.5) + 
    #geom_text_repel(aes(x = mid, y = x, label = SYMBOL, color = strand, size = 2.5)) + 
    scale_color_manual(values = c("#D2D7D9", "#D2D7D9")) + 
    theme(legend.position = "none", 
          axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 10, face = "bold", color = c("#D2D7D9", "#D2D7D9")), 
          plot.background = element_rect(fill = "transparent", color = NA)) + 
    xlab("") + ylab("") + 
    annotate("text", x = (params$chromstart+params$chromend)/2, y = 0.5, label = params$chrom, size = 2.5) + 
    annotate("text", x = params$chromstart, y = 0.5, label = params$chromstart, size = 2.5) +
    annotate("text", x = params$chromend, y = 0.5, label = params$chromend, size = 2.5)
  
  ## combine plots
  combinedPlots <- plot_grid(gwasMan, NULL, mpraMan, NULL, geneTrack, gwasBox, NULL, mpraBox, NULL, rel_heights = c(1, -0.25, 1, -0.35, 0.7, 1, -0.05, 1, -0.1), nrow=5, byrow = FALSE, align = "hv", axis = "trbl")
  
  ## save to list
  plots[[i]] <- combinedPlots
  
}

## Resulting figures
plots[[1]] #rs1151105 from fig 1h
plots[[2]] #rs6979218 from fig 1e


