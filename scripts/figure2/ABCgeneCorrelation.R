### Script to use all ABC pairs in at least one sample at DEGs specifically to figure out if Up and Down DEGs have diff ABC scores: 
## Generates figure 2h

## Library loading
library(mariner)
library(plotgardener)
library(data.table)
library(GenomicRanges)
library(InteractionSet)
library(nullranges)

## Bigwigs from Reed et al 2022
atacBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/atac/LIMA_ATAC_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]
enhBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/chip/H3K27ac/LIMA_h3k27ac_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]

#Hi-C files from Reed et al 2022
hicfiles <- c("/proj/phanstiel_lab/Data/processed/LIMA/hic/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap/aligned/LIMA_THP1_WT_LPIF_0000_S_0.0.0_megaMap_inter.hic", 
              "/proj/phanstiel_lab/Data/processed/LIMA/hic/LIMA_THP1_WT_LPIF_1440_S_0.0.0_megaMap/aligned/LIMA_THP1_WT_LPIF_1440_S_0.0.0_megaMap_inter.hic")

## ABC output from the supplemental information in this paper 
abc <- list.files("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/Neighborhoods", "Pred", full.names = TRUE)[c(2,4)] |> 
  lapply(fread)
## Reformat
for (i in 1:length(abc)){
  df <- abc[[i]]
  #df <- df[df$seqnames1=="chr22"]
  df <- df[,c(1:3,6:8,11:13,15:21)] |> 
    na.omit() |> 
    as_ginteractions() |> 
    snapToBins(5e3)
  df$pair <- paste0(df$anchor1.id, "-", df$anchor2.TargetGene)
  df <- df[df$anchor2.TargetGeneIsExpressed == TRUE]
  ## do more stringent ABC score filtering: 
  df <- df[df$ABCscore>0.05]
  abc[[i]] <- df
}

abc0 <- abc[[1]]
abc24 <- abc[[2]]
abc <- c(abc0, abc24)

pairs <- (unique(abc0$pair, abc24$pair))


## Try subsetting for the genes that are DIFFERENTIAL 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt") #LIMA
rna$ENSEMBL <- rownames(rna)
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
rna <- merge(genes, rna, by = "ENSEMBL")
rna$basemean <- (rna$X0 + rna$X24)/2
upGenes <- rna[rna$class=="gained",]$SYMBOL
downGenes <- rna[rna$class=="lost",]$SYMBOL
deg <- c(upGenes, downGenes)

## Just get the DEGs that are in the abc 
deg <- deg[deg %in% abc$anchor2.TargetGene]

## Select the same number of DEGs but static genes 
set.seed(718)
static <- matchRanges(focal = rna[!rna$class=="static",], 
                      pool = rna[rna$class=="static",], 
                      covar = ~basemean, 
                      method = "stratified", replace = F)
static <- static[static$SYMBOL %in% abc$anchor2.TargetGene,]$SYMBOL |> unique()


## For all up DEGs, calculate ABC for all possible enhancers in both contexts: 

#degRanges <- data.frame()
staticRanges <- data.frame()
for (i in 1:length(static)){
  print(i)
  gene <- static[i]
  enhCoord <- anchors(abc[abc$anchor2.TargetGene %in% gene], "first")
  enhCoord$id <- abc[abc$anchor2.TargetGene %in% gene]$anchor1.id
  geneCoord <- anchors(abc[abc$anchor2.TargetGene %in% gene], "second")
  geneCoord$gene <- gene
  
  ## For each of the enhancers, extract the H3K27ac and ATAC scores from EACH CONTEXT: 
  enhActivity <- GRanges()
  for (x in 1:length(enhCoord)){
    ## Get Activity
    enh <- enhCoord[x]
    enh$atac0 <- readBigwig(atacBW[1], chrom = seqnames(enh), chromstart = start(enh), chromend = end(enh))$score |> sum()
    enh$atac24 <- readBigwig(atacBW[2], chrom = seqnames(enh), chromstart = start(enh), chromend = end(enh))$score |> sum()
    enh$enh0 <- readBigwig(enhBW[1], chrom = seqnames(enh), chromstart = start(enh), chromend = end(enh))$score |> sum()
    enh$enh24 <- readBigwig(enhBW[2], chrom = seqnames(enh), chromstart = start(enh), chromend = end(enh))$score |> sum()
    enhActivity <- c(enhActivity, enh)
  }
  ## Now get the HiC contacts for each enhancer with the gene: 
  pairs <- GInteractions(enhActivity, geneCoord) |> as.data.frame()
  pairs$class <- ""
  if (gene %in% upGenes) {
    pairs$class <- "up"
  } else {
    pairs$class <- "down"
  }
  
  staticRanges <- rbind(staticRanges, pairs)
  
}

## Write the output to a txt file 
#write.table(degRanges, file = "/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/abc_deg_activity.txt", row.names = F, sep = "\t", quote = F)
write.table(staticRanges, file = "/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/abc_static_activity.txt", row.names = F, sep = "\t", quote = F)



### If running at a later date, can read in the output written above 
degRanges <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/abc_deg_activity.txt", sep = "\t", header = T)
staticRanges <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/abc_static_activity.txt", sep = "\t", header = T)
staticRanges$class <- "static"

## Combine all data
ranges <- rbind(degRanges, staticRanges)

## Get the HiC counts at each of these ranges:
pairs <- assignToBins(ranges, binSize = 5e3)
seqlevelsStyle(pairs) <- "ENSEMBL" 
pixelsNorm <- pullHicPixels(x = pairs,
                            files = hicfiles,
                            binSize = 5e3,
                            blockSize = 500e6,
                            norm = "KR")

counts <- counts(pixelsNorm) |> as.data.frame()
pairs$contact0 <- counts[,1]
pairs$contact24 <- counts[,2]
pairs <- as.data.frame(pairs)
pairs <- pairs[,c(1:17,19,20,18)]


## Save the output: 
write.table(pairs, file = "/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/all_genes_activity_contact.txt", quote=F, sep = "\t", row.names = F)


## If running at a later date, can read in the output from above: 
pairs <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/all_genes_activity_contact.txt", header=T)

## So for every single DEG, get all of the pairs, calculate an ABC score for each of the pairs, and then sum those scores together
## That's the score for each gene, then we divide the restingABC and activatedABC to get a diff ABC score
## Partition by which are at resting DEGs and which are at activated DEGs

results <- data.frame()
for (i in 1:length(unique(pairs$anchor2.gene))){
  print(i)
  gene <- unique(pairs$anchor2.gene)[i]
  sub <- pairs[pairs$anchor2.gene==gene,]
  ## separate the contexts:
  sub0 <- sub[,c(11,12,14,18)]
  sub24 <- sub[,c(11,13,15,19)]
  ## Get the sum ABC score for each condition
  if (nrow(sub0)==1){
    score0 <-
      apply(sub0[,2:4], 2, log)
    score0 <- prod(score0[1], score0[2], score0[3])
    score24 <-
      apply(sub24[,2:4], 2, log)
    score24 <- prod(score24[1], score24[2], score24[3])
    ## Get the delta ABC score
    deltaABC <- log(score24/score0)
  } else {
    score0 <-
      apply(sub0[,2:4], 2, log) |>
      apply(1, prod) |>
      sum()
    score24 <-
      apply(sub24[,2:4], 2, log) |>
      apply(1, prod) |>
      sum()
  }
  ## Get the delta ABC score
  deltaABC <- log(score24/score0)
  ## Make the row
  row <- data.frame(gene = gene,
                    class = sub$class[1],
                    deltaABC = deltaABC)
  results <- rbind(results, row)
}


results <- results[!results$deltaABC=="-Inf",]
results <- results[!is.na(results$deltaABC),]

## Density plots from figure 2h
abcRNA_density <- 
  ggplot(results, aes(x = deltaABC, fill = class)) +
  geom_density(alpha = 0.8, color = NA) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) +
  geom_vline(xintercept = 0, lty = 3) +
  scale_fill_manual(values = c("#F0B0AA", "#D7DDE6",  "#EC6F61")) +
  xlab("Differential genes at ABC anchors") +
  ylab("log(activatedABC/restingABC)") + 
  coord_cartesian(xlim = c(-0.5,0.5))

## stats: 
wilcox.test(results[results$class=="up",]$deltaABC, results[results$class=="static",]$deltaABC)$p.value
wilcox.test(results[results$class=="down",]$deltaABC, results[results$class=="static",]$deltaABC)$p.value

## Violin plots: 
abcRNA_violin <- 
  ggplot(results, aes(x = class, y = deltaABC, fill = class)) +
  geom_violin(color = NA) +
  geom_boxplot(outlier.shape = NA, width = 0.1) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) +
  geom_hline(yintercept = 0, lty = 3) +
  scale_fill_manual(values = c("#FAA885", "gray",  "#F17D7E")) +
  xlab("Differential genes at ABC anchors") +
  ylab("log(activatedABC/restingABC)")


## Try looking at the relationship between delta abc and RNA fc 
data <- left_join(results, rna, by = c("gene" ="SYMBOL"))
data <- na.omit(data)

cor <- cor(data$log2FoldChange, data$deltaABC)^2
abcRNA_scatter <- 
  ggplot(data, aes(x = log2FoldChange, y = deltaABC, color = class.x)) + 
  geom_point() + 
  scale_color_manual(values = c("#FAA885", "gray", "#F17D7E")) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  geom_vline(xintercept = 0, lty = 3) + 
  geom_hline(yintercept = 0, lty = 3) + 
  xlab("log(activatedRNA/restingRNA)") + 
  ylab("log(activatedABC/restingABC)") + 
  geom_smooth(method = "lm", color = "black", lwd = 0.5)
save(abcRNA_boxes, abcRNA_density, abcRNA_violin, abcRNA_scatter, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/abcRNA_comparison.rda")

## Alternatively, we should be able to use the full ABC results in order to answer this same question

## ABC FILES LIMA: 
abc <- list.files("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/Neighborhoods", "Full", full.names = TRUE)[c(1,2)] |> 
  lapply(fread)

## reformat 
for (i in 1:length(abc)){
  df <- abc[[i]]
  #df <- df[df$seqnames1=="chr22"]
  df <- df[,c(1:3,6:8,11:13,15:21)] |> 
    na.omit() |> 
    as_ginteractions() |> 
    snapToBins(5e3)
  df$pair <- paste0(df$anchor1.id, "-", df$anchor2.TargetGene)
  df <- df[df$anchor2.TargetGeneIsExpressed == TRUE]
  ## do more stringent ABC score filtering: 
  #df <- df[df$ABCscore>0]
  abc[[i]] <- as.data.frame(df)
}

## Just get our genes of interest
# for (i in 1:length(abc)){
#   df <- abc[[i]]
#   abc[[i]] <- df[df$anchor2.TargetGene %in% c(deg, static)]
# }

## separate
abc0 <- abc[[1]]
abc24 <- abc[[2]]

## Just get the pairs that have at least a 0.05 score in at least one condition
valid0 <- abc0[abc0$ABCscore>0.05,]$pair
valid24 <- abc24[abc24$ABCscore>0.05,]$pair
validpairs <- unique(c(valid0, valid24))


## Make sure they both have the same pairs 
abc0 <- abc0[abc0$pair %in% validpairs,]
abc24 <- abc24[abc24$pair %in% validpairs,]
abc0 <- abc0[abc0$pair %in% abc24$pair,]
abc24 <- abc24[abc24$pair %in% abc0$pair,]
genes <- unique(c(abc0$anchor2.TargetGene, abc24$anchor2.TargetGene))
genes <- c(static, deg)
genes <- genes[genes %in% c(abc0$anchor2.TargetGene, abc24$anchor2.TargetGene)]

## calc abc scores per each gene in both conditions 
results <- data.frame()
for (i in 1:length(genes)){
  print(i)
  gene <- genes[i]
  ## separate 0 and 24 
  #score0 <- abc0[abc0$anchor2.TargetGene==gene,]$ABCnumerator |> sum()
  #score24 <- abc24[abc24$anchor2.TargetGene==gene,]$ABCnumerator |> sum()
  
  sub0 <- abc0[abc0$anchor2.TargetGene==gene,]
  sub24 <- abc24[abc24$anchor2.TargetGene==gene,]
  
  ## Get the sum ABC score for each condition
  if (nrow(sub0)==1){
    score0 <-
      apply(sub0[,c(13,18)], 2, log)
    score0 <- prod(score0[1], score0[2])
    score24 <-
      apply(sub24[,c(13,18)], 2, log)
    score24 <- prod(score24[1], score24[2])
    ## Get the delta ABC score
    deltaABC <- log(score24/score0)
  } else {
    score0 <-
      apply(sub0[,c(13,18)], 2, log) |>
      apply(1, prod) |>
      sum()
    score24 <-
      apply(sub24[,c(13,18)], 2, log) |>
      apply(1, prod) |>
      sum()
  }
  
  ## Get the delta ABC score
  deltaABC <- log(score24/score0)
  ## get the class 
  class <- na.omit(rna[rna$SYMBOL==gene,])$class
  if (length(class)==0){
    class <- NA
  }
  ## Make the row
  row <- data.frame(gene = gene,
                    class = class,
                    deltaABC = deltaABC)
  results <- rbind(results, row)
}

results <- na.omit(results)

## Density plots 
ggplot(results, aes(x = deltaABC, fill = class)) +
  geom_density(alpha = 0.8, color = NA) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_fill_manual(values = c("#F17D7E", "#FAA885", "gray")) +
  xlab("Differential genes at ABC anchors") +
  ylab("log(activatedABC/restingABC)")

