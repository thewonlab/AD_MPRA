## Script to identify target genes for emVars from Macrophage ABC and microglia eQTL: 
## Generates most panels in figure 4

## Load libraries
library(data.table)
library(ggplot2)
library(gprofiler2)

## Mpra allelic variants
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
## Mpra active elements 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
adjusted$active <- FALSE
adjusted[adjusted$id %in% active,]$active <- TRUE
emvar <- adjusted[adjusted$sig==TRUE & adjusted$active == TRUE,]$rsid

## Read in files
eqtl <- read.table("/work/users/m/a/marielle/work/AD3D/data/qtl/sig_BH_eqtls_kosoy_withsymbols.txt", header=TRUE)
abc <- read.table("/work/users/m/a/marielle/work/AD3D/data/abc/LIMA_ABC_enhOverlap.txt", header=TRUE)
## subset for those that overlap emvars 
eqtl <- eqtl[eqtl$rsid %in% emvar,]
abc <- abc[abc$rsid %in% emvar,]

## Get unique genes: 
eqtlGenes <- unique(eqtl$Gene_Symbol)
abcGenes <- unique(abc$anchor2.TargetGene)
allGenes <- unique(c(eqtlGenes, abcGenes)) |> na.omit()

## How many of these variants map to multiple genes? 
results <- data.frame()
for (i in 1:length(emvar)){
  variant <- emvar[i]
  genesA <- eqtl[eqtl$rsid %in% variant,]$Gene_Symbol
  genesB <- abc[abc$rsid %in% variant,]$anchor2.TargetGene
  targets <- unique(c(genesA, genesB))
  results <- rbind(results, 
                   data.frame(variant = variant, 
                              numGenes = length(targets)))
}
## remove the ones that do not have a target
results <- results[results$numGenes>0,]

## figure 4c right 
ggplot(results, aes(x = numGenes)) + 
  geom_histogram(fill = "#58C4B8", binwidth = 1, color = "white", lwd = 0.6) + 
  theme_minimal() + 
  xlab("# Genes per emVar") + 
  scale_x_continuous(breaks=seq(0,10)) + 
  ylab("Frequency") + 
  geom_vline(xintercept = median(results$numGenes)) + 
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  annotate("text", x = 4.25, y = 12, label = paste0("Mean = ", round(mean(results$numGenes), digits = 1)), size = 3, hjust = 0) + 
  annotate("text", x = 4.25, y = 11, label = paste0("Median = ", round(median(results$numGenes), digits = 1)), size = 3, hjust = 0)

## The inverse, how many variants per each emvar? 
results <- data.frame()
for (i in 1:length(allGenes)){
  gene <- allGenes[i]
  genesA <- eqtl[eqtl$Gene_Symbol %in% gene,]$rsid
  genesB <- abc[abc$anchor2.TargetGene %in% gene,]$rsid
  targets <- unique(c(genesA, genesB))
  results <- rbind(results, 
                   data.frame(gene = gene, 
                              numVar = length(targets)))
}
## remove the ones that do not have a target
results <- results[results$numVar>0,]

## figuere 4c left
ggplot(results, aes(x = numVar)) + 
  geom_histogram(fill = "#58C4B8", binwidth = 1, color = "white", lwd = 0.6) + 
  theme_minimal() + 
  xlab("# emVars per gene") + 
  scale_x_continuous(breaks=seq(0,10)) + 
  ylab("Frequency") + 
  geom_vline(xintercept = median(results$numVar)) + 
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  annotate("text", x = 2, y = 40, label = paste0("Mean = ", round(mean(results$numVar), digits = 1)), size = 3, hjust = 0) + 
  annotate("text", x = 2, y = 35, label = paste0("Median = ", round(median(results$numVar), digits = 1)), size = 3, hjust = 0)

## Create a pie chart for the 47 variants and whether they have abc support, qtl support, or both: 
abcVar <- unique(abc$rsid)
eqtlVar <- unique(eqtl$rsid)
shared <- eqtlVar[eqtlVar %in% abcVar]
abcOnly <- abcVar[!abcVar %in% shared]
eqtlOnly <- eqtlVar[!eqtlVar %in% shared]
any <- unique(c(shared, abcOnly, eqtlOnly))
neither <- emvar[!emvar %in% any]
classes <- list(abcOnly, eqtlOnly, shared, neither)

data <- 
  data.frame(class = c("abc", "eqtl", "both", "neither"), 
             vars = unlist(lapply(classes, length)))
data$class <- factor(data$class, levels = c("abc", "eqtl", "both", "neither"))
data$label <- paste0(c("ABC only:", "eQTL only:", "ABC and eQTL:", "Neither:"), "\n", data$vars, " variants")


## Visualize: 
## figure 4b
ggplot(data, aes(x="", y=vars, fill=class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  theme_void() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#7BDBAB", "#6ACDDE", "#58C4B8", "#D3DAE0")) + 
  geom_text(aes(y = c(45, 34, 15, 3), label = label), color = "white", size = 1.15, fontface="bold")

## Look at the gene expression of the genes that are uniquely identified by abc/qtl or both
## Read in data: 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt") #LIMA
rna$ENSEMBL <- rownames(rna)
## If RNA has transcripts: 
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]
rna$basemean <- (rna$X0+rna$X24)/2 #only do this for merged reps
## Get RNA coordinate information: 
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
rna <- merge(genes, rna, by = "ENSEMBL")
## Alternate gene level counts: 
counts <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaCounts.txt", header=T)
colnames(counts)[7] <- "ENSEMBL"
counts$ENSEMBL <- unlist(strsplit(counts$ENSEMBL, "[.]"))[seq(1,(2*nrow(counts)),by=2)]
rna <- merge(rna, counts, by = "ENSEMBL")

## Granges: 
rna <- GRanges(seqnames = Rle(rna$seqnames), 
               ranges = IRanges(start = rna$start, end = rna$end), 
               lfc = rna$log2FoldChange, 
               SYMBOL = rna$SYMBOL, 
               class = rna$class, 
               baseMean = rna$basemean, 
               padj = rna$padj, 
               counts0_1 = rna$LIMA_THP1_WT_0_3_1_1, 
               counts0_2 = rna$LIMA_THP1_WT_0_2_1_1,
               counts24_1 = rna$LIMA_THP1_WT_1440_3_1_1, 
               counts24_2 = rna$LIMA_THP1_WT_1440_2_1_1)
promoters <- promoters(rna)


## Determine distance to TSS for each group: 
abcDistance <- data.frame()
for (i in 1:length(unique(abc$rsid))){
  var <- unique(abc$rsid)[i]
  dist <- abs(abc[abc$rsid==var,]$start2 - abc[abc$rsid==var,]$pos)
  df <- data.frame(dist = dist, 
                   class = "abc")
  abcDistance <- rbind(abcDistance, df)
}

eqtlDistance <- data.frame()
for (i in 1:length(unique(eqtl$rsid))){
  var <- unique(eqtl$rsid)[i]
  genes <- eqtl[eqtl$rsid==var,]$Gene_Symbol |> unique()
  genes <- genes[!is.na(genes)]
  dist <- abs(unique(eqtl[eqtl$rsid==var,]$start) - start(rna[rna$SYMBOL %in% genes]))
  if(length(dist)==0){
    dist <- NA
  }
  df <- data.frame(dist = dist, 
                   class = "eqtl")
  eqtlDistance <- rbind(eqtlDistance, df)
}

## Combine: 
distances <- rbind(abcDistance, eqtlDistance) |> na.omit()
distances$dist <- distances$dist/1e3

## figure s6b
ggplot(distances, aes(x = dist)) + 
  geom_histogram(fill = "gray", alpha = 0.75, color = "black", lwd = 0.25) + 
  theme_minimal() + 
  ylab("Frequency") + xlab("Distance to TSS (kb)") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  geom_vline(xintercept = mean(distances$dist), lty = 2, color = "firebrick") + 
  annotate("text", label = paste0("Mean distance \n", mean(distances$dist), "Kb"), x = 450, y = 40, color = "firebrick", size = 2)

### USE GPROFILER TO DO GO/KEGG ENRICHMENT
eqtlGenes <- unique(eqtl$Gene_Symbol)
abcGenes <- unique(abc$anchor2.TargetGene)
allGenes <- unique(c(eqtlGenes, abcGenes))

genes <- allGenes ##76 genes

data <- data.frame()
sets <- c("GO", "KEGG") 
for (i in 1:length(sets)){
  goterm <- gost(query = genes, organism="hsapiens", ordered_query=F, significant=T, 
                 user_threshold=0.05, correction_method="fdr", sources=sets[i])
  goterm = goterm$result
  goterm = goterm[goterm$term_size<600 & goterm$intersection_size>4, ]
  goterm = goterm[!(goterm$term_id %in% goterm$parents),]
  goterm = goterm[order(goterm$p_value),]
  goterm2plot = goterm[1:10,c("term_name","p_value")]
  goterm2plot = goterm2plot[complete.cases(goterm2plot),]
  goterm2plot$log10P = -log10(goterm2plot$p_value)
  goterm2plot$term.name = factor(goterm2plot$term_name, levels=rev(goterm2plot$term_name))
  goterm2plot$class <- sets[i]
  data <- rbind(data, goterm2plot)
}
data$class <- factor(data$class, levels = c("GO", "KEGG"))
update_geom_defaults("text", list(size = 3))

## figure 5d
ggplot(filter(data, class=="GO"), aes(y = term.name, x = log10P, fill = class)) + 
  geom_bar(stat = "identity", alpha = 0.75) + 
  theme_minimal() + xlab("-log10(p-value)") + ylab("") + 
  geom_text(label = filter(data, class=="GO")$term.name, hjust = "left", x = 0.1, color = "black", size = 1.5) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "none") + 
  scale_fill_manual(values = c("#66CCC0"))

## figure s6a
ggplot(filter(data, class=="KEGG"), aes(y = term.name, x = log10P, fill = class)) + 
  geom_bar(stat = "identity", alpha = 0.75) + 
  theme_minimal() + xlab("-log10(p-value)") + ylab("") + 
  geom_text(label = filter(data, class=="KEGG")$term.name, hjust = "left", x = 0.1, color = "black", size = 1.5) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
        legend.position = "none") + 
  scale_fill_manual(values = c("#66CCC0"))




