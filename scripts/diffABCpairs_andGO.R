## Script to perform downstream ABC analysis 
## Generates figures 2g and 2i 

library(ggvenn)
library(GenomicRanges)
library(InteractionSet)
library(mariner)
library(data.table)
library(ggplot2)
library(VennDiagram)
library(rrvgo)
library(cowplot)
library(ggpubr)
library(ggsignif)
library(gprofiler2)

## ABC FILES LIMA: 
## Read in ABC files, this is available as supplemental information 
abc <- list.files("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/Neighborhoods", "Predictions.txt", full.names = TRUE)[1:2] |> 
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

## Separate 
abc0 <- abc[[1]]
abc24 <- abc[[2]]

### RUN THIS FOR THE SHARED ABC 
shared0 <- abc0[abc0$pair %in% abc24$pair]
unique0 <- abc0[!abc0$pair %in% abc24$pair]

shared24 <- abc24[abc24$pair %in% abc0$pair]
unique24 <- abc24[!abc24$pair %in% abc0$pair]

resting <- unique0
activated <- unique24
shared <- shared0
##total pairs
pairs <- unique(c(resting$pair, activated$pair, shared$pair))

## Vennn diagram to show pairs 
venndata <- data.frame(pairs)
venndata$resting <- ""
venndata$activated <- ""
venndata[venndata$pairs %in% resting$pair,]$resting <- TRUE
venndata[venndata$pairs %in% shared$pair,]$resting <- TRUE
venndata[venndata$pairs %in% activated$pair,]$activated <- TRUE
venndata[venndata$pairs %in% shared$pair,]$activated <- TRUE
venndata[venndata$resting=="",]$resting <- FALSE
venndata[venndata$activated=="",]$activated <- FALSE
venndata$resting <- as.logical(venndata$resting)
venndata$activated <- as.logical(venndata$activated)

## Figure 2g
abcVenn <- 
  ggplot(venndata, aes(A = resting, B = activated)) + 
  geom_venn(auto_scale = TRUE, fill_color = c("#F0B0AA", "#EC6F61"), fill_alpha = 0.8, stroke_size = 0, text_size = 2, set_name_size = 3, show_percentage=FALSE) + 
  theme_void() + 
  coord_fixed()

##### GENE EXPRESSION #####

## Which genes are unique to each timepoint? 
genes0 <- unique0[unique0$anchor2.TargetGeneExpression > 0]$anchor2.TargetGene |> unique()
genes24 <- unique24[unique24$anchor2.TargetGeneExpression > 0]$anchor2.TargetGene |> unique()

## Shared genes: 
shared <- genes0[genes0 %in% genes24]
genesUniq0 <- genes0[!genes0 %in% genes24]
genesUniq24 <- genes24[!genes24 %in% genes0]
sharedOnly <- shared[!shared %in% c(genesUniq0, genesUniq24)]
background <- unique(c(genesUniq0, sharedOnly, genesUniq24))
#background <- unique(c(abc0$anchor2.TargetGene, abc24$anchor2.TargetGene))

## Try subsetting for the genes that are DIFFERENTIAL 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt") #Differential gene expression information from Reed et al 2022
rna$ENSEMBL <- rownames(rna)
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
rna <- merge(genes, rna, by = "ENSEMBL")
upGenes <- rna[rna$class=="gained",]$SYMBOL
downGenes <- rna[rna$class=="lost",]$SYMBOL
staticGenes <- rna[rna$class=="static",]$SYMBOL
deg <- c(upGenes, downGenes)

## Differential genes that are uniquely at each condition of ABC pairs
genesUniq0 <- genesUniq0[genesUniq0 %in% deg] 
genesUniq24 <- genesUniq24[genesUniq24 %in% deg]

## Perform GO enrichment with gprofiler on deg unique to each category 
cats <- c("resting", "activated")
sets <- list(genesUniq0, genesUniq24)

godata <- data.frame()
for (i in 1:length(cats)){
  g <- sets[[i]]
  goterm <- gost(query = g, organism="hsapiens", ordered_query=F, significant=F,
                 user_threshold=0.05, correction_method="fdr", sources=c("GO"), 
                 custom_bg = background) # running GO
  goterm = goterm$result
  goterm = goterm[goterm$term_size<600 & goterm$intersection_size>4, ]
  goterm = goterm[!(goterm$term_id %in% goterm$parents),]
  goterm = goterm[order(goterm$p_value),]
  goterm2plot = goterm[1:10,c("term_name","p_value")]
  goterm2plot = goterm2plot[complete.cases(goterm2plot),]
  goterm2plot$log10P = -log10(goterm2plot$p_value)
  goterm2plot$term.name = factor(goterm2plot$term_name, levels=rev(goterm2plot$term_name))
  goterm2plot$class <- cats[[i]]
  godata <- rbind(godata, goterm2plot)
}
godata$class <- factor(godata$class, levels = c("resting", "activated"))
abcExpressedGO <- 
  ggplot(godata, aes(y = term.name, x = abs(log10P), fill = class, group = class)) + 
  geom_bar(stat = "identity", alpha = 0.85, width = 0.5) + 
  theme_minimal() +
  scale_fill_manual(values = c("#FAC5C1", "#F7948A")) + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        aspect.ratio = 1) + 
  ylab("GO Terms enriched in\nresting/activated ABC pairs") + xlab("-log10(p-val)") + 
  geom_text(label = godata$term.name, hjust = "left", x = 0.1, color = "black", size = 2)





