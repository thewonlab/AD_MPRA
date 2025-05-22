## Script to compare AD emVar target genes to various datasets
## Generates figure 5e

## Load libraries: 
library(readxl)
library(ggforce)
library(cowplot)

## Read in data: 
targets <- read.csv("/work/users/m/a/marielle/work/AD3D/LIMA/targetGenes_qtlBH_correction_newactivelements.txt", header=F)$V1

## Background genes 
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
backgroundGenes <- rna[rna$basemean>10,]$SYMBOL #all

## Read in DAM data
dam <- read_xlsx("/work/users/m/a/marielle/external/hasselman2019/NIHMS1563396-supplement-Table_S7.xlsx", sheet = 4)
colnames(dam) <- 1:ncol(dam)
damUp <- 
  dam[8:nrow(dam),1] |> 
  as.data.frame() |> 
  na.omit()
damUp <- damUp$`1`
damDown <- 
  dam[8:nrow(dam),9] |> 
  as.data.frame() |> 
  na.omit()
damDown <- damDown$`9`
dam <- list(damUp, damDown)
names(dam) <- c("dam\nUp", "dam\nDown")

## Sun et al scRNA clusters
markers <- read_xlsx("/work/users/m/a/marielle/external/sun2023/1-s2.0-S0092867423009716-mmc1.xlsx", sheet = 3) |> 
  as.data.frame()
states <- unique(markers$microgliaState)
sun <- list()
for (i in 1:length(states)){
  mk <- states[i]
  sun[[i]] <- markers[markers$microgliaState==mk,]$gene
}
names(sun) <- states

## Gazestani et al response to abeta and tau
gazestani <- read_xlsx("/work/users/m/a/marielle/external/gazestani2023/1-s2.0-S0092867423008590-mmc4.xlsx", sheet = 5, skip = 7)
gazAbeta <- gazestani[,1:16]
gazTau <- gazestani[,19:34]
## significant
abeta <- gazAbeta[gazAbeta$adj.P.Val...6<0.05,]$gene...1
tau <- gazTau[gazTau$adj.P.Val...24<0.05,]$gene...19
gazestani <- list(abeta, tau)
names(gazestani) <- c("abeta", "tau+abeta")

## Mouse aging genes 
sex <- read.csv("/work/users/m/a/marielle/external/li2023/SEX_Human.txt", header=F)$V1

## Li with just the sex genes
li <- list(sex)
names(li) <- "aging"

## Olah et al scRNA clusters
markers <- read_xls("/work/users/m/a/marielle/external/olah2020/41467_2020_19737_MOESM7_ESM.xls", sheet = 1) |> 
  as.data.frame()
markers$up_type <- paste0("MG", markers$up_type)
states <- unique(markers$up_type)
olah <- list()
for (i in 1:length(states)){
  mk <- states[i]
  olah[[i]] <- markers[markers$up_type==mk,]$gene
}
names(olah) <- states

## Dolan et al scRNA clusters
markers <- fread("/work/users/m/a/marielle/external/dolan2023/41590_2023_1558_MOESM3_ESM.csv") |> 
  as.data.frame()
states <- unique(markers$cluster)
dolan <- list()
for (i in 1:length(states)){
  mk <- states[i]
  dolan[[i]] <- markers[markers$cluster==mk,]$gene
}
names(dolan) <- states

## Make a list of the different gene sets: 
sets <- list(dam, sun, gazestani, li, olah, dolan)
names(sets) <- c("Hassleman", "Sun", "Gazestani", "Li", "Olah", "Dolan")

enrichments <- data.frame()
overlapGenes <- list()
for (i in 1:length(sets)){
  paper <- sets[[i]]
  df <- data.frame()
  overlapGenes[[i]] <- list()
  for (x in 1:length(paper)){
    set <- unlist(paper[x])
    ## Contingency table: 
    a <- length(targets[targets %in% set])
    b <- length(targets[!targets %in% set])
    c <- length(backgroundGenes[backgroundGenes %in% set])
    d <- length(backgroundGenes[!backgroundGenes %in% set])
    
    ## Fishers test: 
    fisher <- fisher.test(matrix(c(a,b,c,d), nrow = 2, byrow = T))
    
    ## Row of data.frame
    row <- 
      data.frame(set = names(paper)[x], 
                 pvalue = fisher$p.value, 
                 or = fisher$estimate, 
                 overlap = a, 
                 setsize = length(set))
    df <- rbind(df, row)
    #save the actual overlap genes
    overlapGenes[[i]][[x]] <- targets[targets %in% set]
  }
  df$study <- names(sets)[i]
  enrichments <- rbind(enrichments, df)
}
enrichments$log10p <- -log10(enrichments$pvalue)
## Adjust pvalue
enrichments$padj <- p.adjust(enrichments$pvalue, method = "BH")
enrichments$log10p <- -log10(enrichments$padj)

# Calculate the number of unique x values for each facet
x_counts <- aggregate(set ~ study, enrichments, function(x) length(unique(x)))
total_x <- sum(x_counts$set)
x_counts$width <- x_counts$set / total_x

## Get the label to put under the points for set size
enrichments$label <- paste0("n = ", enrichments$setsize)

## Visualize the results! 
plots <- list()
size_scale <- scale_size_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 4), range = c(1, 5))
for (i in 1:length(sets)){
  x <- names(sets)[i]
  sub <- enrichments[enrichments$study==x,]
  if (i==1){
    pal <- c("#F3756F", "#F3756F")
    title <- paste(names(sets)[[1]], "et al 2019")
  } else if (i==2){
    sub$set <- factor(sub$set, levels = (c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG10", "MG11", "MG12")))
    pal <- (c("#EE2329", "#0094C8", "#24B355", "#BA65A8", "#F68A1F", "#FFEA16", "#995D2D", "#F391BC", "#AAAAA9", "#3FBDAC", "#F7986F", "#90B0D5"))
    title <- paste(names(sets)[[2]], "et al 2023")
  } else if (i==3){
    pal <- c("#6DB144", "#425EAB")
    title <- paste(names(sets)[[3]], "et al 2023")
  } else if (i==4){
    pal <- c("#F0AEB3")
    title <- paste(names(sets)[[4]], "et al 2023")
  } else if (i==5) {
    sub$set <- factor(sub$set, levels = (c("MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG9")))
    pal <- (c("#454EAF", "#97DCE9", "#9143B8", "#AB362D", "#F4E310", "#F0A205", "#1B6A31", "#71C403", "#DD2D15"))
    title <- paste(names(sets)[[5]], "et al 2020")
  } else if (i==6){
    sub$set <- factor(sub$set, levels = (c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)))
    pal <- (c("#B6B7B9", "#EC6C34", "#0F8D48", "#0F8D48", "#33479A", "#EEE430", "#0F8D48", "#EC6C34", "#EEE430", "#EEE430", "#B91C36"))
    title <- paste(names(sets)[[6]], "et al 2023")
  }
  
  if (i==1){
    plots[[i]] <- 
      ggplot(sub, aes(x = set, y = log10p, fill = set, size = overlap)) + 
      #geom_bar(stat = "identity", width = 0.008) +
      geom_point(pch = 21, color = "black") + 
      geom_hline(yintercept = 1.3, lty = 3) + 
      theme_minimal() + 
      xlab("") + ylab("-log10(Fisher adjusted p-value)") + 
      scale_fill_manual(values = pal) + 
      size_scale + 
      ggtitle(title) + 
      theme(plot.title = element_text(size = 6), 
            legend.position = "none", 
            axis.title.y = element_text(size = 8), 
            axis.text.x = element_text(size = 6), 
            plot.margin = unit(c(0, -0.1, 0, 0), "cm")) + 
      #annotate("text", label = sub$label, y = -1.25, x = 1:length(sub$label), size = 2) + 
      coord_cartesian(ylim = c(0, 10), clip = "off")
  } else if (i %in% c(2,3,4,5)){
    plots[[i]] <- 
      ggplot(sub, aes(x = set, y = log10p, fill = set, size = overlap)) + 
      #geom_bar(stat = "identity", width = 0.008) +
      geom_point(pch = 21, color = "black") + 
      geom_hline(yintercept = 1.3, lty = 3) + 
      theme_minimal() + 
      xlab("") + ylab("") + 
      scale_fill_manual(values = pal) + 
      size_scale + 
      ggtitle(title) + 
      theme(plot.title = element_text(size = 6), 
            legend.position = "none", 
            axis.text.y = element_blank(), 
            axis.text.x = element_text(size = 6), 
            plot.margin = unit(c(0, -0.1, 0, -0.1), "cm")) + 
      #annotate("text", label = sub$label, y = -1.25, x = 1:length(sub$label), size = 2) + 
      coord_cartesian(ylim = c(0, 10), clip = "off")
  } else if (i==6){
    plots[[i]] <- 
      ggplot(sub, aes(x = set, y = log10p, fill = set, size = overlap)) + 
      #geom_bar(stat = "identity", width = 0.008) +
      geom_point(pch = 21, color = "black") + 
      geom_hline(yintercept = 1.3, lty = 3) + 
      theme_minimal() + 
      xlab("") + ylab("") + 
      scale_fill_manual(values = pal) + 
      size_scale + 
      ggtitle(title) + 
      guides(fill="none") + 
      theme(plot.title = element_text(size = 6), 
            legend.key.size = unit(0.001, "cm"), 
            axis.text.y = element_blank(), 
            axis.text.x = element_text(size = 6), 
            legend.title = element_text(size = 6), 
            legend.text = element_text(size = 6), 
            plot.margin = unit(c(0, 0, 0, 0.5), "cm")) + 
      #annotate("text", label = sub$label, y = -1.25, x = 1:length(sub$label), size = 2) + 
      coord_cartesian(ylim = c(0, 10), clip = "off")
  }
  
}

## Use cowplot to figure out the correct sizing: 
targetGeneEnrich <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], rel_widths = c(3, 10, 3, 2, 10, 10), nrow = 1)

### SUPPLEMENTARY TABLE FOR THE GENE SETS THAT WERE USED: 
flat <- unlist(sets, recursive = FALSE)
maxLength <- max(sapply(flat, length))
df <- as.data.frame(lapply(flat, function(x) c(x, rep(NA, maxLength - length(x)))))
names <- colnames(df)

write.table(df, file = "/work/users/m/a/marielle/work/AD3D/tables/enrichmentGeneSets.txt", row.names = FALSE, sep = "\t", quote = F)

## Repeat this for the overlapping data: 
flat <- unlist(overlapGenes, recursive = FALSE)
maxLength <- max(sapply(flat, length))
df <- as.data.frame(lapply(flat, function(x) c(x, rep(NA, maxLength - length(x)))))
colnames(df) <- names

write.table(df, file = "/work/users/m/a/marielle/work/AD3D/tables/enrichmentOverlapGeneSets.txt", row.names = FALSE, sep = "\t", quote = F)

