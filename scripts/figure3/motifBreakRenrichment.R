## Motif Breaker downstream analaysis, putting it all together: 

## Load libraries
library(data.table)

## Get the list of files
files <- list.files("/work/users/m/a/marielle/work/AD3D/data/motifbreakr", "txt", full.names = TRUE)
backFiles <- list.files("/work/users/m/a/marielle/work/AD3D/data/motifbreakr/allBackground", "txt", full.names = TRUE)

## Combine all data into one dataframe
data <- data.frame()
for (i in 1:length(files)){
  data <- rbind(data, fread(files[i], header=TRUE))
}

background <- data.frame()
for (i in 1:length(backFiles)){
  print(i)
  background <- rbind(background, fread(backFiles[i], header=TRUE))
}

## Note to anyone poking around in the depths of this github (hi)
## We switched how we were identifying MPRA-active elements after this analysis, and instead of re-running it I just did this to shift the emvars that became inactive back to the background 
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
### LMER VERSION: 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
adjusted$active <- FALSE
adjusted[adjusted$id %in% active,]$active <- TRUE
emvar <- adjusted[adjusted$sig==TRUE & adjusted$active==TRUE,]$rsid
back <- adjusted[!adjusted$rsid %in% emvar,]$rsid
combined <- rbind(data, background)
data <- combined[combined$SNP_id %in% emvar,]
## save the data for table s4
write.table(data, file = "/work/users/m/a/marielle/work/AD3D/tables/Table_S4_motifbreakR_lmer.txt", row.names = F, sep = "\t", quote = F)
background <- combined[combined$SNP_id %in% back,]

## Figure out which TF motifs are present in either the regular files or the background files, we can only calculate p-values for those 
tf_a <- data$geneSymbol |> unique()
tf_b <- background$geneSymbol |> unique()

tf <- tf_a[tf_a %in% tf_b]

enrich <- data.frame()
for (i in 1:length(tf)){
  
  a <- length(unique(data[data$geneSymbol == tf[i],]$SNP_id))
  b <- (length(files))-a
  c <- length(unique(background[background$geneSymbol == tf[i],]$SNP_id))
  d <- (length(backFiles))-c
  
  fisher <- fisher.test(matrix(c(a,b,c,d), nrow = 2, byrow = T))
  
  df <- data.frame(tf = tf[i], 
                   pval = fisher$p.value, 
                   oddsratio = fisher$estimate, 
                   numvar = a)
  enrich <- rbind(enrich, df)
  
}

## Which TFs are enriched? 
enrich <- enrich[enrich$pval < 0.03,]
enrich$log10pval <- -log10(enrich$pval)
enrich <- enrich[order(enrich$log10pval, decreasing = TRUE),]
enrich$tf <- factor(enrich$tf, levels = enrich$tf)

## Get the corresponding RNA data
## Rather than showing all of the top motifs per group, pick a few top ones and then show their family member RNA expression: 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt") #LIMA
rna$ENSEMBL <- rownames(rna)
## If RNA has transcripts: 
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]
## Get RNA coordinate information: 
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
rna <- merge(genes, rna, by = "ENSEMBL")
rna <- na.omit(rna)
rna$baseMean <- (rna$X0+rna$X24)/2

## Get the genes of interest from the TF list: 
df <- data.frame()
remove <- c()
goi <- enrich$tf
for(i in 1:length(goi)){
  g <- as.character(goi[i])
  if(length(grep(":", g)>=1)){
    g <- unlist(strsplit(g, ":"))[1]
  }
  df <- rbind(df, rna[rna$SYMBOL %in% g,])
  ## save which didn't match up
  if(nrow(rna[rna$SYMBOL %in% g,])==0){
    remove <- c(remove, i)
  }
  
}

## Remove the TFs that do not have corresponding genes: 
enrich <- enrich[-remove,]

## Figure 3d, TF enrichments 
motifBreakRPlot <- 
  ggplot(enrich, aes(y = rev(tf), x = log10pval, size = oddsratio, fill = numvar)) + 
  geom_point(pch = 21, color = "black") + 
  theme_minimal() + 
  scale_size_continuous(range = c(1, 4)) + 
  xlab("") + ylab("-log10(P)") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 8, angle = 90, hjust = -0.25), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  scale_fill_gradient(high = "navy", low = "white")

## Add plot below showing the RNA expression of these TFs: 
for (i in 1:nrow(df)){
  x <- df$SYMBOL[i]
  if(nrow(df[df$SYMBOL==x,])>1){
    df$SYMBOL[i] <- paste0(df$SYMBOL[i], "_1")
  }
}
df$SYMBOL <- factor(df$SYMBOL, levels = df$SYMBOL)
df$y <- ""

## figure 3d, gene expression
expressionForMotifBreakR <- 
  ggplot(df, aes(y = rev(SYMBOL), x = y, color = log2FoldChange, size = baseMean)) + 
  geom_point(pch = 15) + 
  scale_size_continuous(range = c(1, 3)) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        legend.key.size = unit(0.4, "cm"), 
        legend.title = element_text(size = 10)) + 
  ylab("") + xlab("") + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6)) + 
  guides(size = "none")


expLegend <- 
  ggplot(df, aes(x = SYMBOL, y = y, color = log2FoldChange)) + 
  geom_point(pch = 15) + 
  theme_minimal() + 
  theme(aspect.ratio = 0.05, 
        axis.text.x = element_blank(), 
        legend.key.size = unit(0.2, "cm"), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6)) + 
  ylab("") + xlab("") + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-3, 6), name = "RNA\nlog(FC)") + 
  guides(size = "none")
expLegend <- legend <- cowplot::get_legend(expLegend)

colorOnly <- 
  ggplot(enrich, aes(x = tf, y = log10pval, color = numvar)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("") + ylab("-log10(P)") + 
  theme(axis.text.x = element_text(size = 6, angle = 90), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), 
        legend.key.size = unit(0.3, "cm"), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5)) + 
  scale_color_continuous(high = "navy", low = "gray", name = "Number of variants", guide = guide_colorbar(direction = "horizontal"))
colLeg <- legend <- cowplot::get_legend(colorOnly)
plot_grid(colLeg, nrow = 1)

sizeOnly <- 
  ggplot(enrich, aes(x = tf, y = log10pval, size = oddsratio)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("") + ylab("-log10(P)") + 
  theme(axis.text.x = element_text(size = 6, angle = 45), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6), 
        legend.key.size = unit(0.2, "cm"), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5), 
        legend.direction = "horizontal",
        aspect.ratio = 1) + 
  scale_size_continuous(range = c(1, 4))
sizeLeg <- legend <- cowplot::get_legend(sizeOnly)
plot_grid(sizeLeg, nrow = 1)
