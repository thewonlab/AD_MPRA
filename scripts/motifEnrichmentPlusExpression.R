## Motif enrichment analysis for peaks that are activated or deactivated in microglia: 

###Initialize###
library(ggplot2)
library(ggforce)
library(cowplot)
library("rvest")

###Load this Function to plot Known Motifs and save the plots to files
plotKnownMotifs <- function(dirPath){
  
  #PHASE 1: read in each knownResults.txt file: 
  gained <- read.delim(paste0(dirPath, "gained/output/knownResults.txt"))
  lost <- read.delim(paste0(dirPath, "lost/output/knownResults.txt"))
  
  #PHASE2: put all motif files into a list
  motifs <- list(gained, lost)
  
  #PHASE3: clean up each motif file (new colnames, remove % signs)
  newcolnames <- c("Motif Name", "Consensus", "P.value", "Log P.value", "q value", 
                   "Number target sequences", "Percent_target_sequences", 
                   "Number background sequences", "Percent_background_sequences")
  
  abridge <- function(longname){
    temp <- strsplit(longname, '[(]')
    shortname <- unlist(temp)[1]
  }
  
  for (i in 1:length(motifs)){
    
    #Subset by the motif list
    motif <- motifs[[i]]
    #Assign new colnames
    colnames(motif) <- newcolnames
    #Assign abreviated names
    motif$`Motif Name` <- sapply(1:nrow(motif), function(x) abridge((as.character(motif[,1]))[x]))
    #Remove % signs from percent target sequences and percent background sequences
    motif$Percent_target_sequences <- as.numeric(sub("%", "", motif$Percent_target_sequences))
    motif$Percent_background_sequences <- as.numeric(sub("%", "", motif$Percent_background_sequences))
    #Calculate log2 enrichment and log10pval
    motif$log2enrichment <- log2(motif$Percent_target_sequences/motif$Percent_background_sequences)
    motif$log10pval <- -log10(exp(motif$`Log P.value`))
    #Set blank family column for later colorization
    motif$family <- ""
    
    #Reassign to the list 
    motifs[[i]] <- motif
    
  }
  
  #PHASE 4: Colorize motifs by their family 
  #Set motifs of interest (change this accordingly)
  ap1Names <- c("fos", "jun", "fra", "atf", "maf", "ap1", "ap-1")
  klfNames <- c("klf")
  runxNames <- c("runx")
  ctcfNames <- c("ctcf", "boris")
  gataNames <- c("gata")
  statNames <- c("stat")
  irfNames <- c("irf")
  
  family <- list(ap1Names, klfNames, runxNames, ctcfNames, gataNames, statNames, irfNames)
  familyLab <- c("AP1", "KLF", "RUNX", "CTCF", "GATA", "STAT", "IRF")
  
  for (x in 1:length(family)){
    for(i in 1:length(motifs)){
      motif <- motifs[[i]]
      name <- unlist(family[x])
      pos <- unique(sort(unlist(lapply(name, grep, motif$`Motif Name`, ignore.case = TRUE))))
      motif[pos,]$family <- familyLab[x]
      motifs[[i]] <- motif
    }
  }
  
  #PHASE 5: take care of any inf outliers in any of the motif sets 
  
  #Determine maximum for each log2enrichment column
  list <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    enrichment <- motif$log2enrichment[!(motif$log2enrichment==Inf)]
    maxE <- as.list(max(enrichment))
    maxE$i <- i
    list[[i]] <- maxE
  }
  maxE <- unlist(do.call(rbind,list)[,1])
  
  #Determine maximum for each log10pval column
  list <- list()
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    p <- motif$log10pval[!(motif$log10pval==Inf)]
    maxP <- as.list(max(p))
    maxP$i <- i
    list[[i]] <- maxP
  }
  maxP <- unlist(do.call(rbind,list)[,1])
  
  #Replace any Inf or -Inf with max + 1.1X max for enrichment and p value
  for (i in 1:length(motifs)){
    motif <- motifs[[i]]
    
    if(length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]>0)){
      motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)] <- rep((maxE[i]+(maxE[i]*1.1)), length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]))
    } else {
      print(paste0("There are no values of Inf for ", i))
    }
    
    if(length(motif$log2enrichment[motif$log2enrichment==-Inf])>0){
      motif$log2enrichment[motif$log2enrichment==-Inf] <- rep((-maxE[i]+(-maxE[i]*1.1)), length(motif$log2enrichment[motif$log2enrichment==-Inf]))
    } else {
      print(paste0("There are no values of -Inf for ", i))
    }
    
    if(length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]>0)){
      motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)] <- rep((maxP[i]+(maxP[i]*1.1)), length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]))
    } else {
      print(paste0("There are no values of Inf for ", i))
    }
    
    if(length(motif$log10pval[motif$log10pval==-Inf]>0)){
      motif$log10pval[motif$log10pval==-Inf] <- rep(0, length(motif$log10pval[motif$log10pval==-Inf]))
    } else {
      print(paste0("There are no values of -Inf for ", i))
    }
    
    motifs[[i]] <- motif
    
  }
  
}

## Results of Homer on differential peaks
gained <- read.delim("/work/users/m/a/marielle/work/AD3D/LIMA/atac/homer/gained/output/knownResults.txt")
lost <- read.delim("/work/users/m/a/marielle/work/AD3D/LIMA/atac/homer/lost/output/knownResults.txt")

motifs <- list(gained, lost)
names <- c("gained", "lost")

## Clean up each motif file (new colnames, remove % signs)
newcolnames <- c("MotifName", "Consensus", "P.value", "Log P.value", "q value", 
                 "Number target sequences", "Percent_target_sequences", 
                 "Number background sequences", "Percent_background_sequences")

abridge <- function(longname){
  temp <- strsplit(longname, '[(]')
  shortname <- unlist(temp)[1]
}

for (i in 1:length(motifs)){
  
  #Subset by the motif list
  motif <- motifs[[i]]
  #Assign new colnames
  colnames(motif) <- newcolnames
  #Assign abreviated names
  motif$MotifName <- sapply(1:nrow(motif), function(x) abridge((as.character(motif[,1]))[x]))
  #Remove % signs from percent target sequences and percent background sequences
  motif$Percent_target_sequences <- as.numeric(sub("%", "", motif$Percent_target_sequences))
  motif$Percent_background_sequences <- as.numeric(sub("%", "", motif$Percent_background_sequences))
  #Calculate log2 enrichment and log10pval
  motif$log2enrichment <- log2(motif$Percent_target_sequences/motif$Percent_background_sequences)
  motif$log10pval <- -log10(exp(motif$`Log P.value`))
  motif$class <- names[i]
  #Reassign to the list 
  motifs[[i]] <- motif
}

## Take care of any inf outliers in any of the motif sets 
#Determine maximum for each log2enrichment column
list <- list()
for (i in 1:length(motifs)){
  motif <- motifs[[i]]
  enrichment <- motif$log2enrichment[!(motif$log2enrichment==Inf)]
  maxE <- as.list(max(enrichment))
  maxE$i <- i
  list[[i]] <- maxE
}
maxE <- unlist(do.call(rbind,list)[,1])
#Determine maximum for each log10pval column
list <- list()
for (i in 1:length(motifs)){
  motif <- motifs[[i]]
  p <- motif$log10pval[!(motif$log10pval==Inf)]
  maxP <- as.list(max(p))
  maxP$i <- i
  list[[i]] <- maxP
}
maxP <- unlist(do.call(rbind,list)[,1])
#Replace any Inf or -Inf with max + 1.1X max for enrichment and p value
for (i in 1:length(motifs)){
  motif <- motifs[[i]]
  if(length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]>0)){
    motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)] <- rep((maxE[i]+(maxE[i]*1.1)), length(motif$log2enrichment[motif$log2enrichment==Inf | is.na(motif$log2enrichment)]))
  } else {
    print(paste0("There are no values of Inf for ", i))
  }
  
  if(length(motif$log2enrichment[motif$log2enrichment==-Inf])>0){
    motif$log2enrichment[motif$log2enrichment==-Inf] <- rep((-maxE[i]+(-maxE[i]*1.01)), length(motif$log2enrichment[motif$log2enrichment==-Inf]))
  } else {
    print(paste0("There are no values of -Inf for ", i))
  }
  
  if(length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]>0)){
    motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)] <- rep((maxP[i]+(maxP[i]*1.01)), length(motif$log10pval[motif$log10pval==Inf | is.na(motif$log10pval)]))
  } else {
    print(paste0("There are no values of Inf for ", i))
  }
  
  if(length(motif$log10pval[motif$log10pval==-Inf]>0)){
    motif$log10pval[motif$log10pval==-Inf] <- rep(0, length(motif$log10pval[motif$log10pval==-Inf]))
  } else {
    print(paste0("There are no values of -Inf for ", i))
  }
  motifs[[i]] <- motif
}


## Gained Motifs
gData <- motifs[[1]] 

## Lost Motifs 
lData <- motifs[[2]] 

## Rather than showing all of the top motifs per group, pick a few top ones and then show their family member RNA expression: 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt") #LIMA
rna$ENSEMBL <- rownames(rna)
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]

## Get RNA coordinate information: 
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
rna <- merge(genes, rna, by = "ENSEMBL")
rna <- na.omit(rna)

topLost <- 
  lData |> 
  head(10) |> 
  dplyr::pull(MotifName)

topGained <- 
  gData |> 
  head(10) |> 
  dplyr::pull(MotifName)


tf <- 
  list(AP1 = c("FOSL1", "FOSB", "FOSL2", "FOS", "JUN", "JUNB", "JUND"), 
       NFKB = c("NFKB2", "RELB", "NFKB1", "RELA", "REL"), 
       IRF = c("IRF2", "IRF3", "IRF4", "IRF7", "IRF8", "IRF9"),
       RUNX = c("RUNX1", "RUNX2", "RUNX3"), 
       ETS = c("ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "GABPA", "ERG1", "FLI1", "ERF", "ETV1", "ETV2", "ETV3", "ETV4", "ETV5", "ETV6", "ETV7", "EHF", "ETS1", "ETS2", "SPI1", "SPIB", "SPIC"), 
       ELK = c("ELK1", "ELK3", "ELK4"))

rna$family <- ""
for (i in 1:length(tf)){
  t <- tf[[i]]
  name <- names(tf)[i]
  ## Find the genes that match EXACTLY with the genes in the list, add the family name 
  rna[rna$SYMBOL %in% t,]$family <- name
}

## Pivot for timepoints to be separate: 
rnaSub <- 
  rna |> 
  dplyr::filter(!family=="")


rnaSub$SYMBOL <- factor(rnaSub$SYMBOL, levels = unlist(tf))
rnaSub$family <- factor(rnaSub$family, levels = c("ETS", "RUNX", "ELK", "AP1", "NFKB", "IRF"))

## Alternatively, try a lollipop style plot that combines LFC and RPKM 
rnaSub$baseMean <- (rnaSub$X0 + rnaSub$X24)/2


## Pick the highest three expressed for each category 
top3 <- 
  rnaSub |>
  dplyr::group_by(family) |>
  dplyr::arrange(desc(baseMean)) |> 
  dplyr::group_by(family) |> 
  dplyr::slice_max(order_by = baseMean, n = 3)
top3$x <- ""
top3$box <- ceiling(max(top3$baseMean)) #to add box around each 


## Flip this plot to make it short and wide 
tfRNAdots <- 
  ggplot(top3[top3$family %in% c("AP1", "NFKB", "IRF"),], aes(y = x, x = SYMBOL, color = log2FoldChange, size = baseMean)) + 
  geom_point(shape = 15) + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5), 
                        name = "RNA log2(fc)") + 
  scale_size_continuous(range = c(1, 8)) + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_classic() + 
  ylab("") + xlab("") + 
  guides(color = guide_legend(title.position = "top"), 
         size = guide_legend(title = "RNA counts",
                             override.aes = list(size = c(1,2,3)))) + 
  theme(legend.title = element_text(size = 8, color = "black"), 
        legend.text = element_text(size = 6, color = "black"),
        legend.position = "none",
        legend.box.margin = margin(4, 0, 0, 0, "cm"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        legend.key.size = unit(0.35, "cm"),
        aspect.ratio = 1, 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank()) +
  geom_point(aes(y = x, x = SYMBOL, size = box), shape = 0, color = "black") 

## To get the legends, make plots that only change one variable
colorOnly <- 
  ggplot(top3[top3$family %in% c("AP1", "NFKB", "IRF"),], aes(y = x, x = SYMBOL, color = log2FoldChange)) + 
  geom_point(shape = 15) + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-2.5, 2.5), 
                        name = "FC") + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_minimal() + 
  ylab("") + xlab("") + 
  theme(legend.title = element_text(size = 5, color = "black"), 
        legend.text = element_text(size = 4, color = "black"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        legend.key.size = unit(0.1, "cm"),
        aspect.ratio = 1) 
colLeg <- cowplot::get_legend(colorOnly)
plot_grid(colLeg, nrow = 1)

sizeOnly <- 
  ggplot(top3[top3$family %in% c("AP1", "NFKB", "IRF"),], aes(y = x, x = SYMBOL, size = baseMean)) + 
  geom_point(shape = 15) + 
  scale_size_continuous(range = c(1, 8)) + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_minimal() + 
  ylab("") + xlab("") + 
  guides(size = guide_legend(title = "RNA counts",
                             override.aes = list(size = c(0.5,1,1.5,2,2.5)))) + 
  theme(legend.title = element_text(size = 5, color = "black"), 
        legend.text = element_text(size = 4, color = "black"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        legend.key.size = unit(0.1, "cm"),
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        aspect.ratio = 1) 
sizeLeg <- cowplot::get_legend(sizeOnly)
plot_grid(sizeLeg, nrow = 1)

## Try splitting the motifs into two plots--one for gained and one for lost (to go to supplement)
tfRNAdotsLost <- 
  ggplot(top3[top3$family %in% c("ETS", "RUNX", "ELK"),], aes(y = x, x = SYMBOL, color = log2FoldChange, size = baseMean)) + 
  geom_point(shape = 15) + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5), 
                        name = "RNA log2(fc)") + 
  scale_size_continuous(range = c(1, 8)) + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_classic() + 
  ylab("") + xlab("") + 
  guides(color = guide_legend(title.position = "top"), 
         size = guide_legend(title = "RNA counts",
                             override.aes = list(size = c(1,2,3)))) + 
  theme(legend.title = element_text(size = 8, color = "black"), 
        legend.text = element_text(size = 6, color = "black"),
        legend.position = "none",
        legend.box.margin = margin(4, 0, 0, 0, "cm"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        legend.key.size = unit(0.35, "cm"),
        aspect.ratio = 1, 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background = element_blank()) +
  geom_point(aes(y = x, x = SYMBOL, size = box), shape = 0, color = "black") 

## To get the legends, make plots that only change one variable
colorOnlyLost <- 
  ggplot(top3[top3$family %in% c("ETS", "RUNX", "KLF"),], aes(y = x, x = SYMBOL, color = log2FoldChange)) + 
  geom_point(shape = 15) + 
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-2.5, 2.5), 
                        name = "FC") + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_minimal() + 
  ylab("") + xlab("") + 
  theme(legend.title = element_text(size = 5, color = "black"), 
        legend.text = element_text(size = 4, color = "black"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        legend.key.size = unit(0.1, "cm"),
        aspect.ratio = 1) 
colLegLost <- cowplot::get_legend(colorOnlyLost)
plot_grid(colLegLost, nrow = 1)

sizeOnlyLost <- 
  ggplot(top3[top3$family %in% c("ETS", "RUNX", "KLF"),], aes(y = x, x = SYMBOL, size = baseMean)) + 
  geom_point(shape = 15) + 
  scale_size_continuous(range = c(1, 8)) + 
  facet_row(vars(family), scales = "free_x", space = "free", strip.position = "bottom") + 
  theme_minimal() + 
  ylab("") + xlab("") + 
  guides(size = guide_legend(title = "RNA counts",
                             override.aes = list(size = c(0.5,1,1.5, 2, 2.5)))) + 
  theme(legend.title = element_text(size = 5, color = "black"), 
        legend.text = element_text(size = 4, color = "black"),
        axis.text.x = element_text(size = 5.5, vjust = 15), 
        panel.spacing = unit(0.1, "lines"), 
        legend.key.size = unit(0.1, "cm"),
        strip.text.x = element_text(size = 8, face = "bold", vjust = -2), 
        aspect.ratio = 1) 
sizeLegLost <- cowplot::get_legend(sizeOnlyLost)
plot_grid(sizeLegLost, nrow = 1)






