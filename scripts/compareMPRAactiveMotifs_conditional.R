## Script to identify what TF motifs are enriched in MPRA-active that are active in either the resting or LPS/IFNy conditions
## Figures from this script: 1d

## Load libraries
library(data.table)

## Read in data: 
data <- fread("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_rescond.txt") ## Results from MPRA_active.R
## Further restrict to FDR < 0.05
data <- data[data$fdr_interaction<0.05,]

## Get a list of the variants that are significantly more active in each condition: 
resting <- data[data$condition=="resting",]$var |> unique()
treated <- data[data$condition=="treated",]$var |> unique()

## Regular MPRA data: 
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt") ## results from MPRA_allelic.R
mpra$class <- ""
mpra[mpra$id %in% resting,]$class <- "control"
mpra[mpra$id %in% treated,]$class <- "treated"
mpra[mpra$class=="",]$class <- "background"
mpra$id <- paste0("variant", 1:nrow(mpra))

## Prepare each set for Homer: 
cats <- unique(mpra$class)

homer <- list()
for (i in 1:length(cats)){
  data <- mpra[mpra$class==cats[i],]
  data <- data[,c(1,10,11,12)]
  data$start <- data$pos-75
  data$end <- data$pos+75
  data$score <- 1 #for bed file
  #data <- data[,c(1,2,5,6,4)] #homer file
  data <- data[,c(2,5,6,4,7,1)] #bed file
  #colnames(data) <- c("id", "chrom", "start", "end", "strand") #homer file
  colnames(data) <- c("chrom", "start", "end", "strand", "score", "name") #bed file
  homer[[i]] <- data
}
names(homer) <- cats

## Write outputs: 
write.table(homer$control, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_ttest/control/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer$treated, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_ttest/treated/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer$background, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_ttest/background.txt", row.names = F, sep = "\t", quote = F)

## I then move to the command line to launch homer with the following command from the directories where the background.txt files were writen to
## findMotifs.pl input/input.txt hg19 ./output -bg background.txt -size given -nomotif (you can run without nomotif if you want denovo motifs it just takes longer)

## After running homer, read in the output files for visualiziation 
## Read in data: 
control <- read.delim("/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM_FDR05/control/output/knownResults.txt")
treated <- read.delim("/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM_FDR05/treated/output/knownResults.txt")
motifs <- list(control, treated)
names <- c("control", "treated")

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

## Gained Motifs
data <- motifs[[1]] 
data <- data[1:5,]
data$MotifName <- factor(data$MotifName, levels = rev(data$MotifName))
motifsCon <- 
  ggplot(data, aes(x = MotifName, y = log10pval)) + 
  geom_bar(stat = "identity", fill = "#ECABA5", alpha = 0.75) + 
  theme_classic() + coord_flip() + 
  labs(title = "Resting") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 8)) + 
  xlab("") + ylab("-log10(p-val)") + 
  geom_text(label = data$MotifName, hjust = "left", y = 0.2, color = "white")


## Lost Motifs 
data <- motifs[[2]] 
data <- data[1:5,]
data$MotifName <- factor(data$MotifName, levels = rev(data$MotifName))
motifsTrt <- 
  ggplot(data, aes(x = MotifName, y = log10pval)) + 
  geom_bar(stat = "identity", fill = "#EB6F60", alpha = 0.75) + 
  theme_classic() + coord_flip() + 
  labs(title = "Activated") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 8)) + 
  xlab("") + ylab("-log10(p-val)") + 
  geom_text(label = data$MotifName, hjust = "left", y = 0.2, color = "white")

## Combine: (Figure 1d)
plot_grid(motifsCon, motifsTrt, nrow = 1)

## Save the plots themselves for later: 
save(motifsCon, motifsTrt, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/activityTFEnrich_lmer.rda")