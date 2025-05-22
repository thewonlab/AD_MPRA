## Script to identify what TF motifs are enriched in MPRA-active vs MPRA-inactive elements: 
## Figures from this script: 1c, S2c

## Load libraries: 
library(nullranges)
library(data.table)

## Read in data:

### LMER VERSION: 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t") ## This is the output from MPRA_active.R
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

## Regular MPRA data: 
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt") ## This is the output from MPRA_allelic.R
mpra$class <- ""
mpra[mpra$id %in% active,]$class <- "active"
mpra[mpra$class=="",]$class <- "background" ## aka inactive
mpra$id <- paste0("variant", 1:nrow(mpra))

## Get list of active variants: 
active <- filter(mpra, !class=="background")

## Get list of background variants: 
background <- filter(mpra, class=="background")

## Use matchranges to select 379 background variants on the same chromosomes as our active variants 
set.seed(100)
matched <- matchRanges(focal = filter(mpra, !class=="background"),
                       pool = filter(mpra, class=="background"),
                       covar = ~seqnames, method = "stratified", replace = FALSE) |> 
  as.data.frame()

## Background should be the vars that are NOT in either set
background <- mpra[!mpra$id %in% c(active$id, matched$id),]


## Prepare each set for Homer: 
cats <- list(active, matched, background)

homer <- list()
for (i in 1:length(cats)){
  data <- cats[[i]]
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
names(homer) <- c("active", "matched", "background")

## Write outputs: 
write.table(homer[[1]], "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/active_vs_inactive_GLM/active/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer[[2]], "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/active_vs_inactive_GLM/inactive/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer[[3]], "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/active_vs_inactive_GLM/background.txt", row.names = F, sep = "\t", quote = F)

## I then move to the command line to launch homer with the following command from the directories where the background.txt files were writen to
## findMotifs.pl input/input.txt hg19 ./output -bg background.txt -size given -nomotif (you can run without nomotif if you want denovo motifs it just takes longer)

## Read in homer results
active <- read.delim("/work/users/m/a/marielle/work/AD3D/data/mpra/homer/active_vs_inactive_GLM/active/output/knownResults.txt")
inactive <- read.delim("/work/users/m/a/marielle/work/AD3D/data/mpra/homer/active_vs_inactive_GLM/inactive/output/knownResults.txt")
motifs <- list(active, inactive)
names <- c("active", "inactive")

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

## Combine into one large dataframe: 
data <- do.call(rbind, motifs)
sigMotifs <- 
  data |> 
  filter(P.value < 0.05) |> 
  pull(class) |> 
  table() |> 
  as.data.frame()
colnames(sigMotifs) <- c("class", "numMotifs")

## Figure out which TFs are active in at least one condition and subset the data for them only 
tfs <- unique(data[data$P.value < 0.05,]$MotifName)
data <- data[data$MotifName %in% tfs,]
data$sig <- FALSE
data[data$P.value<0.05,]$sig <- TRUE
data$sig <- factor(data$sig, levels = c(TRUE, FALSE))
data$x <- "x"

## Figure S2c
tfBars <- 
  ggplot(sigMotifs, aes(x = class, y = numMotifs, fill = class))  + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  xlab("") + ylab("Number of Significantly\nEnriched Motifs") + 
  scale_fill_manual(values = c("#9FD141", "#CEE6A1")) + 
  geom_text(label = sigMotifs$numMotifs, vjust = -0.15, size = 3)

## Figure 1c
tfenrich_Density <- 
  ggplot(data, aes(x = log10pval, fill = class)) + 
  geom_density(alpha = 0.5, color = NA) + 
  theme_minimal() + 
  xlab("TF Motif Enrichment -log10(p-value)") + ylab("Density") +  
  scale_fill_manual(values =  c("#9FD141", "#CEE6A1")) + 
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  annotate("text", x = 9, y = 0.15, label = "Active Elements", size = 3, color = "#9FD141") + 
  annotate("text", x = 9, y = 0.5, label = "Inactive Elements", size = 3, color = "#CEE6A1")

save(tfBars, tfenrich_Density, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/activeInactiveTFenrich_lmer.rda")

