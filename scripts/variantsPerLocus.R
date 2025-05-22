## Script to show the number of variants at each GWAS locus that are acive and/or allelic 
## Generates figure 1g and s2b

## Load libraries
library(ggplot2)
library(ggpubr)
library(ggtext)
library(ggbreak)
library(GenomicRanges)
library(cowplot)

## Read in the rest of the MPRA data for variant info: ## results from mpra_allelic.R
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T, sep = "\t")
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t") ## results from mpra_active.R
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

## Add activity: 
adjusted$activity <- ""
adjusted[adjusted$id %in% active,]$activity <- TRUE
adjusted[adjusted$activity=="",]$activity <- FALSE

## Convert to GRanges: 
mpraGR <- GPos(seqnames = Rle(adjusted$seqnames), 
               pos = adjusted$pos, 
               id = adjusted$id, 
               rsid = adjusted$rsid, 
               active = adjusted$activity, 
               allelic = adjusted$sig)

## Generate venn Diagram from fig s2b
data <- as.data.frame(mpraGR)
data <- data[,c(5:7)]
data <- data[data$active==TRUE | data$allelic==TRUE,]
data$active <- as.logical(data$active)
putativeVenn <- 
  ggplot(data, aes(A = active, B = allelic)) + 
  geom_venn(auto_scale = TRUE, fill_color = c("#9BCC3F", "#66ABF2"), fill_alpha = 0.8, stroke_size = 0, text_size = 2, set_name_size = 3, show_percentage=FALSE) + 
  theme_void() + 
  coord_fixed()


## Overlap with GWAS loci 
loci <- read_xlsx("/work/users/m/a/marielle/work/MPRA/align/jansenLonelySNPs.xlsx") |> 
  as.data.frame() ## downloaded from Jansen 2019 
loci <- loci[!is.na(loci$pos),] #these are the 29 loci
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
loci <- sort(loci) ## should be the 29 loci, regardless of if they're lonely SNPs or not 

## Make the seqlevels actually in chr1:22 order
seqlevels(loci) <- paste0("chr", 1:22)
loci <- sort(loci)

## Add an id column that has the chromosome and start position 
loci$id <- paste0(seqnames(loci), ":", start(loci))
## And a numerical locus 1:29 as well 
loci$locus <- as.factor(1:29)

## Assign each variant to the closest GWAS loci: 
mpraGR$locus <- nearest(mpraGR, loci)
mpraGR$locus <- as.numeric(mpraGR$locus)
mpraGR <- mpraGR[-which(is.na(mpraGR$locus))]
mpraGR$emvar <- FALSE
mpraGR[mpraGR$active==TRUE & mpraGR$allelic==TRUE,]$emvar <- TRUE
mpraGR$active <- as.logical(mpraGR$active)
mpraGR$emvar <- as.logical(mpraGR$emvar)

## Categorize each variant as to whether it's only active, only allelic, or both: 
mpraGR$group <- case_when(
  mpraGR$emvar == TRUE ~ "emvar",
  mpraGR$active == TRUE & mpraGR$allelic == FALSE ~ "active",
  mpraGR$active == FALSE & mpraGR$allelic == TRUE ~ "allelic",
  TRUE ~ "inactive"
)

## Then for each locus, figure out if it has an active element, allelic variant, or both
results <- data.frame()
for (i in 1:length(loci)){
  vars <- mpraGR[mpraGR$locus==i]
  df <- 
    data.frame(locus = i, 
               rsid = loci$rsid[i], 
               vars = length(vars), 
               chr = as.character(seqnames(loci)[i]),
               active = sum(vars$active), 
               allelic = sum(vars$allelic), 
               emvar = sum(vars$emvar), 
               numActive = length(vars[vars$group=="active"]), 
               numAllelic = length(vars[vars$group=="allelic"]), 
               numEmVar = length(vars[vars$group=="emvar"]))
  results <- rbind(results, df)
}

results$locus <- factor(results$locus, levels = 1:29)
results$chr <- factor(results$chr, levels = paste0("chr", c(1,2,4,6,7,8,10,11,14,15,16,17,18,19,20)))
results <- results[order(results$chr),]
results <- results[order(results$chr),]

## Add back the identifier column from the locus GR object
results <- 
  left_join(results, as.data.frame(loci), by = "locus") |> 
  dplyr::select(-c(seqnames, start, end, width, strand, rsid.y))

## QC Checks--mark which loci were not tested in MPRA
results$tested <- "tested"
results[results$vars==0,]$tested <- "untested"
colnames(results)[2] <- "rsid"

## Reformat the data to make a barplot per locus and group them by which loci have how many active/allelic/emvar
data <- 
  results |> 
  pivot_longer(-c(locus, rsid, vars, chr, allelic, active, emvar, id, tested), names_to = "feature", values_to = "value") |> 
  mutate(feature = factor(feature, levels = c("numEmVar", "numAllelic", "numActive")))

## Top part of figure 1g
lociBarplot <- 
  ggplot(data, aes(x = locus, y = value, fill = feature)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme_minimal() + 
  ylab("Number of variants") + xlab("") + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none") + 
  scale_fill_manual(values = c("#158C7E", "#66ABF2", "#9BCC3F"))

## Make a heatmap to annotate the barplot Z
## Annotate below: 
annotate <- results
colnames(annotate)[7]
annotate[annotate$active>0,]$active <- 1
annotate[annotate$allelic>0,]$allelic <- 1
annotate[annotate$emvar>0,]$emvar <- 1
## Add blank columns for the id 
colnames(annotate)[11] <- "locusName"
annotate$id <- 0
## mutate
annotate <- 
  annotate |> 
  dplyr::select(-c(rsid, vars, chr, numActive, numAllelic, numEmVar, tested)) |> 
  mutate(across(everything(), as.factor)) |> 
  pivot_longer(cols = -"locusName", names_to = "feature", values_to = "overlap")
annotate$overlap <- as.character(annotate$overlap)
annotate$feature <- factor(annotate$feature, levels = c("id", "locus", "emvar", "allelic", "active"))
annotate$locusName <- factor(annotate$locusName, levels = results$id)

## fill with specific colors for groups: 
annotate$fill <- "#D3DAE0"
annotate[annotate$feature=="active" & annotate$overlap==1,]$fill <- "#9BCC3F"
annotate[annotate$feature=="allelic" & annotate$overlap==1,]$fill <- "#66ABF2"
annotate[annotate$feature=="emvar" & annotate$overlap==1,]$fill <- "#158C7E"
annotate[annotate$feature=="locus",]$fill <- "white"
annotate[annotate$feature=="id",]$fill <- "white"

annoPlot <- 
  ggplot(annotate, aes(x = locusName, y = feature, fill = fill)) + 
  geom_tile(color = "white", linewidth = 0.5) +  
  rotate_x_text(angle = 45) + 
  scale_fill_manual(values = c("#158C7E", "#66ABF2", "#9BCC3F", "#D3DAE0", "white", "white")) + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("") + ylab("") + 
  theme(legend.position = "none",
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.line = element_blank(), 
        panel.background = element_rect(fill="transparent"), 
        aspect.ratio = 0.1, 
        axis.text.x = element_blank()) + 
  geom_text(data = annotate[annotate$feature=="locus",], label = annotate[annotate$feature=="locus",]$overlap, size = 2.5) + 
  geom_text(data = annotate[annotate$feature=="id",], label = annotate[annotate$feature=="id",]$locusName, size = 2.5, angle = 90, hjust = 0.95) + 
  coord_cartesian(clip = "off")

## combine the plots
plot_grid(lociBarplot, annoPlot, ncol = 1, align = "v")

