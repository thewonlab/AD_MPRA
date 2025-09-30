## Script to look at the relationship between ABC predicted variant-gene pairs and CRISPRi to generate figure 6d

## Load libraries
library(readxl)
library(data.table)
library(cowplot)
library(edgeR)

## Read in CRISPR data
crispr <- fread("/work/users/m/a/marielle/work/AD3D/CRISPRi/DESeq2_all_results_200kbCutoff_NewPadj.csv")
crispr$pair <- paste0(crispr$gRNA_name, "_", crispr$Gene)

## Read in ABC data
abc <- fread("/work/users/m/a/marielle/work/AD3D/LIMA/abc_sharedEnhancers/Neighborhoods/0h_EnhancerPredictions_Full.txt")

## Read in THP1 MPRA results 
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T, sep = "\t")
## Active variants
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
## Add activity: 
mpra$activity <- ""
mpra[mpra$id %in% active,]$activity <- TRUE
mpra[mpra$activity=="",]$activity <- FALSE
emvar <- mpra[mpra$sig==TRUE & mpra$activity==TRUE,]

## Convert to GRanges: 
mpraGR <- GPos(seqnames = Rle(mpra$seqnames), 
               pos = mpra$pos, 
               id = mpra$id, 
               rsid = mpra$rsid, 
               active = mpra$activity, 
               allelic = mpra$sig)
emvar <- mpraGR[mpraGR$active==TRUE & mpraGR$allelic==TRUE,]

## Intersect ABC enhancer anchor with emvars
abcGR <- 
  GRanges(seqnames = abc$seqnames1, 
          ranges = IRanges(start = abc$start1, end = abc$end1), 
          name = abc$anchor1.id, 
          abc = abc$ABCscore,
          gene = abc$anchor2.TargetGene)

abcOlap <- subsetByOverlaps(abcGR, emvar)
abcOlap$rsid <- emvar[nearest(abcOlap, emvar)]$rsid
abcOlap$pair <- paste0(abcOlap$rsid, "_", abcOlap$gene)
abcOlap <- as.data.frame(abcOlap)


## make a contingency table for 
### abc predicted vs not abc predicted
## upregulated vs downregulated
pairs <- unique(crispr$pair)
pairs <- pairs[grep("rs", pairs)]
df <- data.frame(pair = pairs)
df$abcPredicted <- FALSE
validabc <- abcOlap[abcOlap$abc>0.02,]
df[df$pair %in% validabc$pair,]$abcPredicted <- TRUE
df$downregulated <- FALSE
downreg <- crispr[crispr$log2FoldChange<0,]
df[df$pair %in% downreg$pair,]$downregulated <- TRUE

## Visualize it
df <- 
  data.frame(group = c("ABC Predicted", "ABC Predicted", "Not ABC Predicted", "Not ABC Predicted"), 
             direction = c("upregulated", "downregulated", "upregulated", "downregulated"), 
             frequency = c(7, 20, 79, 86))
df <- 
  df |> 
  group_by(group) |> 
  mutate(prop = frequency/sum(frequency))
df$prop <- round(df$prop, digits = 2) * 100

df$label <- paste0(df$frequency, " \n(", df$prop, "%)")

ggplot(df, aes(x = prop, y = group, fill = direction)) + 
  geom_bar(stat = "identity", position = "stack") + 
  geom_text(label = df$label, x = c(10, 60, 25, 75)) + 
  theme_minimal() + 
  scale_fill_manual(values = c("#EFB3AC", "#DB5849"))









