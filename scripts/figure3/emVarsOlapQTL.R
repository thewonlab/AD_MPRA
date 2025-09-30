## Script to compare enrichment of HEK/THP1 emvars to microglia eQTLs from Kosoy 2022
## Generates figure 3g

## Load libraries
library(data.table)
library(dplyr)
library(GenomicRanges)
library(readxl)

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

## Comparison with eQTL 
eqtl = fread("/work/users/m/a/marielle/work/AD3D/data/qtl/sig_BH_eqtls_kosoy.txt")
mpra_eqtl = mpra %>%
  left_join(eqtl, by=c("rsid"="Variant"))
allelevec = (mpra_eqtl$a1==mpra_eqtl$`non-effect_allele`) + (mpra_eqtl$a1==mpra_eqtl$effect_allele) + (mpra_eqtl$a2==mpra_eqtl$`non-effect_allele`) + (mpra_eqtl$a2==mpra_eqtl$effect_allele)
table(allelevec) # All should be 2
mpra_eqtl$qtl_Beta_cor = ifelse(mpra_eqtl$a1==mpra_eqtl$`non-effect_allele`, mpra_eqtl$Beta, mpra_eqtl$Beta*-1)

emVar_eqtl = mpra_eqtl[mpra_eqtl$sig==TRUE & mpra_eqtl$activity==TRUE, ]
length(unique(emVar_eqtl$rsid)) # 47
emVar_eqtl_olap = emVar_eqtl[!is.na(emVar_eqtl$Gene), ]
length(unique(emVar_eqtl_olap$rsid)) # 32
length(unique(emVar_eqtl_olap$rsid))/length(unique(emVar_eqtl$rsid)) # 0.6808511
emVar_eqtl_IDE = emVar_eqtl_olap[emVar_eqtl_olap$logFC*emVar_eqtl_olap$qtl_Beta_cor>0,]
length(unique(emVar_eqtl_IDE$rsid)) # 27
length(unique(emVar_eqtl_IDE$rsid))/length(unique(emVar_eqtl_olap$rsid)) # 0.84375
length(unique(emVar_eqtl_IDE$rsid))/length(unique(emVar_eqtl$rsid)) # 0.5744681

## Compare with HEK cell MPRA data
hek_all = read_excel("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s1.xlsx", sheet = 2, skip=4) # allele
hek_act = read_excel("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s5.xlsx", sheet = 2, skip=11) # activity

hek_AD = hek_all[hek_all$Disorder=="AD",]
hek_act = hek_act[hek_act$label.fdr=="active",]

hekadrange = GRanges(hek_AD$chr, IRanges(hek_AD$pos, hek_AD$pos))
mcols(hekadrange) = hek_AD[,c(1:4,7:ncol(hek_AD))]

hekactrange = GRanges(hek_act$chr, IRanges((hek_act$start+hek_act$end)/2,(hek_act$start+hek_act$end)/2))

olap = findOverlaps(hekadrange, hekactrange)
hek_allact = hekadrange[unique(queryHits(olap))] 
hek_emvar = data.frame(hek_allact[hek_allact$q<0.05,]) # 76 emVars according to their standards

mpra_eqtl = hek_emvar %>%
  left_join(eqtl, by=c("rsID"="Variant"))
allelevec = (mpra_eqtl$A0==mpra_eqtl$`non-effect_allele`) + (mpra_eqtl$A0==mpra_eqtl$effect_allele) + (mpra_eqtl$A1==mpra_eqtl$`non-effect_allele`) + (mpra_eqtl$A1==mpra_eqtl$effect_allele)
table(allelevec) # All should be 2
mpra_eqtl$qtl_Beta_cor = ifelse(mpra_eqtl$A0==mpra_eqtl$`non-effect_allele`, mpra_eqtl$Beta, mpra_eqtl$Beta*-1) # MPRA logFC = alt/ref = A1/A0 vs. eQTL Beta = effect/non-effect

length(unique(mpra_eqtl$rsID)) # 76
emVar_eqtl_olap = mpra_eqtl[!is.na(mpra_eqtl$Gene), ]
length(unique(emVar_eqtl_olap$rsID)) # 28
length(unique(emVar_eqtl_olap$rsID))/length(unique(mpra_eqtl$rsID)) # 0.3684211
emVar_eqtl_IDE = emVar_eqtl_olap[emVar_eqtl_olap$`Log2.FC`*emVar_eqtl_olap$qtl_Beta_cor>0,]
length(unique(emVar_eqtl_IDE$rsID)) # 20
length(unique(emVar_eqtl_IDE$rsID))/length(unique(emVar_eqtl_olap$rsID)) # 0.7142857
length(unique(emVar_eqtl_IDE$rsID))/length(unique(mpra_eqtl$rsID)) # 0.2631579


## Combine into a dataframe 

df <- 
  data.frame(celltype = c("THP-1", "THP-1", "HEK", "HEK"), 
           factor = c("eQTL", "IDE", "eQTL", "IDE"), 
           value = c(0.6808511, 0.84375, 0.3684211, 0.7142857))
df$factor <- factor(df$factor, levels = c("IDE", "eQTL"))
df$celltype <- factor(df$celltype, levels = c("THP-1", "HEK"))

## to complete the stacked bar plot, make the inverse barplot 
inv <- df
inv$value <- 1-inv$value

df <- rbind(df, inv)
df$fill <- c(rep("olap", 4), rep("not", 4))
## adjut this so eqtl and idea have different colors: 
df$fill <- paste0(df$fill, "_", df$factor)
df$label <- round(df$value*100)
df[grep("not", df$fill),]$label <- ""

ggplot(df, aes(y = factor, x = value, fill = fill)) + 
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(~celltype, nrow = 2, strip.position = "left") + 
  theme_minimal() + 
  scale_fill_manual(values = c("#C2C2D1", "#C2C2D1", "#F2D18F", "#C6E3CA")) + 
  theme(legend.position = "none") + 
  xlab("") + ylab("") + 
  geom_text(aes(x = 0.1, label = label), position = position_stack(vjust = 0.5), na.rm = TRUE)


