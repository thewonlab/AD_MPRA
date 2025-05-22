## Script to look at the GWAS p-values and betas for MPRA-active and MPRA-allelic variants (vs MPRA-inactive and MPRA-nonallelic respectively)
## Generates figures s2e-h 

## Script to look at the GWAS p-values and odds ratios of the allelic variants: 

## MPRA allelic data, generated through the script MPRA_allelic.R
adjusted <- 
  fread("/work/users/m/a/marielle/qtl/adMPRAresults.txt") |>
  mutate(start_pos = pos - 1) |> 
  as.data.frame()

## GWAS, Jansen 2019 summary statistics
gwas <- fread("/work/users/m/a/marielle/gwas/jansen/AD_sumstats_Jansenetal_2019sept.txt.gz")

### MPRA acitve data, generated from the script MPRA_active.R 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)
#reassess adjusted 
adjusted$activity <- "inactive"
adjusted[adjusted$id %in% active,]$activity <- "active"

## What are the GWAS p-vals for sig and nonsig allelic variants? 
pos <- filter(adjusted, allelic==TRUE) |> 
  pull(rsid)
neg <- filter(adjusted, allelic==FALSE) |> 
  pull(rsid)


active <- filter(adjusted, activity %in% c("shared", "control", "treated")) |> 
  pull(rsid)
active <- filter(adjusted, activity %in% c("active")) |> 
  pull(rsid)
inactive <- filter(adjusted, activity=="inactive") |> 
  pull(rsid)

posGWAS <- gwas[gwas$SNP %in% pos,]
posGWAS$class <- "MPRA-allelic"
negGWAS <- gwas[gwas$SNP %in% neg,]
negGWAS$class <- "MPRA-nonallelic"

actGWAS <- gwas[gwas$SNP %in% active,]
actGWAS$activity <- "MPRA-active"
inactGWAS <- gwas[gwas$SNP %in% inactive,]
inactGWAS$activity <- "MPRA-inactive"

data <- rbind(posGWAS, negGWAS)
dataAct <- rbind(actGWAS, inactGWAS)

data <- left_join(data, dataAct)
data$class <- factor(data$class, levels = c("MPRA-nonallelic", "MPRA-allelic"))

## Visualize
## figure s2e
gwasPboxplot <- 
  ggplot(data, aes(x = class, y = -log10(P), fill = class)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  ylim(0,25) + 
  scale_fill_manual(values = c("#ADD8E6", "#66ABF2")) + 
  theme_minimal() + 
  theme(legend.position="none", 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  xlab("") + 
  ylab("GWAS log10(p-value)")

## figure s2f
gwasBetaboxplot <- 
  ggplot(data, aes(x = class, y = BETA, fill = class)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  scale_fill_manual(values = c("#ADD8E6", "#66ABF2")) + 
  ylim(-0.1, 0.1) +
  theme_minimal() + 
  theme(legend.position="none", 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  xlab("") + 
  ylab("GWAS Beta")

## Any difference? no! 
wilcox.test(posGWAS$P, negGWAS$P)
wilcox.test(posGWAS$BETA, negGWAS$BETA)

## Now do activity: 
## Visualize
## figure s2g
data$activity <- factor(data$activity, levels = c("MPRA-inactive", "MPRA-active"))
gwasPboxplotActive <- 
  ggplot(data, aes(x = activity, y = -log10(P), fill = activity)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  ylim(0,25) + 
  scale_fill_manual(values = c("#CEE6A1", "#9FD141")) + 
  theme_minimal() + 
  theme(legend.position="none", 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  xlab("") + 
  ylab("GWAS log10(p-value)")

## figure s2h
gwasBetaboxplotActive <- 
  ggplot(data, aes(x = activity, y = BETA, fill = activity)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  scale_fill_manual(values = c("#CEE6A1", "#9FD141")) + 
  ylim(-0.1, 0.1) +
  theme_minimal() + 
  theme(legend.position="none", 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6)) + 
  xlab("") + 
  ylab("GWAS Beta")

## Any difference? no! 
wilcox.test(actGWAS$P, inactGWAS$P)
wilcox.test(actGWAS$BETA, inactGWAS$BETA)

