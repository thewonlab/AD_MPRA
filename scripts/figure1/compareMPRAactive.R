## Script to visualize the density of mpra active/inactive elements for N=9 pool vs negative controls and positive controls: 
## Figures from this script: 1b, S2a 

## Load libraries: 
library(ggvenn)


## Read in data: ###

### This is available as processed data at GEO GSE273887 (AD_MPRA_aggregated_counts.txt)
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batchvar_N9N6_combined_withControls.txt", header=T, sep = "\t")

## MPRA-active elements: 

### LMER VERSION: 
## Read in the activity output from the MPRA_active_lmer.R script
## Specifically we want the res_var object that describes activity at the variant level
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t") ## generated from MPRA_active.R script
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)


## Add variant category to batchvar: 
batchvar$class <- ""
batchvar[grep("chr", batchvar$SNP),]$class <- "tested"
batchvar[grep("Scrambled", batchvar$SNP),]$class <- "negative"
batchvar[batchvar$class=="",]$class <- "positive"

## How many elements for each class? 
table(batchvar$class) ## Note: this is shown in Fig S1b

## For the batchvar, get the control and treated data separately; 
ratios <- batchvar[seq(5, (ncol(batchvar)), by = 4)]
control <- apply(ratios[,1:9], 1, median)
treated <- apply(ratios[,10:15], 1, median)

## Add variant info back on to control and treated: 
control <- cbind(batchvar[,c(1,62:64)], control)
colnames(control) <- c("variant", "SNP", "allele", "class", "lfc")
control <- na.omit(control)
control$lfc <- scale(control$lfc, center = TRUE, scale = FALSE)
control$condition <- "control"

treated <- cbind(batchvar[,c(1,62:64)], treated)
colnames(treated) <- c("variant", "SNP", "allele", "class", "lfc")
treated <- na.omit(treated)
treated$lfc <- scale(treated$lfc, center = TRUE, scale = FALSE)
treated$condition <- "treated"

## Change tested to be either active or inactive, depending on whether it's in the active variants or not: 
control[control$SNP %in% active,]$class <- "active"
control[control$class=="tested",]$class <- "inactive"

treated[treated$SNP %in% active,]$class <- "active"
treated[treated$class=="tested",]$class <- "inactive"

## Combine into one dataframe: 
merged <- rbind(control, treated)

## VISUALIZATIONS ###

## For violin plots, use the following code: 
## Identify which positive controls are not Cmv/e1f and remove them 
merged <- merged[-c(grep("BAR", merged$variant), grep("TK", merged$variant)),] ## Note, if you're interested these were included in our construct but were for another project, so we remove them here
merged$class <- factor(merged$class, levels = c("negative", "inactive", "active", "positive"))


## Figure 1b
activityViolin <- 
  ggplot(merged, aes(x = class, y = lfc, fill = class)) + 
  geom_violin(alpha = 0.8) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, alpha = 0.25) + 
  facet_wrap(~condition, strip.position = "bottom") + 
  theme_minimal() + 
  scale_fill_manual(values = c("#D3DAE0", "#A1BD6C", "#7AAB1C", "#EBBC2E")) + 
  xlab("") + 
  ylab("log(RNA/DNA)") + 
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 5.5, face = "bold", vjust = 20))
wilcox.test(merged[merged$class=="active" & merged$condition=="control",]$lfc, merged[merged$class=="inactive" & merged$condition=="control",]$lfc)$p.value
wilcox.test(merged[merged$class=="active" & merged$condition=="treated",]$lfc, merged[merged$class=="inactive" & merged$condition=="treated",]$lfc)$p.value


## Another density plot that combines active/inactive to show theres no difference overall:
merged$class2 <- as.character(merged$class)
merged[merged$class2=="active" | merged$class2=="inactive",]$class2 <- "AD"
merged$class2 <- factor(merged$class2, levels = c("negative", "AD", "positive"))

## Figure S2a
densitySupplement <- ggplot(merged, aes(x = class2, y = lfc, fill = class2)) + 
  geom_violin(alpha = 0.8) + 
  geom_boxplot(outlier.shape = NA, width = 0.1, alpha = 0.25) + 
  facet_wrap(~condition, strip.position = "bottom") + 
  theme_minimal() + 
  scale_fill_manual(values = c("#D3DAE0", "#8EC427", "#FFAC31")) + 
  xlab("") + ylab("log(RNA/DNA)") + 
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 5.5, face = "bold", vjust = 20), 
        axis.text.x = element_text(size = 6), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 6))
wilcox.test(merged[merged$class2=="AD",]$lfc, merged[merged$class2=="negative",]$lfc)

save(activityViolin, densitySupplement, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/activityViolin_lmer.rda")



