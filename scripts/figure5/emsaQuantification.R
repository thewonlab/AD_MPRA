## Script to generate boxplots showing the quantification of MPRA binding
## This script uses the quantification for the emvar at the HLA locus and the 3 variants at the BIN1 locus and generates the figures shown in 5a, S8a,b,d,f

## Load libraries
library(readxl)
library(plotgardener)
library(ggtext)
library(tidyverse)

## Read in data for the HLA variant: 
hla <- read_xlsx("/work/users/m/a/marielle/work/AD3D/data/EMSA/HLA_EMSA_quant.xlsx")
hla$Group <- factor(hla$Group, levels = c("risk", "protect", "risk_competitor", "protect_competitor"))

## alternatively, do a binding score that uses the competitor to normalize? 
binding <- 
  hla |> 
  select(-`Raw Signal`) |> 
  na.omit()
binding$Probe <- factor(binding$Probe, levels = c("risk", "protect"))

## Definte statistical comparisons to be made
comp <- list(c("risk", "protect"))
## Visualize: 
ggplot(binding, aes(x = Probe, y = (`Binding score`), fill = Probe)) + 
  geom_boxplot() + 
  geom_point(shape = 21, size = 3) + 
  facet_wrap(~Condition, labeller = as_labeller(c("Control" = "<span style='color:#DD8D84;'>resting</span>", "Treated" = "<span style='color:#DE685B;'>LPS+IFNγ</span>"))) +
  theme_minimal() + 
  scale_x_discrete(labels = c("risk" = "<span style='color:#9653A0;'>A</span>", "protect" = "<span style='color:#B9ABD3;'>C</span>")) +
  theme(legend.position = "none", 
        strip.text = element_markdown(size = 12, face = "bold"), 
        axis.text.x = element_markdown(size = 12, face = "bold"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", size = 0.4)) + 
  xlab("") + 
  scale_fill_manual(values = c("#9653A0", "#B9ABD3")) + 
  ylab("EMSA Binding Score \n(Probe / Probe + Competitor)") + 
  stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 3)


## Variants tested at the bin1 locus: 
vars <- data.frame(chrom = "chr2", 
                   pos = c(127886157, 127892809, 127892955), 
                   snp = c("rs13025717", "rs6733839", "rs72838287"), 
                   p = c(0.5, 0.65, 0.4),
                   color = c("#366BFF", "#5532A0", "#6CB7AE"))

## bin1 quantification
data <- read_xlsx("/work/users/m/a/marielle/work/AD3D/data/EMSA/BIN1_Quant.xlsx") |> 
  select(-`Raw data`) |> 
  na.omit()
data$name <- paste0(data$Variant, "_", data$Allele)
data$name <- factor(data$name, levels = c("rs13025717_risk", "rs13025717_protect", "rs6733839_risk", "rs6733839_protect", "rs72838287_risk", "rs72838287_protect"))

## Plot each of the variants separately
comp <- list(c("rs13025717_risk", "rs13025717_protect"))
var1 <- 
  ggplot(data[data$Variant=="rs13025717",], aes(x = name, y = `Binding score`, fill = name)) + 
  geom_boxplot() +  
  geom_point() + 
  scale_fill_manual(values = c("#366BFF", "#94A5F8")) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 7)) + 
  xlab("") + 
  ylab("EMSA Binding Score\n(Probe / Probe + Competitor)") + 
  ylim(0,45) + 
  stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 3)

comp <- list(c("rs6733839_risk", "rs6733839_protect"))
var2 <- 
  ggplot(data[data$Variant=="rs6733839",], aes(x = name, y = `Binding score`, fill = name)) + 
  geom_boxplot() +  
  geom_point() + 
  scale_fill_manual(values = c("#5532A0", "#9582BE")) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 7)) + 
  xlab("") + 
  ylab("EMSA Binding Score\n(Probe / Probe + Competitor)") + 
  ylim(0,45) + 
  stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 3)

comp <- list(c("rs72838287_risk", "rs72838287_protect"))
var3 <- ggplot(data[data$Variant=="rs72838287",], aes(x = name, y = `Binding score`, fill = name)) + 
  geom_boxplot() +  
  geom_point() + 
  scale_fill_manual(values = c("#6CB7AE", "#BEDAD6")) + 
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 7)) + 
  xlab("") + 
  ylab("EMSA Binding Score\n(Probe / Probe + Competitor)") + 
  ylim(0,45) + 
  stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 3)

## make a schematic that shows where these variants fall in the BIN1 locus:
pageCreate(width = 5, height = 5, showGuides = FALSE)
genes <- plotGenes(chrom = "chr2", chromstart = 127847955, chromend = 127897955, assembly = "hg19", 
                   x = 0.5, y = 0.75, width = 3, height = 1, 
                   geneHighlights = data.frame(gene = "BIN1", 
                                               color = "#929090"))
annoGenomeLabel(plot = genes, 
                x = 0.5, y = 1.8)
plotManhattan(vars, chrom = "chr2", chromstart = 127847955, chromend = 127897955, assembly = "hg19", 
              x = 0.5, y = 0.5, height = 0.75, width = 3, 
              fill = colorby("snp", palette = colorRampPalette(vars$color)))
annoHighlight(genes, chrom = "chr2", chromstart = 127886157-150, chromend = 127886157+150, assembly = "hg19", 
              y = 0.75, height = 1.05, fill = "#366BFF")
annoHighlight(genes, chrom = "chr2", chromstart = 127892809-150, chromend = 127892809+150, assembly = "hg19", 
              y = 0.75, height = 1.05, fill = "#5532A0")
annoHighlight(genes, chrom = "chr2", chromstart = 127892955-150, chromend = 127892955+150, assembly = "hg19", 
              y = 0.75, height = 1.05, fill = "#6CB7AE")
## snp labels: 
plotText(label = "rs13025717\n(Cooper 2023)", fontcolor = "#366BFF", fontsize = 8, x = 2.45, y = 1.1)
plotText(label = "rs72838287 (Bond 2025)", fontcolor = "#6CB7AE", fontsize = 8, x = 3.85, y = 1.0)
plotText(label = "rs6733839 (Nott 2019)", fontcolor = "#5532A0", fontsize = 8, x = 3.85, y = 1.15)
## add the box plots 
plotGG(var1, x = 0.05, y = 2.15, width = 1.5, height = 3)
plotGG(var2, x = 1.75, y = 2.15, width = 1.5, height = 3)
plotGG(var3, x = 3.45, y = 2.15, width = 1.5, height = 3)
