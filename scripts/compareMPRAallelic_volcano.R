## Script to generate a volcano plot to visualize the difference in allelic activity
## Generates figure 1f

## Load libraries
library(data.table)
library(ggplot2)

## Read in data
### MPRA-active elements: 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

## MPRA-allelic elements: 
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt")

## Assign their activity 
mpra$class <- ""
mpra[mpra$id %in% active,]$class <- "active"
mpra[mpra$class=="",]$class <- "inactive"

## If the variant is active in at least one condition, highlight as putative causal
mpra$emvar <- ""
mpra[mpra$sig==TRUE & !mpra$class=="inactive",]$emvar <- TRUE
mpra[mpra$emvar=="",]$emvar <- FALSE
mpra$emvar <- factor(mpra$emvar, levels = c(TRUE, FALSE))

## Make volcano plot: 
mpraVolcano <- 
  ggplot(mpra, aes(x = corrected_logFC, y = -log10(adj.P.Val))) + 
  geom_point(data = mpra[mpra$sig==FALSE,], inherit.aes = TRUE, color = "#D3DAE0", alpha = 0.6) + 
  geom_point(data = mpra[mpra$sig==TRUE,], inherit.aes = TRUE, color = "#66ABF2", alpha = 0.6) + 
  geom_point(data = mpra[mpra$emvar==TRUE,], inherit.aes = TRUE, color = "#158C7E") + 
  geom_hline(yintercept = -log10(0.05), lty = 3) + 
  theme_minimal() +
  xlim(c(-2.5,2.5)) +
  ylab("-log10(padj)") + xlab("log2(risk/protective)") + 
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10))

## Create barplot showing the number of risk/protective upregulated variants to put below: 
sig <- mpra[mpra$sig==TRUE,]
sig[sig$corrected_logFC>0,]$class <- "risk"
sig[sig$corrected_logFC<0,]$class <- "protective"
sig$group <- ""

data <- table(sig$class, sig$emvar) |> 
  as.data.frame()
colnames(data) <- c("class", "emvar", "freq")
data$group <- ""
data$label <- data$freq
data[data$class=="protective",]$freq <- data[data$class=="protective",]$freq * -1


## Make two separate plots instead of one that combines them 
riskBar <- 
  ggplot(data, aes(x = group, fill = emvar)) + 
  geom_bar(data = data[data$class=="risk",], aes(y = freq), stat = "identity", position = "stack", width = 0.4) + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.line = element_blank(), 
        panel.background = element_rect(fill="transparent"), 
        plot.margin = unit(c(-1,-1,-1,-1), "cm")) + 
  xlab("") + ylab("") + 
  scale_fill_manual(values = c("#158C7E", "#66ABF2")) +
  geom_text(data = data[data$class=="risk",], aes(label = rev(label), y = label), vjust = 0.5, hjust = -1.5, color = "white", size = 2) + 
  coord_flip()

protBar <- 
  ggplot(data, aes(x = group, fill = emvar)) + 
  geom_bar(data = data[data$class=="protective",], aes(y = freq), stat = "identity", position = "stack", width = 0.4) + 
  theme_classic() + 
  coord_flip() + 
  scale_fill_manual(values = c("#158C7E", "#66ABF2")) + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.line = element_blank(), 
        panel.background = element_rect(fill="transparent"), 
        plot.margin = unit(c(-1,-1,-1,-1), "cm")) + 
  xlab("") + ylab("") + 
  geom_text(data = data[data$class=="protective",], aes(label = rev(label), y = freq), vjust = 0.5, hjust = 3, color = "white", size = 2)

## Save plots: 
save(mpraVolcano, riskBar, protBar, file = "/work/users/m/a/marielle/work/AD3D/savedPlotsforFigs/mpraVolcanowithBars_lmer.rda")


# ggplot(data, aes(x = group, fill = hit)) + 
#   geom_bar(data = data[data$class=="risk",], aes(y = freq), stat = "identity", position = "stack", color = "white", width = 0.4) +
#   geom_bar(data = data[data$class=="protective",], aes(y = -freq), stat = "identity", position = "stack", color = "white", width = 0.4) +
#   theme_classic() + 
#   coord_flip() + 
#   theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank(), 
#         axis.ticks.y = element_blank(), axis.ticks.x = element_blank()) + 
#   xlab("") + ylab("") + 
#   scale_fill_manual(values = c("#227B7F", "gray70"))




