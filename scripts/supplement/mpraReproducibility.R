## Full reproducibility analysis for supplement: 
## Generates figures s1d-e (reproducibility heatmaps)

## Load libraries
library(ggcorrplot)
library(corrplot)
library(RColorBrewer)

## Load the data: 
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch_N9N6_combined_batchvar.txt", header=T) ## Available from GEO: GSE273887_AD_MPRA_aggregated_counts.txt.gz

## Separate the control and treated samples: 
con <- batchvar[grep("con", colnames(batchvar))]
con <- con[,seq(3, ncol(con), by = 4)]
con <- na.omit(con)

treat <- batchvar[grep("treat", colnames(batchvar))]
treat <- treat[,seq(3, ncol(treat), by = 4)]
treat <- na.omit(treat)
treat <- treat[!is.infinite(rowSums(treat)),]


## Correlation 
conCor <- cor(con, method = "pearson")
colnames(conCor) <- paste0("Rep", 1:ncol(conCor))
conVals <- apply(conCor, 2, list) |> unlist() |> unique()
conMax <- max(conVals[-1])
conMin <- min(conVals[-1])
conCor[conCor==1] <- conMax

treatCor <- cor(treat, method = "pearson")
colnames(treatCor) <- paste0("Rep", 1:ncol(treatCor))
treatVals <- apply(treatCor, 2, list) |> unlist() |> unique()
treatMax <- max(treatVals[-1])
treatMin <- min(treatVals[-1])
treatCor[treatCor==1] <- treatMax

## Try getting upper/lower triangles first 
## Plot without the help of the corplot 
conCor <- cor(con, method = "pearson")
upper <- conCor[upper.tri(conCor)]
conCor[upper.tri(conCor)] <- 0
conCor <- as.data.frame(conCor)
conCor$name <- rownames(conCor)
conCor <- pivot_longer(conCor, -name, names_to = "rep") |> 
  as.data.frame()
conCor$x <- c(rep(1, 9), rep(2, 9), rep(3, 9), 
              rep(4, 9), rep(5, 9), rep(6, 9), 
              rep(7, 9), rep(8, 9), rep(9, 9))
conCor$y <- rep(1:9, 9)
conCor <- na.omit(conCor)
conCor <- conCor[!conCor$value==1,]
conCor$color <- c(conCor[!conCor$value==0,]$value, conCor[!conCor$value==0,]$value)

label <- data.frame(rep = paste0("Rep", 1:9), 
                    x = 1:9, 
                    y = 1:9, 
                    value = 1)

## figure s1d
conHeatmap <- 
  ggplot(conCor, aes(x = x, y = y, fill = value)) + 
  geom_tile(color = "gray") + 
  theme_classic() + 
  xlab("") + ylab("") + 
  xlim(0,10) + ylim(0,10) + 
  ggtitle("Resting") + 
  scale_fill_gradientn(colors = c("white", brewer.pal(8, "Blues")), limits = c(0, 1), guide = guide_colorbar(title = "Pearson's R", title.position = "right")) + 
  geom_text(data = conCor[conCor$value==0,], aes(color = color), label = round((conCor[!conCor$value==0,]$value), 2), size = 2) + 
  scale_color_gradientn(colors = c("white", brewer.pal(8, "Blues")), limits = c(0, 1), guide = "none") + 
  geom_text(data = label, aes(x = x, y = y), label = label$rep, color = "black", size = 2) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        line = element_blank(), 
        legend.key.height = unit(1, "cm"), 
        legend.margin = margin(l = -0.1, unit='cm'), 
        legend.title = element_text(angle = -90), 
        legend.title.align = 0.5)

#Treated
treatCor <- cor(treat, method = "pearson")
upper <- treatCor[upper.tri(treatCor)]
treatCor[upper.tri(treatCor)] <- 0
treatCor <- as.data.frame(treatCor)
treatCor$name <- rownames(treatCor)
treatCor <- pivot_longer(treatCor, -name, names_to = "rep") |> 
  as.data.frame()
treatCor$x <- c(rep(1, 6), rep(2, 6), rep(3, 6), 
                rep(4, 6), rep(5, 6), rep(6, 6))
treatCor$y <- rep(1:6, 6)
treatCor <- na.omit(treatCor)
treatCor <- treatCor[!treatCor$value==1,]
treatCor$color <- c(treatCor[!treatCor$value==0,]$value, treatCor[!treatCor$value==0,]$value)

label <- data.frame(rep = paste0("Rep", 1:6), 
                    x = 1:6, 
                    y = 1:6, 
                    value = 1)

## figure s1e
treatHeatmap <- 
  ggplot(treatCor, aes(x = x, y = y, fill = value)) + 
  geom_tile(color = "gray") + 
  theme_classic() + 
  xlab("") + ylab("") + 
  xlim(0,10) + ylim(0,10) + 
  ggtitle("Activated") + 
  scale_fill_gradientn(colors = c("white", brewer.pal(8, "Blues")), limits = c(0, 1), guide = guide_colorbar(title = "Pearson's R", title.position = "right")) + 
  geom_text(data = treatCor[treatCor$value==0,], aes(color = color), label = round((treatCor[!treatCor$value==0,]$value), 2), size = 2) + 
  scale_color_gradientn(colors = c("white", brewer.pal(8, "Blues")), limits = c(0, 1), guide = "none") + 
  geom_text(data = label, aes(x = x, y = y), label = label$rep, color = "black", size = 2) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        line = element_blank(), 
        legend.key.height = unit(1, "cm"), 
        legend.margin = margin(l = -2, unit='cm'), 
        legend.title = element_text(angle = -90), 
        legend.title.align = 0.5)

