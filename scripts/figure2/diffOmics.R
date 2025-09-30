## Script to summarize all omic data with pie charts: 

## Load libraries
library(ggplot2)

## Read in data: Note--all of this data is from the supplementary material from Reed et al 2022 Cell Reports 
hic <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/hic/hicLFC.txt", header=TRUE)
atac <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/atac/LIMA_atacLFC.txt", header=TRUE)
enh <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/cnr/LIMA_enhLFC.txt", header=TRUE)
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt", header=TRUE)

## Combine into list
data <- list(hic, atac, enh, rna)
names(data) <- c("Loops", "ATAC", "H3K27ac", "RNA")
labels <- c("loops", "peaks", "peaks", "genes")

## For each dataset, make a pie chart showing differential vs non differential: 
pies <- list()
for (i in 1:length(data)){
  df <- data[[i]]
  df[!df$class=="static",]$class <- "differential"
  df <- 
    table(df$class) |>
    as.data.frame()
  percent <- (df[df$Var1=="differential",]$Freq)/(sum(df$Freq))
  percent <- round(percent, digits = 4)*100
  total <- sum(df$Freq)
  plot <- 
    ggplot(df, aes(x = "", y = Freq, fill = Var1)) + 
    geom_bar(stat = "identity") + 
    coord_polar("y", start = 0) + 
    theme_void() + 
    scale_fill_manual(values = c("#56A2CC", "#D7DDE6")) + 
    ggtitle(label = names(data)[i]) + 
    labs(caption = paste0(total, " total ", labels[i], "\n", df[df$Var1=="differential",]$Freq, " DE ", labels[i], "\n(", percent, "%)")) + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 8, hjust = 0.5, vjust = -1), 
          plot.caption = element_text(size = 3.5, hjust = 0.5, vjust = 6))
  pies[[i]] <- plot
}
