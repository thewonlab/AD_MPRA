## Compare the results of our RNA-seq to the scRNA-seq clusters from Sun et al 2023
## Generates figure 2d

## Load libraries
library(readxl)

## Compare to Sun 2023 scRNA-seq microglia data: 
markers <- read_xlsx("/work/users/m/a/marielle/external/sun2023/1-s2.0-S0092867423009716-mmc1.xlsx", sheet = 3) |> 
  as.data.frame()
states <- unique(markers$microgliaState)

## Read in data: 
rna <- read.table("/work/users/m/a/marielle/work/AD3D/LIMA/rna/LIMA_rnaLFC.txt")
rna$ENSEMBL <- rownames(rna)

## If RNA has transcripts: 
rna$ENSEMBL <- unlist(strsplit(rna$ENSEMBL, "[.]"))[seq(1,(2*nrow(rna)),by=2)]

## Get RNA coordinate information: 
genes <- read.table("/work/users/m/a/marielle/ref/hg19_annotation_txdb.txt", header = T)
colnames(genes)[7] <- "gene"
rna <- merge(genes, rna, by = "ENSEMBL")

## Combine the markers with RNA information: 
data <- merge(rna, markers, by = "gene")

## Reformat
data <- data[,c(1,3,4,5,8,9,10,11,12,18)]
data$microgliaState <- factor(data$microgliaState, levels = ((c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG10", "MG11", "MG12"))))
pal <- (c("#EE2329", "#0094C8", "#24B355", "#BA65A8", "#F68A1F", "#FFEA16", "#995D2D", "#F391BC", "#AAAAA9", "#3FBDAC", "#F7986F", "#90B0D5"))

annoCluster <- data.frame(microgliaState = states, 
                          name = rev(c("Cycling", "Antiviral", "Inflammatory III", "Inflammatory II", "Glycolytic", "Stress Signature", "Phagocytic", "Lipid Processing", "Ribosome Biogenesis", "Inflammatory I", "Neuronal Surveillance", "Homeostatic")))
annoCluster$pvalue <- ""
#annoCluster$sig <- rev(c("*", "*", "*", "*", "*", "*", "*", "*", "ns", "*", "ns", "ns"))


## Determine which are statistically different than zero: 
for (i in 1:12){
  cluster <- annoCluster$microgliaState[[i]]
  tmp <- data[data$microgliaState==cluster,]
  annoCluster$pvalue[[i]] <- as.numeric(wilcox.test(tmp$log2FoldChange)$p.value)
}
annoCluster$pvalue <- as.numeric(annoCluster$pvalue)
annoCluster$pvalue <- round(annoCluster$pvalue, digits = 8)
## adjust pvalue
annoCluster$FDR <- p.adjust(annoCluster$pvalue, method = "BH")
## Simplify: 
annoCluster$sig <- "ns"
annoCluster[annoCluster$pvalue<0.05,]$sig <- "*"

data <- merge(data, annoCluster, by = "microgliaState")
data$name <- factor(data$name, levels = (c(annoCluster$name)))
data$microgliaState <- factor(data$microgliaState, levels = rev(c("MG0", "MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG10", "MG11", "MG12")))


## Visualize: 
## Figure 2d
mgStatesRNA <- 
  ggplot(data, aes(y = microgliaState, x = log2FoldChange, fill = name)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() + 
  geom_hline(yintercept = 0, lty = 3) + 
  ylab("log2(Activated/Resting)") + xlab("") + 
  theme(axis.title.x = element_text(size = 6), 
        axis.text.x = element_text(size = 5),         
        axis.title.y = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 5), 
        legend.key.size = unit(0.1, "cm"), 
        legend.position = "right") + 
  scale_fill_manual(values = pal) +
  annotate("text", label = rev(annoCluster$sig), y = 1:12, x = 7, size = 3) + 
  geom_vline(xintercept = 0, lty = 3)


