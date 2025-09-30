## Script to read in the output from epiRMextractParallel to get DNase/ATAC for all cell types and visualize their activity per quantile
## Generates figure 3b

## Load libraries
library(preprocessCore)
library(data.table)

## Read in data
files <- list.files("/work/users/m/a/marielle/work/AD3D/data/epiRM/", "txt", full.names = TRUE)
macMic <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/mac_mgl_dipScore.txt", sep = "\t", header=T)
data <- lapply(files, fread, header=T)

## ATAC VERSION: 
files <- list.files("/work/users/m/a/marielle/work/AD3D/data/epiRM/atac/", "txt", full.names = TRUE)
data <- lapply(files, fread, header=T)
mac <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/mac_mgl_atac.txt", sep = "\t", header=T)
mac <- mac[mac$celltype=="THP1",] |> as.data.table()

## Read in the microglia data from Yin's paper here 
mic <- read.table("/work/users/m/a/marielle/work/AD3D/data/epiRM/atac/IMGL_ATAC_data.txt", sep = "\t", header = TRUE)
## put mac and mic into list 
macMic <- list(mac, mic)

## add to the list 
data <- c(data, macMic)


## First, for each cell type we need to get an average value for each of the 10 quantiles (collapse the boxplots by their mean)
means <- list()
for (i in 1:length(data)){
  d <- data[[i]]
  d <- 
    d |> 
    group_by(quantile) |> 
    summarise(mean = mean(atac), name = unique(d$celltype)) |> 
    as.data.frame()
  means[[i]] <- d
}
means <- do.call(rbind, means) |> 
  as.data.frame()

## First quantile normalize across cell types 
bycell <- pivot_wider(means, id_cols = name, names_from = quantile, values_from = mean, values_fn = mean)
bycell <- as.matrix(bycell[,2:11])
rownames(bycell) <- unique(means$name)
## Normalizing for read depth due to the fact that there's such large differences in the THP1/MGL datasets: 
scaleFactors <- apply(bycell, 1, sum)
## For each row, divide by its respective scaling factor
for (i in 1:nrow(bycell)){
  bycell[i,] <- bycell[i,]/scaleFactors[i]
}

## Now that we've read depth normalized and quantile normalized, can account for differences within each cell type across the quantiles: 
bycell <- as.data.frame(bycell)
colnames(bycell) <- 1:10
bycell$name <- rownames(bycell)

## Pivot the data 
data <- pivot_longer(bycell, !name, names_to = "quantile", values_to = "qnorm")
## normalize by quantile
data <- 
  data |> 
  group_by(quantile) |> 
  mutate(zscore = (qnorm - mean(qnorm))/sd(qnorm))
## normalize by celltype
data <- 
  data |> 
  group_by(name) |> 
  mutate(zscore = (qnorm - mean(qnorm))/sd(qnorm))

## Determine order of data
## Figure out another order where for each cell type we figure out where the max quantile is, order them like that: 
data <- 
  data |> 
  group_by(name) |> 
  mutate(cluster = which.max(zscore))
## Within each cluster, order from max to min: 
cluster <- 
  data |> 
  group_by(cluster) |> 
  mutate(cluster = which.max(zscore))

clusterOrder <- 
  data |> 
  arrange(desc(cluster), desc(zscore)) |> 
  as.data.frame()
data$name <- factor(data$name, levels = unique(clusterOrder$name))

## Meta data to get actual names for cell types: 
meta <- read.csv("/proj/phanstiel_lab/External/consortium/epiRM/EGRM_meta.txt", sep = "\t") |> 
  as.data.frame()
meta <- meta[meta$Epigenome.ID..EID. %in% means$name,]
colnames(meta)[2] <- "name"
## Replace the EpiRM names with actual names
merged <- left_join(data, meta, by = "name")
## Factor quantiles: 
merged$quantile <- as.character(merged$quantile)
merged$quantile <- factor(merged$quantile, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
merged$name <- factor(merged$name, levels = unique(clusterOrder$name))

## Ordering when the data is just read depth normalized
order <- 
  data |> 
  filter(quantile==10) |> 
  group_by(name) |> 
  summarise(max = max(zscore)) |> 
  as.data.frame()
order <- order[order(order$max, decreasing = FALSE),]$name
data$name <- factor(data$name, levels = order)

## Figure out another order where for each cell type we figure out where the max quantile is, order them like that: 
data <- pivot_longer(bycell, !name, names_to = "quantile", values_to = "norm")
data <- 
  data |> 
  group_by(name) |> 
  mutate(cluster = which.max(norm))
## Within each cluster, order from max to min: 
cluster <- 
  data |> 
  group_by(cluster) |> 
  mutate(cluster = which.max(norm))

clusterOrder <- 
  data |> 
  arrange(desc(cluster), desc(norm)) |> 
  as.data.frame()
data$name <- factor(data$name, levels = unique(clusterOrder$name))

## Add colors to macrophages and microglia 
merged[merged$name %in% c("THP1"),]$COLOR <- "#227B7F"
merged[merged$name %in% c("THP1"),]$GROUP <- "Macrophage"
merged[merged$name %in% c("IMGL"),]$COLOR <- "pink"
merged[merged$name %in% c("IMGL"),]$GROUP <- "Microglia"

## Play with sizes
merged$size <- ""
merged[merged$quantile==10 & merged$GROUP=="Macrophage",]$size <- "Large"
merged[merged$quantile==10 & merged$GROUP=="Microglia",]$size <- "Large"
merged[merged$quantile==10 & merged$name=="E081",]$size <- "Large"
merged[merged$quantile==10 & merged$name=="E082",]$size <- "Large"
merged[merged$quantile==10 & merged$name=="E082",]$size <- "Large"
merged[merged$quantile==10 & merged$name=="E124",]$size <- "Large"

## Make the ENCODE2012 samples less intense 
merged[merged$GROUP=="ENCODE2012",]$COLOR <- "#AFAFAF"
pal <- unique(merged$COLOR)
names(pal) <- unique(merged$GROUP)
#Try sending encode to back: 
order <- unique(merged$GROUP)
order <- order[c(14,1:13,15:16)]
merged$GROUP <- factor(merged$GROUP, levels = order)

## figure 3b
epiRM_dotPlot <- 
  ggplot(merged, aes(x = quantile, y = qnorm, color = GROUP, size = size)) + 
  geom_point(data = merged[merged$GROUP=="ENCODE2012",]) +
  geom_point(data = merged[!merged$GROUP=="ENCODE2012",]) +
  scale_color_manual(values = pal) + 
  theme_minimal() + 
  scale_size_manual(values = c(0.75, 2)) + 
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.key.size = unit(0.2, "cm"), 
        legend.text = element_text(size = 5), 
        legend.title = element_blank(), 
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  xlab("MPRA Activity Quantile") + ylab("Normalized Accessibility") + 
  guides(size = "none", 
         color = guide_legend(nrow = 4))

