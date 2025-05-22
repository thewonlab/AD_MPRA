### Show activity quantile and dip scores for conditions separately
### Generates figures 3a, s4a-d

## Read in data
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch_N9N6_combined_batchvar.txt", header=T)
batchvar <- na.omit(batchvar)
## Add coordinates: 
batchvar$chr <- unlist(strsplit(batchvar$SNP, "_"))[seq(1, (2*nrow(batchvar)), by = 2)]
batchvar$pos <- unlist(strsplit(batchvar$SNP, "_"))[seq(2, (2*nrow(batchvar)), by = 2)]
batchvar$allele <- rep(c("a1", "a2"),(nrow(batchvar))/2)
## Separate logratios in each condition: 
logratio <- batchvar[grep("log", colnames(batchvar))]
con <- logratio[grep("con", colnames(logratio))]
treat <- logratio[grep("treat", colnames(logratio))]
## Get median across replicates
con$median <- apply(con, 1, median)
treat$median <- apply(treat, 1, median)

## Combine into dataframe: 
con <- data.frame(id = batchvar$variant,
                  logratio = con$median,
                  chr = batchvar$chr, 
                  pos = batchvar$pos, 
                  allele = batchvar$allele,
                  class = "control")
treat <- data.frame(id = batchvar$variant,
                    logratio = treat$median,
                    chr = batchvar$chr, 
                    pos = batchvar$pos, 
                    allele = batchvar$allele,
                    class = "treat")

## Scale log ratio around 0 
con$logratio <- scale(con$logratio, center = TRUE, scale = FALSE)
treat$logratio <- scale(treat$logratio, center = TRUE, scale = FALSE)

## Add  quantile information
con <- con[order(con$logratio, decreasing = TRUE),]
treat <- treat[order(treat$logratio, decreasing = TRUE),]
n <- 10

con <- 
  con %>%
  mutate(quantile = ntile(con$logratio, n))
con$quantile <- factor(con$quantile, levels = c(1:n))
treat <- 
  treat %>%
  mutate(quantile = ntile(treat$logratio, n))
treat$quantile <- factor(treat$quantile, levels = c(1:n))

## Combine to plot in one plot: 
data <- rbind(con, treat)
data$quantile <- as.numeric(data$quantile)
data$quantile <- factor(data$quantile, levels = c(1:n))

## Plot and facet wrap by condition: 
pal <- colorRampPalette(c("#CEE6A1", "#9FD141"))
## figure s4a
mpraQuantilesbyCondition <- 
  ggplot(data, aes(x = quantile, y = logratio, fill = quantile)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() + 
  theme(legend.position = "none", panel.spacing = unit(0.5, "lines"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  facet_wrap(~class) + 
  xlab("MPRA Activity Quantile") + 
  ylab("log(RNA/DNA)") + 
  scale_fill_manual(values = pal(10)) + 
  geom_hline(yintercept = 0, lty = 3) + 
  coord_cartesian(ylim = c(-0.5, 1.1))

## Enhancer activity: 
enhBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/chip/H3K27ac/LIMA_h3k27ac_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]

## ATAC activity: 
atacBW <- list.files("/proj/phanstiel_lab/Data/processed/LIMA/atac/LIMA_ATAC_THP1_WT_LPIF_S/signal/", "MERGE", full.names = TRUE)[c(1,8)]

## For each enhancer sequence, read in the bigwig for the enhancers: 
i <- 1
buffer <- 200

mpra <- list(con, treat)

for (x in 1:length(mpra)){
  mpraCounts <- mpra[[x]]
  mpraCounts$dipScore <- ""
  mpraCounts$atac <- ""
  mpraCounts$pos <- as.numeric(mpraCounts$pos)
  
  for (i in 1:nrow(mpraCounts)){
    print(i)
    data <- readBigwig(file = enhBW[[x]],
                       chrom = mpraCounts$chr[i],
                       chromstart = mpraCounts$pos[i]-buffer,
                       chromend = mpraCounts$pos[i]+buffer)
    scores <- data$score
    upstream <- scores[1]
    dnstream <- scores[length(scores)]
    
    dip <- scores[round(length(scores)/2)]
    dipScore <- (upstream+dnstream)-(2*dip)
    if (length(dipScore)==0){
      dipScore <- 0
    }
    if (dipScore < 0){
      dipScore <- 0
    }
    mpraCounts$dipScore[i] <- as.numeric(dipScore)
    
    ## Also add on ATAC score: 
    atac <- readBigwig(file = atacBW[x],
                       chrom = mpraCounts$chr[i],
                       chromstart = mpraCounts$pos[i]-buffer,
                       chromend = mpraCounts$pos[i]+buffer)
    atacScore <- max(atac$score)
    if (length(atacScore)==0){
      atacScore <- 0
    }
    mpraCounts$atac[i] <- as.numeric(atacScore)
  }
  mpra[[x]] <- mpraCounts
}

## Recombine both components of list into one dataframe: 
data <- do.call(rbind, mpra)
data$dipScore <- as.numeric(data$dipScore)
data$atac <- as.numeric(data$atac)

## Save this so output to save time later
write.table(data, file = "/work/users/m/a/marielle/work/AD3D/data/mpra/dipScorePerCondition.txt", sep = "\t", quote = F)

## And re-read it for running later 
data <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/dipScorePerCondition.txt", header=T)
data$quantile <- as.character(data$quantile)
data$quantile <- factor(data$quantile, levels = as.character(1:10))
data[data$class=="control",]$class <- "resting"
data[data$class=="treat",]$class <- "active"
data$class <- factor(data$class, levels = c("resting", "active"))

## Plot: 
## figure s4b
dipScorebyCondition <- 
  ggplot(data, aes(x = quantile, y = log(dipScore), fill = class)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() + 
  theme(legend.position = "none", panel.spacing = unit(0.5, "lines"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  facet_wrap(~class) + 
  xlab("MPRA Activity Quantile") + 
  ylab("log(Dip Score)") + scale_fill_manual(values = c("#F0B0AA", "#EC6F61"))

## figure 3a
atacByCondition <- 
  ggplot(data, aes(x = quantile, y = log(atac), fill = class)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() + 
  theme(legend.position = "none", panel.spacing = unit(0.5, "lines"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  facet_wrap(~class) + 
  xlab("MPRA Activity Quantile") + 
  ylab("log(ATAC)") + 
  scale_fill_manual(values = c("#F0B0AA", "#EC6F61"))

## Look at dip scores and atac over the inactive mpra variants vs the active mpra variants (and context specific mpra variants) 
## add unique identifier to data
data[data$allele=="a1",]$allele <- "alt"
data[data$allele=="a2",]$allele <- "ref"
#data$id <- paste0(data$chr, "_", data$pos, "_", data$allele)
data$id <- paste0(data$chr, "_", data$pos)
## Classify what group the allele fits into: 
### LMER VERSION: 
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

data$activity <- "inactive"
data[data$id %in% active,]$activity <- "active"

## figure s4c
activeVsInactive <- 
  ggplot(data, aes(x = activity, y = log(dipScore), fill = activity)) + 
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~class) + 
  theme_minimal() + 
  theme(legend.position = "none", panel.spacing = unit(0.5, "lines"), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  scale_fill_manual(values = c("#9FD141", "#CEE6A1")) + 
  xlab("") + ylab("log(Dip Score)")

## figure s4d
activeVsInactiveAtac <- 
  ggplot(data, aes(x = activity, y = log(atac), fill = activity)) + 
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~class) + 
  theme_minimal() + 
  theme(legend.position = "none", panel.spacing = unit(0.5, "lines"), 
        axis.title.y = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6)) + 
  scale_fill_manual(values = c("#9FD141", "#CEE6A1")) + 
  xlab("") + ylab("log(ATAC)")
