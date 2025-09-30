## Script to compare our MPRA data to the HEK MRPA published in Cooper 2022 

## Load libraries
library(readxl)

## Read in data
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T, sep = "\t")
cooper <- read_xlsx("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s5.xlsx", sheet = 8, col_names = TRUE, skip = 3) |>
  as.data.frame()
colnames(cooper)[2] <- "rsid"
cooper$sig <- ""
cooper[cooper$q.value<0.01,]$sig <- TRUE
cooper[cooper$sig=="",]$sig <- FALSE
cooper$sign <- sign(cooper$Log2.FC)
## Apparently they calculated their log2.fc by doing log(alt minus ref), not divided by???
## Try reading in the raw activity score so we can calculate mpra activity the same way that we do in our MPRA
cooperRaw <- read_xlsx("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s1.xlsx", sheet = 2, col_names = TRUE, skip = 3) |>
  as.data.frame()
colnames(cooperRaw)[2] <- "rsid"
## Order by the same order as the other df 
cooperRaw <- cooperRaw[order(match(cooperRaw$rsid, cooper$rsid)),]
## which are AD variants?
cooperAD <- cooperRaw[cooperRaw$Disorder=="AD",]$rsid
cooper <- cooper[cooper$rsid %in% cooperAD,]


## For supplemental figure--first make a venn diagram showing which variants were even tested in both of our datasets: 
adjusted$tested <- TRUE
cooper$tested <- TRUE
merged <- full_join(adjusted, cooper, by = "rsid")
merged <- merged[,c(1,19,34)]
merged[is.na(merged$tested.x),]$tested.x <- FALSE
merged[is.na(merged$tested.y),]$tested.y <- FALSE
colnames(merged) <- c("rsid", "THP1", "HEK")

## figure s5a
ggplot(merged, aes(A = THP1, B = HEK)) + 
  geom_venn(auto_scale = TRUE, fill_color = c("seagreen3", "darkturquoise"), fill_alpha = 0.8, stroke_size = 0, text_size = 2, set_name_size = 3, show_percentage=FALSE) + 
  theme_void() + 
  coord_fixed() + 
  ggtitle("Variants tested in MPRA") + 
  theme(plot.title = element_text(size = 10, hjust = 0.4, vjust = -8))

## Of the 1155 variants tested in both, how many are MPRA sig in each? 
merged <- merge(adjusted, cooper, by = "rsid")
## We need to re-calculate padj for just these tested variants in each cell type: 
merged$thp1_padj <- p.adjust(merged$P.Value, method = "BH")
merged$hek_padj <- p.adjust(merged$p.value, method = "BH")
#subset certain columns only: 
subset <- merged[,grep("padj", colnames(merged))]
subset$rsid <- merged$rsid
subset$thp1_sig <- FALSE
subset$hek_sig <- FALSE
subset[subset$thp1_padj<0.05,]$thp1_sig <- TRUE
subset[subset$hek_padj<0.05,]$hek_sig <- TRUE
subset <- subset[,3:5]
colnames(subset) <- c("rsid", "THP1", "HEK")
## Make barplots showing the number of mpra allelic vars in each celltype: 
bars <- apply(subset[2:3], 2, table) |> 
  as.data.frame()
bars <- pivot_longer(bars, cols = colnames(bars), names_to="celltype") |> 
  as.data.frame()
bars$group <- c("Neg", "Neg", "THP1_sig", "HEK_sig")
bars$group <- factor(bars$group, levels = c("Neg", "THP1_sig", "HEK_sig"))
bars$celltype <- factor(bars$celltype, levels = c("THP1", "HEK"))

## figure s5b
ggplot(bars, aes(x = celltype, y = value, fill = group)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.5) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        aspect.ratio = 2) + 
  xlab("") + ylab("Number of MPRA allelic\n variants (padj<0.05)") + 
  scale_fill_manual(values = c("gray", "seagreen3", "darkturquoise")) + 
  geom_text(label = bars$value, vjust = 1.5, size = 3)

## figure s5c
ggplot(subset[subset$THP1==TRUE | subset$HEK==TRUE,], aes(A = THP1, B = HEK)) + 
  geom_venn(auto_scale = TRUE, fill_color = c("seagreen3", "darkturquoise"), fill_alpha = 0.8, stroke_size = 0, text_size = 2, set_name_size = 3, show_percentage=FALSE) + 
  theme_void() + 
  coord_fixed() + 
  ggtitle("Significant Allelic Variants\n (padj<0.05)") + 
  theme(plot.title = element_text(hjust = 0.45, vjust = -10, size = 10))

## For those 7 variants significant in both sets: 
shared <- subset[subset$THP1==TRUE & subset$HEK==TRUE,]$rsid
directions <- merged[merged$rsid %in% shared,]
directions <- directions[grep("sign", colnames(directions))]
ide <- apply(directions, 1, prod) |> 
  table() |> 
  as.data.frame()
ide$x <- "x"
ide$ide <- c(FALSE, TRUE)

## figure s5d
ggplot(ide, aes(x = Freq, y = x, fill = ide)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.5) + 
  xlab("Identical Direction of Effect\nbetween THP1/HEK") + ylab("") + 
  theme_minimal() + 
  theme(legend.position = "bottom", 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.key.size = unit(0.5, "cm")) + 
  labs(fill = "") + 
  geom_text(label = rev(ide$Freq), hjust = -3, size = 3) + 
  scale_fill_manual(values = c("gray70", "lightblue"))


## Our data: 
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch_N9N6_combined_batchvar.txt", header=T)
batchvar <- na.omit(batchvar)
## Add coordinates: 
batchvar$chr <- unlist(strsplit(batchvar$SNP, "_"))[seq(1, (2*nrow(batchvar)), by = 2)]
batchvar$pos <- unlist(strsplit(batchvar$SNP, "_"))[seq(2, (2*nrow(batchvar)), by = 2)]
batchvar$allele <- rep(c("a1", "a2"),(nrow(batchvar))/2)
## Separate logratios in each condition: 
logratio <- batchvar[grep("log", colnames(batchvar))]
## Get median across replicates
logratio$median <- apply(logratio, 1, median)
## Combine into dataframe: 
thp <- data.frame(id = batchvar$variant,
                  logratio = logratio$median,
                  chr = batchvar$chr, 
                  pos = batchvar$pos, 
                  allele = batchvar$allele)

## reformat the id column
thp$variant <- unlist(strsplit(thp$id, "_ad"))[seq(1,2*nrow(thp), by=2)]
cooperRaw$variant <- paste0(cooperRaw$chr, "_", cooperRaw$pos)

## Subset for just the variants that are in both assays: 
assayed <- merged[merged$tested.x==TRUE & merged$tested.y==TRUE,]$id
hek <- cooperRaw[cooperRaw$variant %in% assayed,]
thp <- thp[thp$variant %in% assayed,]

## Scale log ratio around 0 
hek$activity <- scale(hek$Log2.FC, center = TRUE, scale = FALSE)
thp$logratio <- scale(thp$logratio, center = TRUE, scale = FALSE)

## Assign variants to quantiles
hek <- hek[order(hek$activity, decreasing = TRUE),]
thp <- thp[order(thp$logratio, decreasing = TRUE),]
n <- 10

hek <-
  hek %>%
  mutate(quantile = ntile(hek$activity, n))
hek$quantile <- factor(hek$quantile, levels = c(1:n))
thp <- 
  thp %>%
  mutate(quantile = ntile(thp$logratio, n))
thp$quantile <- factor(thp$quantile, levels = c(1:n))

## Determine how many of the hek top are thp top: 
topHek <- hek[hek$quantile==10,]$variant
topThp <- thp[thp$quantile==10,]$variant

## Read in signal data: 
macBW <- "/work/users/m/a/marielle/work/AD3D/data/lima/LIMA_atac_mergedSignal.bigwig"

## merge the hek and thp into one df so we only have to extract counts once per celltype
mergedData <- merge(thp, hek, by = "variant")
colnames(mergedData)[c(7,23)] <- c("THP1_quantile", "HEK_quantile")

## Just look at quant10 variants 
mergedData <- mergedData[mergedData$HEK_quantile==10 | mergedData$THP1_quantile==10,]

buffer <- 200
mergedData$THP1_atac <- ""
i <- 1
buffer <- 200
for (i in 1:nrow(mergedData)){
  print(i)
  pos <- as.numeric(mergedData$pos.x[i])
  #THP
  atacTHP <- readBigwig(file = macBW,
                        chrom = mergedData$chr.x[i],
                        chromstart = pos-buffer,
                        chromend = pos+buffer)
  atacScore <- max(atacTHP$score)
  if (length(atacScore)==0){
    atacScore <- 0
  }
  mergedData$THP1_atac[i] <- atacScore
}
## Make sure that ATAC is numeric
mergedData$THP1_atac <- as.numeric(mergedData$THP1_atac)
mergedData <- mergedData[!mergedData$THP1_atac=="-Inf",]
mergedData <- pivot_longer(mergedData, cols = c("THP1_atac"), names_to = "celltype", values_to = "atac") |> as.data.frame()

## Plotting the quant 10 for each set of variants: 
topthp <- mergedData[mergedData$THP1_quantile==10,]
topHek <- mergedData[mergedData$HEK_quantile==10,]
topthp$group <- "THP1"
topHek$group <- "HEK"
##combine 
data <- rbind(topthp[,c(1,25,26)], 
              topHek[,c(1,25,26)])

## figure s5e
ggplot(data, aes(x = group, y = log(atac), fill = group)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  scale_fill_manual(values = c("darkturquoise", "seagreen3")) + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  xlab("") + ylab("log(THP1 ATAC)") + 
  geom_hline(yintercept = 2.56, lty = 3)

