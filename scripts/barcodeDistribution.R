## Script to determine how many barcodes we recovered from each variant: 
## Generates figure s1c

## Barcode map 
load("/work/users/m/a/marielle/proj/AD/MPRA/bcmap/fastq/unique_bcdat.rda") ## available from GEO (GSE273887_AD_MPRA_Barcode_Map.txt.gz)
head(bcdat)
bcdat <- table(bcdat$variant)
bcdat <- as.data.frame(bcdat)

## Data
load("/work/users/m/a/marielle/work/batch3_NYGC/mpraOutput/mergedCounts_noOutliers_atLeastFiveBCs_perfectmatch.rda") ## available from GEO (GSE273887_AD_MPRA_countTable.txt.gz)
batch <- as.data.frame(batch)
batch <- table(batch$variant)
batch <- as.data.frame(batch)

## Plot: 
## This is the distribution of all barcodes in the construct, not currently in any figure
allBcHist <- 
  ggplot(bcdat, aes(x = Freq)) + 
  geom_histogram(color = "black", lwd = 0.25, fill = "gray") + 
  theme_minimal() + 
  xlab("Number of Unique Barcodes \nPer Variant (in BC map)") + 
  ylab("Frequency") + 
  geom_vline(xintercept = median(bcdat$Freq)) + 
  annotate("text", x = 1200, y = 800, label = paste0("Median = ", median(bcdat$Freq)), color = "firebrick", size = 3) + 
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8))

## This is the distribution of all barcodes detected in the data (figure s1c)
batchBCHist <- 
  ggplot(batch, aes(x = Freq)) + 
  geom_histogram(color = "black", lwd = 0.25, fill = "gray") + 
  theme_minimal() + 
  xlab("Number of Unique Barcodes \nPer Variant (in batch3)") + 
  ylab("Frequency") + 
  geom_vline(xintercept = median(batch$Freq)) + 
  annotate("text", x = 250, y = 700, label = paste0("Median = ", median(batch$Freq)), color = "firebrick", size = 3) + 
  theme(axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8))


