## Script to compare the epigenetic profiles of THP-1, Primary Macrophages, induced-pluripotent stem cell derived Microglia, and HEK293T cels
## Generates the figures S5f-k

## Load libraries: 
library(ggvenn)
library(data.table)
library(ggplot2)
library(GenomicRanges)

## Read in data: 
## Our THP-1 macrophages (resting only):
mac <- read.table("/proj/phanstiel_lab/Data/processed/LIMA/atac/LIMA_ATAC_THP1_WT_LPIF_S/peaks/MERGE_LIMA_ATAC_ATAC_THP1_WT_LPIF_0000_S_peaks.narrowPeak")
mac <- mac[,1:3]
colnames(mac) <- c("chr", "start", "end")
mac <- GRanges(mac)
mac$peak <- paste0("macPeak", 1:length(mac))

## iMGLs from Yang 2023 (resting only)
mic <- fread("/work/users/m/a/marielle/external/yang2023/GSE206479_WTC11_rest_idr_optimal_peak.narrowPeak.gz")
mic <- mic[,1:3]
colnames(mic) <- c("chr", "start", "end")
mic <- GRanges(mic)
mic$peak <- paste0("MicPeak", 1:length(mic))

## Alasoo primary macrophages
priMac <- fread("/work/users/m/a/marielle/external/alasoo2018/ATAC_peak_metadata.txt.gz") |>
  as.data.frame()
priMac <- priMac[,c(4,5,6,7,1)]
priMac <- GRanges(priMac)
seqlevelsStyle(priMac) <- "UCSC"
## liftover to hg19
chain <- import.chain("/work/users/m/a/marielle/ref/liftOver/hg38ToHg19.over.chain")
priMac <- unlist(liftOver(priMac, chain))

## hek atac
hek <- fread("/work/users/m/a/marielle/external/dong2024/GSE266490_Combined_ATAC_rep9-16_peak_counts.txt.gz")
counts <- hek[,c(4:13)] |> 
  as.data.frame()
counts$median <- apply(counts[,3:10], 1, median)
## order them by the qvalue
counts <- counts[order(counts$qvalue, decreasing=TRUE),]
#pick the ~140k strongest since there's lots more HEK peaks that microglia/macrophage
peaks <- counts[1:140000,]$name
hek <- hek[hek$name %in% peaks,]
hek <- hek[,1:3]
colnames(hek) <- c("chr", "start", "end")
hek <- GRanges(hek)
hek$peak <- paste0("HEKPeak", 1:length(hek))
## liftover to hg19
chain <- import.chain("/work/users/m/a/marielle/ref/liftOver/hg38ToHg19.over.chain")
hek <- unlist(liftOver(hek, chain))

## data frame for which comparisons we want to look at: 
df <- data.frame(celltype1 = c("THP1", "THP1", "HEK", "HEK"), 
           celltype2 = c("priMac", "iMGL", "priMac", "iMGL"),
           color1 = c("#ADD8E6", "#ADD8E6", "#2E8B57", "#2E8B57"), 
           color2 = c("#4682B4", "#FFC0CB", "#4682B4", "#FFC0CB"))
df$cell2overlap <- ""

for (i in 1:nrow(df)){
  print(i)
  row <- df[i,]
  if (row$celltype1=="THP1"){
    gr1 <- mac
  } else {
    gr1 <- hek
  }
  
  if (row$celltype2=="priMac"){
    gr2 <- priMac
  } else {
    gr2 <- mic
  }
  
  gr_union <- GenomicRanges::reduce(c(gr1, gr2))
  
  # Check for overlaps with each set
  overlap1 <- countOverlaps(gr_union, gr1) > 0  # TRUE if the range is in gr1
  overlap2 <- countOverlaps(gr_union, gr2) > 0  # TRUE if the range is in gr2
  
  # Create a data frame from the GRanges object and the overlap info
  result_df <- data.frame(
    seqnames = seqnames(gr_union),
    start = start(gr_union),
    end = end(gr_union),
    present_in_gr1 = overlap1,
    present_in_gr2 = overlap2
  )
  
  ggplot(result_df, aes(A = present_in_gr1, B = present_in_gr2)) + 
    geom_venn(auto_scale = TRUE, fill_color = c(row$color1, row$color2), fill_alpha = 0.8, show_percentage = FALSE, stroke_size = 0, text_size = 4, set_name_size = 4, set_names = c(row$celltype1, row$celltype2)) + 
    theme_void() + 
    coord_fixed()
  
  ## Get the proportion of celltype2 peaks that are shared: 
  df$cell2overlap[[i]] <- nrow(result_df[result_df$present_in_gr2==TRUE & result_df$present_in_gr1==TRUE,]) / nrow(result_df[result_df$present_in_gr2==TRUE,])
  
}

df$cell2overlap <- as.numeric(df$cell2overlap)
df$celltype1 <- factor(df$celltype1, levels = c("THP1", "HEK"))
ggplot(df, aes(x = celltype2, y = cell2overlap, fill = celltype1)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  scale_fill_manual(values = c("#ADD8E6", "#2E8B57")) + 
  theme_minimal() + 
  xlab("") + 
  ylab("Proportion of overlapping iMGL or primary macrophage peaks")


## Finally, read in our MPRA tested variants and see the overlap with each dataset: 
## Read in the allelic output from the MPRA_allelic.R script
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)
mpra <- GPos(seqnames = Rle(mpra$seqnames), 
             pos = mpra$pos, 
             rsid = mpra$rsid)

length(subsetByOverlaps(mpra, mac))/length(mpra) *100 ## 7.67%
length(subsetByOverlaps(mpra, mic))/length(mpra) *100 ## 3.72%
length(subsetByOverlaps(mpra, hek))/length(mpra) *100 ## 2.97%
