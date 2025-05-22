## Script to look at the enrichment of MPRA elements from THP1 and HEK 

## Load libraries
library(data.table)
library(dplyr)
library(GenomicRanges)
library(readxl)
library(rtracklayer)
library(tibble)
library(ggplot2)
library(pheatmap)


## Read in data: 
## Read in the allelic output from the MPRA_allelic.R script
adjusted <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt", header=T)

## Read in the activity output from the MPRA_active_lmer.R script
## Specifically we want the res_var object that describes activity at the variant level
active <- read.table("/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", header=T, sep = "\t")
active <- active[active$active==TRUE,]
active <- paste0(active$chr, "_", active$pos)

## Combine 
adjusted$active <- FALSE
adjusted[adjusted$id %in% active,]$active <- TRUE

## Define emVars
emvar <- adjusted[adjusted$sig==TRUE & adjusted$active == TRUE,] ## THP1 MPRA emvars from Bond 2024

## Compare with HEK cell MPRA data
hek_all = read_excel("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s1.xlsx", sheet = 2, skip=4) # allele
hek_act = read_excel("/work/users/m/a/marielle/external/cooper2022/science.abi8654_data_s5.xlsx", sheet = 2, skip=11) # activity

hek_AD = hek_all[hek_all$Disorder=="AD",]
hek_act = hek_act[hek_act$label.fdr=="active",]

hekadrange = GRanges(hek_AD$chr, IRanges(hek_AD$pos, hek_AD$pos))
mcols(hekadrange) = hek_AD[,c(1:4,7:ncol(hek_AD))]

hekactrange = GRanges(hek_act$chr, IRanges((hek_act$start+hek_act$end)/2,(hek_act$start+hek_act$end)/2))

olap = findOverlaps(hekadrange, hekactrange)
hek_allact = hekadrange[unique(queryHits(olap))] 
hek_emvar = data.frame(hek_allact[hek_allact$q<0.05,]) # 76 emVars according to their standards

## Comparison with Nott et al cell type specific H3K27ac data

bigwig_depth = fread("/proj/hyejunglab/HiC_resource/HiC_celltype_Nott/bigwig_depth.tab")
colnames(bigwig_depth) = c("chr","start","end","astro","neuro","oligo","micro")
bigwig_depth = bigwig_depth[!(bigwig_depth$chr %in% c("chrX","chrY","chrM")),]
bigwig_normfactor = colSums(bigwig_depth[,4:ncol(bigwig_depth)], na.rm=TRUE)

### EMVARS
## FOR MACROPHAGE
emvaranges = GRanges(unlist(lapply(strsplit(emvar$id,split="_"),'[[',1)), IRanges(as.numeric(unlist(lapply(strsplit(emvar$id,split="_"),'[[',2)))),
                     rsid=emvar$rsid)
## FOR HEK
emvaranges = GRanges(Rle(hek_emvar$seqnames), IRanges(start = hek_emvar$start, end = hek_emvar$end),
                     rsid=hek_emvar$rsID)


snpranges = resize(emvaranges, width = 500, fix = "center")

### ALL TESTED ELEMENTS
tested <- GRanges(seqnames = Rle(adjusted$seqnames),
                  ranges = IRanges(start = adjusted$pos, end = adjusted$pos+1), 
                  rsid = adjusted$rsid, 
                  padj = adjusted$adj.P.Val)
tested <- GRanges(seqnames = Rle(hek_AD$chr),
                  ranges = IRanges(start = hek_AD$pos, end = hek_AD$pos+1), 
                  rsid = hek_AD$rsID, 
                  padj = hek_AD$q)

snpranges <- resize(tested, width = 500, fix = "center")


### ALLELIC VARIANTS: 
tested <- tested[tested$padj<0.05,]
snpranges <- resize(tested, width = 500, fix = "center")

## Extract cell type information 
celltype = c("LHX2","NEUN","OLIG2","PU1")
for(i in 1:length(celltype)){
  celli = celltype[i]
  atac_import = import.bw(paste0("/proj/hyejunglab/HiC_resource/HiC_celltype_Nott/human_",celli,"nuclei_H3K27ac_epilepsy_pooled_hg19.ucsc.bigWig"),which=snpranges)
  olap = findOverlaps(atac_import, snpranges)
  olaprange = atac_import[queryHits(olap)]
  mcols(olaprange) = cbind(mcols(olaprange), mcols(snpranges[subjectHits(olap)]))
  olapdf = data.frame(olaprange)
  olapsummary = olapdf %>% group_by(rsid) %>% summarise(atac_score=sum(score*width))
  
  if(i==1){
    snp_atac=olapsummary
  }else{
    snp_atac = left_join(snp_atac, olapsummary, by="rsid")
  }
}
colnames(snp_atac) = c("snp",celltype)
snp_atac[is.na(snp_atac)] <- 0
snp_atac = data.frame(snp_atac)
rownames(snp_atac) = snp_atac$snp
snp_atac$LHX2 = snp_atac$LHX2/bigwig_normfactor[1]*100000
snp_atac$NEUN = snp_atac$NEUN/bigwig_normfactor[2]*100000
snp_atac$OLIG2 = snp_atac$OLIG2/bigwig_normfactor[3]*100000
snp_atac$PU1 = snp_atac$PU1/bigwig_normfactor[4]*100000

snpdf = data.frame(snpranges)
snp_cluster_atac = left_join(snpdf, snp_atac, by=c("rsid" = "snp"))
snp_cluster_atac[is.na(snp_cluster_atac)] <- 0

data <- snp_cluster_atac[,7:ncol(snp_cluster_atac)]
## run this for allelic variants 
data <- snp_cluster_atac[,8:ncol(snp_cluster_atac)]

## set colnames
colnames(data)<- c("Astrocyte", "Neuron", "Oligodendrocyte", "Microglia")

## remove rows where the max = 0
data$max <- apply(data, 1, max)
data <- data[data$max>0,]
data <- data[,1:4]

## remove rows where there's NA
data <- na.omit(data)

pheatmap(data,
         show_rownames=F,
         scale = "row",
         cluster_cols=FALSE, 
         cluster_rows=TRUE, 
         treeheight_row=0, 
         border_color = NA, 
         main = "HEK293T MPRA-allelic variants")

## save the data to look at to figure out issues with clustering
macdata <- data
hekdata <- data

macEnrich <- apply(macdata, 1, which.max) |> table() |> as.data.frame()
colnames(macEnrich) <- c("celltype", "frequency")
macEnrich$celltype <- colnames(data)
macEnrich$percentage <- macEnrich$frequency/nrow(macdata)
macEnrich$group <- "THP1"

hekEnrich <- apply(hekdata, 1, which.max) |> table() |> as.data.frame()
colnames(hekEnrich) <- c("celltype", "frequency")
hekEnrich$celltype <- colnames(data)
hekEnrich$percentage <- hekEnrich$frequency/nrow(hekdata)
hekEnrich$group <- "HEK"

## combine
enrichdata <- rbind(macEnrich, hekEnrich)
enrichdata$group <- factor(enrichdata$group, levels = c("THP1", "HEK"))
## change factors
enrichdata$celltype <- factor(enrichdata$celltype, levels = c("Astrocyte", "Neuron", "Oligodendrocyte", "Microglia"))
ggplot(enrichdata, aes(x = group, y = percentage, fill = celltype)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#e59f01", "#00a074", "#0072b1", "#ba3b34")) + 
  theme_minimal() + 
  xlab("") + 
  ylab("emVars with highest H3K27ac Signal (%)")

## fishers tests
a <- enrichdata[enrichdata$celltype=="Microglia" & enrichdata$group=="THP1",]$frequency
b <- enrichdata[!enrichdata$celltype=="Microglia" & enrichdata$group=="THP1",]$frequency |> sum()
c <- enrichdata[enrichdata$celltype=="Microglia" & enrichdata$group=="HEK",]$frequency
d <- enrichdata[!enrichdata$celltype=="Microglia" & enrichdata$group=="HEK",]$frequency |> sum()

counts <- matrix(c(a, b, c, d), 
                 nrow = 2, 
                 byrow = TRUE)
colnames(counts) <- c("Microglia", "Not Microglia")
rownames(counts) <- c("THP1", "HEK")

# Perform Fisher's Exact Test
fisher.test(counts)


## Also write a supplementary table that shows which emVars have what enrichment from the heatmap in figure 
data
data$rsid <- snp_cluster_atac$rsid
## add a column saying which cell type is the most enriched 
data$mostEnriched <- apply(macdata[,1:4], 1, which.max)
data[data$mostEnriched==1,]$mostEnriched <- "Astrocyte"
data[data$mostEnriched==2,]$mostEnriched <- "Neuron"
data[data$mostEnriched==3,]$mostEnriched <- "Oligodendrocyte"
data[data$mostEnriched==4,]$mostEnriched <- "Microglia"

macEnrich <- apply(macdata, 1, which.max) |> table() |> as.data.frame()

