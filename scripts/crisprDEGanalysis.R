## Script to perform differential gene expression analysis from the CRISPRi RNA-seq data
## This script generates the count matrix and generates the figures 6b-c and Sup 9a-c

## Load libraries:
library(data.table)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(ggrepel)
library(rtracklayer)
library(DESeq2)
library(GenomicRanges)

#Read in list of file paths and convert them into a list
files=as.list(readLines("/cluster/projects/ADFTD/Phanstiel_Microglia/outs/path_to_counts.tsv"))

#Get sample name from each file path
#assumes the sample name is the 8th field when the path is split by "/"
n=sapply(strsplit(unlist(files),"[/]"),"[[",8)

#Reach each HTSeq count file and store it in a list
counts<-lapply(files, function(x){
  tab<-read.table(x,header=F)
})
#Assign sample names to the "counts" list elements
names(counts)<-n

#assign column names to each count data frame
for(i in 1:length(counts)){
  colnames(counts[[i]])<-c("gene_id", names(counts)[i])
}

#combining all counts into a single matrix
mat<-do.call(cbind, counts)
#set row names of the matrix to the "gene_id" column
rownames(mat)<-mat[,1]
#remove every second column starting from the frist (gene_id columns), leaves only the count columns in matrix
#assumes that gene_id columns are in odd positions (1,3,5,...,263)
mat2<-mat[,-seq(1,264,by=2)]

#add gene names
#Import GTF file containing gene annotations
gtf<-import("/cluster/home/aanderson/myers/BrainTF/RNA/STAR/gencode.v42.primary_assembly.annotation.gtf")
#Filter FTD to only include entries of the type "gene"
gtf<-gtf[gtf$type=="gene"]
#extract gene ID and gene name columns from the metadata columns in columns 5 and 7
key<-mcols(protein_coding)[,c(5,7)]

#Convert count matrix to data frame for merging
mat2<-as.data.frame(mat2)
#add gene_id as new column to df
mat2$gene_id<-rownames(mat2)
#Merge count data frame with gene annotations by "gene_id"
mat3<-merge(mat2, key, by="gene_id")
#make row names of the merged data frame the gene name
rownames(mat3)<-mat3$gene_name

#data cleaning
#remove gene_id and gene_name columns from data frame
mat3<-mat3[,-c(1,118)]
#change sample names in columns to match sample sheet make QJ1.QJ1 -> QJ1
colnames(mat3)<-sapply(strsplit(colnames(mat3),".",fixed = TRUE),`[`,1)

#read in sample sheet with metadata
samplesheet<- read.csv("/cluster/projects/ADFTD/Phanstiel_Microglia/MicrogliaCRISPRiSampleSheet.csv",header = TRUE)
#filter to only include samples that are present in the count matrix
samplesheet<-samplesheet[which(samplesheet$Sample %in% colnames(mat3)),]

#reorder matrix to same order as sample sheet
mat3<-mat3[,samplesheet$Sample]

#Calculate CPM for combined and cleaned count matrix
cpm<-cpm(mat3)

#Filter lowly expressed genes out
cpm2<-cpm[rowSums(cpm)>1,]

#Perform PCA on transposed CPM matrix
#PCA expected variable in columns and observations in rows
#Scaling standardizes the variable to have unit variance before PCA
pca<-prcomp(t(cpm2), scale=T)

#Extract first two principal components & convert to data frame
pca_df<-as.data.frame(pca$x[,1:2])

#Add sample identifiers to the PCA data frame
pca_df$sample<-colnames(mat3)

#Integrate sample metadata into PCA, assigns samples where controls, promoters, and variant regions were targeted
pca_df$target<-samplesheet$Type

#convert mat3 back to matrix if not already
mat3<- as.matrix(mat3)

#Add the total number of reads per sample to the PCA data frame
pca_df$NumReads<-colSums(mat3)

#Integrate Marker Gene Expression into PCA Data Frame
#iPSC marker genes
pca_df$LIN28A<-as.numeric(as.list(cpm["LIN28A",])) #avg cpm 2.382261
pca_df$POU5F1<-as.numeric(as.list(cpm["POU5F1",])) #avg cpm  0.6802256
pca_df$ESRG<-as.numeric(as.list(cpm["ESRG",])) #avg cpm 0
pca_df$NANOG<-as.numeric(as.list(cpm["NANOG",])) #avg cpm 0.02066019
pca_df$SOX2<-as.numeric(as.list(cpm["SOX2",]))#avg cpm 0
pca_df$NODAL<- as.numeric(as.list(cpm["NODAL",])) #avg cpm 0.2829955
pca_df$FOXD3<- as.numeric(as.list(cpm["FOXD3",])) #avg cpm 0.01289008

#Compile iPSC marker genes into separate matrix for downstream analysis
#exclude ESRG and SOX2 because 0 cpm values
iPSCgenes<-cpm[c("LIN28A","POU5F1","NANOG","NODAL","FOXD3"),]

#Microglia marker genes
pca_df$HEXB<-as.numeric(as.list(cpm["HEXB",])) #highly expressed avg cpm 514.8478
pca_df$SELPLG<-as.numeric(as.list(cpm["SELPLG",])) #expressed avg cpm 19.7335
pca_df$CD74<-as.numeric(as.list(cpm["CD74",])) #avg cpm 1305.647
pca_df$CSF1R<-as.numeric(as.list(cpm["CSF1R",])) #avg cpm 1046.967
pca_df$CD14<-as.numeric(as.list(cpm["CD14",])) #avg cpm 250.7919
pca_df$CD40<-as.numeric(as.list(cpm["CD40",]))#avg cpm 846.4192
pca_df$CD86<-as.numeric(as.list(cpm["CD86",])) #avg cpm  21.75886
pca_df$CD163<-as.numeric(as.list(cpm["CD163",])) #avg cpm 35.94298
pca_df$C3<-as.numeric(as.list(cpm["C3",])) #avg cpm 411.484
pca_df$AIF1<-as.numeric(as.list(cpm["AIF1",])) #avg cpm 137.401
pca_df$CX3CR1<-as.numeric(as.list(cpm["CX3CR1",])) #expressed but avg of 9.61664
pca_df$A2M<-as.numeric(as.list(cpm["A2M",])) #avg cpm 465.0429
pca_df$GAS6<-as.numeric(as.list(cpm["GAS6",])) #avg cpm 17.42895

#Compile iPSC marker genes into separate matrix for downstream analysis
MicrogliaGenes<-cpm[c("HEXB","SELPLG", "CD74", "CSF1R","CD14", "CD40", "CD86", "CD163", "C3",  "AIF1", "CX3CR1","A2M", "GAS6"),]


#Perform t tests of cpm control vs test regions and plot 
# Read in the data
counts <- read.csv("/cluster/projects/ADFTD/Phanstiel_Microglia/CombinedCounts.csv", row.names = 1)
samplesheet <- read.csv("/cluster/projects/ADFTD/Phanstiel_Microglia/MicrogliaCRISPRiSampleSheet.csv")

# Calculate Counts Per Million (CPM)
cpm_values <- cpm(counts)

# Remove sample QJ77 from CPM and samplesheet
cpm_values <- cpm_values[, -77]
samplesheet <- samplesheet[-77, ]

# Extract unique genes, excluding the 12th gene
genes <- unique(unlist(strsplit(samplesheet$TargetGene, ",")))
genes <- genes[-12] #removes safe harbor

# Loop through each gene
for (i in 1:length(genes)) {
  
  # Filter samples: Control or those targeting the current gene
  plot_keep <- samplesheet[samplesheet$Type == "Control" | grepl(genes[i], samplesheet$TargetGene),]
  
  # Add gene expression values to the dataframe
  plot_keep$gene <- cpm_values[genes[i], plot_keep$Sample]
  
  # Create a new 'Group' variable:
  plot_keep$Group <- ifelse(plot_keep$Type == "Control", "Control",
                            ifelse(plot_keep$Type == "Promoter", "Promoter",
                                   plot_keep$gRNA_name))
  
  # Define the order of 'Group': Control, Promoter, then Variants
  variant_names <- unique(plot_keep$Group[plot_keep$Group != "Control" & plot_keep$Group != "Promoter"])
  plot_keep$Group <- factor(plot_keep$Group, levels = c("Control", "Promoter", sort(variant_names)))
  
  # Ensure 'Type' is a factor with desired levels for consistent coloring
  plot_keep$Type <- factor(plot_keep$Type, levels = c("Control", "Promoter", "Variant"))
  
  # Define comparisons:
  # - Control vs Promoter
  # - Control vs each Variant
  comparisons <- list(c("Control", "Promoter"))
  
  if (length(variant_names) > 0) {
    # Add comparisons for each Variant vs Control
    variant_comparisons <- lapply(variant_names, function(x) c("Control", x))
    comparisons <- c(comparisons, variant_comparisons)
  }
  
  # Define the PDF output file path
  pdf_file <- paste0(genes[i], "_boxplot.pdf")
  
  # Open the PDF device
  pdf(pdf_file, width = 8, height = 6) # Increased width for multiple variants
  
  # Generate the boxplot
  p <- ggplot(plot_keep, aes(x = Group, y = gene, fill = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Removed outliers for clarity
    geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) + # Added jitter for data points
    theme_classic() +
    scale_fill_manual(values = c("grey30", "darkslategray1", "darkslategray4")) +
    stat_compare_means(comparisons = comparisons,
                       method = "t.test",
                       label = "p.format",
                       hide.ns = TRUE) + # Hide non-significant p-values
    labs(title = paste("Expression of", genes[i]),
         x = "Sample Type",
         y = "CPM") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(size = 12)
    )
  
  # Print the plot to the PDF
  print(p)
  
  # Close the PDF device
  dev.off()
}

#Read gtf from Ensembl Genes 109
biomart_tab = read.table("/cluster/home/jloupe/R/BrainTF_updated/mart_export.txt",sep="\t", header = TRUE,check.names=FALSE)

#Rename columns for consistency and easier reference using dplyr
biomart_tab <- biomart_tab %>% rename(
  "Gene stable ID"="gene_id",
  "Gene name"="gene_name",
  "Transcript length (including UTRs and CDS)"="length" ,
  "RefSeq match transcript (MANE Select)"="MANE", 
  "Transcript support level (TSL)"="TSL")
#Sort gtf by gene_name in ascending order
biomart_tab <- biomart_tab[order(biomart_tab$gene_name),]
#Sort by MAN in descending order to prioritize MAN select transcripts
biomart_tab <- biomart_tab[order(biomart_tab$MANE,decreasing = TRUE),]
#Remove duplicated gene names, keeping first occurrence (which has MANE select transcripts prioritized)
biomart_tab <- biomart_tab[!duplicated(biomart_tab$gene_name),]

#use matrix of counts BEFORE any trimming. no cpm transformation.
matX <- mat3
#add column "gene_name" to "matX" using the row names, which correspond to gene identifiers
matX$gene_name <- rownames(matX)

#merge the transcript length basedd on gene name
gxc <- inner_join(matX, biomart_tab[,c("gene_name","length")], by="gene_name")

#Make new df with gene_name and length columns 
genes <- gxc[,c("gene_name","length")]

#get just raw counts for TPM calculation
countperM <- gxc %>% select(-c(gene_name,length))

#Define a function to calculate TPM from raw counts and gene lengths by dividing counts by gene lengths, normalize by the sum of each sample, and scale by 1e6
tpm <- function(countperM, lengths) {
  rate <- countperM / lengths
  rate / sum(rate) * 1e6
}

#Apply tpm function to each sample in the count matrix
tpms <- apply(countperM, 2, function(x) tpm(x, genes$length))

#Convert to data frame and round all TPM values to one decimal place
tpms <- as.data.frame(tpms) %>%
  mutate_all(~ round(., 1))

#assign the gene names as the row names of TPM df
rownames(tpms) <- gxc$gene_name

#Save as csv 
write.csv(tpms,"/cluster/projects/ADFTD/Phanstiel_Microglia/CRISPRi_TPMcounts.csv")

#Select cell type marker genes
MicrogliaTPM<-tpms[c("HEXB","SELPLG", "CD74", "CSF1R","CD14", "CD40", "CD86", "CD163", "C3",  "AIF1", "CX3CR1","A2M", "GAS6"),]
iPSCtpm<-tpms[c("LIN28A","POU5F1","NANOG","NODAL","FOXD3"),]

#Log10 Transform TPM for plotting
iPSC_log10tpm<-log10(iPSCtpm +1)
Microlgia_log10tpm<-log10(MicrogliaTPM +1)

#convert log10 transformed tpms to df and orer the rows alphabetically by gene name
iPSC_log10tpm_df<-as.data.frame(iPSC_log10tpm)
iPSC_log10tpm_df<-iPSC_log10tpm_df[order(rownames(iPSC_log10tpm_df)),]
Microglia_log10tpm_df<- as.data.frame(Microlgia_log10tpm)
Microglia_log10tpm_df<-Microglia_log10tpm_df[order(rownames(Microglia_log10tpm_df)),]

#Combine iPSC and Microglia log10  tpm data frames
combinedTPM<- rbind(iPSC_log10tpm_df,Microglia_log10tpm_df)

#Convert to matrix and remove sample QJ77 (failed QC --> not enough reads)
combinedTPM_mat<-as.matrix(combinedTPM[,-77])

#Label each gene as either "iPSC" or "Microglia" based on origin
gene_groups<-c(
  rep("iPSC",nrow(iPSC_log10tpm_df)),
  rep("Microglia",nrow(Microglia_log10tpm_df))
)

#Create a row annotation to indicate the gene type (iPSC or Microglia with distinct colors)
row_anno<-rowAnnotation(
  GeneType = gene_groups,
  col = list(
    GeneType = c(
      "iPSC" = "#B2DF8A", 
      "Microglia" = "#33A02C"
    )
  ),
  show_annotation_name = TRUE
)

#Create a column annotaiton to indicate if sample is targeting Control, Promoter, or Variant Region
column_anno<-HeatmapAnnotation(type=samplesheet$Type)

pdf("/cluster/projects/ADFTD/Phanstiel_Microglia/DifferentiationMarkersHeatMap_TPM_removeQJ77_removeTPM0Genes.pdf", width=6, height=6)
Heatmap(
  combinedTPM_mat,
  name = "log10(TPM)",
  top_annotation = column_anno,
  left_annotation = row_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = gene_groups,
  col = colorRamp2(c(min(combinedTPM_mat),max(combinedTPM_mat)),c("blue","red"))
)
dev.off()

#Read in count matrix, first row contains column names, and set first column as row names 
mat3<- read.csv("/cluster/projects/ADFTD/Phanstiel_Microglia/CombinedCounts.csv", header = TRUE, row.names = 1)

#Filter out lowly expressed genes and sample QJ77 that failed QC
cpm<-cpm(mat3)
mat3_filt<-mat3[rowMeans(cpm)>0.5,-77]

#samplesheet
samplesheet<-read.csv("/cluster/projects/ADFTD/Phanstiel_Microglia/MicrogliaCRISPRiSampleSheet.csv", header = TRUE)
#Also remove QJ77 from samplesheet
samplesheet<-samplesheet[-77,]

#Create an empty list to store DESeq2 results for each gRNA
all_results<-list()

#remove Control gRNA from list (comparing Control to Promoter and Variant regions)
keep<-samplesheet[which(samplesheet$Type != "Control"),]$gRNA_name

#Perform DESeq2 analysis comparing each gRNA to the Safe Harbor control 
for (i in keep) {
  dds<-DESeqDataSetFromMatrix(countData=mat3_filt[,which(samplesheet$gRNA_name %in% c(i,"SafeHarbor"))], 
                              colData=samplesheet[which(samplesheet$gRNA_name %in% c(i,"SafeHarbor")),], 
                              design= ~ gRNA_name + Microglia_score + IPSC_score)
  dds<-DESeq(dds) #fit for each test
  results <- results(dds, contrast=c("gRNA_name",i,"SafeHarbor"))
  results$gRNA_name=i
  results$Gene<-rownames(results)
  all_results[[i]]=results
}


#combine together into one large data frame
names(all_results)<-NULL
all_results<-do.call(rbind,all_results)
#filter to just regions with adjusted p value less than 0.05
sig_results<-all_results[which(all_results$padj <= 0.05),]

#Volcano plot of significant differentially expressed genes within 200kb of target
#read in sample metadata
samplesheet<-read.csv("MicrogliaCRISPRiSampleSheet.csv", header = TRUE)

#Create a GenomicRanges object from the first 108 samples in samplesheet to remove safe harbor controls
samplesheetGR<-GRanges(samplesheet[1:108,])

#Import gene annotation from GTF
gtf<-import("/cluster/home/aanderson/myers/BrainTF/RNA/STAR/gencode.v42.primary_assembly.annotation.gtf")
#Filter to include only entries of the type start_codon
gtf<-gtf[gtf$type=="start_codon"]
#filter to only include protein-coding genes with the CCDS tag
gtf_gene<-gtf[which(gtf$gene_type=="protein_coding" & gtf$tag=="CCDS"),]

#Create a GenomicRegions object containing gene coordinates and names
gene_coords <- GRanges(
  seqnames = seqnames(gtf_gene),
  ranges   = ranges(gtf_gene),
  strand   = strand(gtf_gene),
  gene_name = mcols(gtf_gene)$gene_name
)

#Assign gene names as names of the object to allow lookup by name
names(gene_coords) <- mcols(gtf_gene)$gene_name  

#Define a function to determine if a gene is within 200kb of a target region
within_200kb <- function(region, gene_id, gene_coord_map) {
  idx <- which(names(gene_coord_map) == gene_id)
  if (length(idx) == 0) return(FALSE)
  
  # Use only the first matching gene
  gene_gr <- gene_coord_map[idx[1]]
  
  # Check same chromosome
  if (as.character(seqnames(region)) != as.character(seqnames(gene_gr))) return(FALSE)
  
  # Compute distance
  dist <- distance(region, gene_gr)
  return(!is.na(dist) && dist <= 200000)
}

#Convert DESeq2 results to data frame
all_results_df<-as.data.frame(all_results)

#Merge the sample metadata wtih DESeq2 results based on gRNA_name
merge<-merge(samplesheet, all_results_df, by = "gRNA_name")
#Convert to GenomicRanges object
merge<-GRanges(merge)

#Subset gene_coords to include only genes that overlap with the regions within a max gap of 500 kb
genes_keep<-subsetByOverlaps(gene_coords,merge,maxgap = 500000)

#subset to only include DESeq2 results for genes that are kept (within 200kb of region)
merge_sub<-merge[which(merge$Gene %in% genes_keep$gene_name),]

#apply within_200kb function to each region-gene pair to ensure proximity within 200kb
keep_idx <- vapply(seq_along(merge_sub), function(i) {
  within_200kb(merge_sub[i], merge_sub$Gene[i], gene_coords)
}, logical(1))

#Subset merge_sub to include only those genes that are within 200kb of the region
merge200kb<-merge_sub[keep_idx]

#Create a new column that indicates if the gene is the same as the predicted target gene
merge200kb$link<-merge200kb$Gene==merge200kb$TargetGene
#Create a column color2 based on significance
merge200kb$color2<-factor(ifelse(merge200kb$padj<0.05,merge200kb$gRNA_name, 0))
#Create a unique identifier for each gene-guideRNA pair to remove duplicates later
merge200kb$uniq<-paste0(merge200kb$padj, "-", merge200kb$Gene)
#Remove duplicated entries based on 'uniq' column
merge200kb<-merge200kb[!(duplicated(merge200kb$uniq)),]
#Create a new column 'gene_guide' by concatenating 'Gene' and 'gRNA_name' with a hypen
merge200kb$gene_guide<-paste0(merge200kb$Gene, "-",merge200kb$gRNA_name)

#Save results as a csv file
write.csv(merge200kb, "DESeq2allResults_DiffScore_CPMcutoff_fitWithinTest_genesWithin200kb.csv")

#Create volcano plot using ggplot2
pdf("VolcanoPlot_DiffScore_CPMcutoff_fitWithinTest.pdf", width=4, height=4)
ggplot(as.data.frame(merge200kb), aes(x=log2FoldChange, y= -log10(padj), color=as.factor(color2), shape=Type))+
  theme_classic()+
  ylab("-log10(FDR)")+
  geom_vline(xintercept=0, lty="dashed",alpha=0.5,color="grey")+ geom_hline(yintercept= -log10(0.05), lty="dashed",alpha=0.5,color="grey")+
  theme(legend.position="none")+  
  ggtitle("CRISPRi")+
  scale_shape_manual(values=c(1,19))+
  scale_alpha_manual(values=c(0.7,1))+
  geom_text_repel(data=as.data.frame(merge200kb[which(merge200kb$padj<0.05),]), size=4,color="grey30",aes(label=gene_guide), min.segment.length = unit(0, 'lines'))+
  geom_point()
dev.off()

#calculate new padj with only genes +- 200 kb
gRNAs<-unique(merge200kb$gRNA_name)
merge200kb$padj200kb<-NA
for(i in gRNAs){
  subset_df<-merge200kb[which(merge200kb$gRNA_name==i),]
  subset_df$padj200kb<-p.adjust(subset_df$pvalue, method = "BH")
  merge200kb<-rbind(merge200kb, subset_df)
}
merge200kb<-merge200kb[which(is.na(merge200kb$padj200kb)==FALSE),]
head(merge200kb)
write.csv(merge200kb, "DESeq2_all_results_200kbCutoff_NewPadj.csv")

merge200kb$color2<-factor(ifelse(merge200kb$padj200kb<0.01,merge200kb$gRNA_name, 0))
merge200kb$uniq<-paste0(merge200kb$padj200kb, "-", merge200kb$Gene)
merge200kb<-merge200kb[!(duplicated(merge200kb$uniq)),]
merge200kb$gene_guide<-paste0(merge200kb$Gene, "-",merge200kb$gRNA_name)

pdf("VolcanoPlot_200kb_newpadj.pdf", width=4, height=4)
ggplot(as.data.frame(merge200kb), aes(x=log2FoldChange, y= -log10(padj200kb), color=as.factor(color2), shape=Type))+theme_classic()+ylab("-log10(FDR)")+
  geom_vline(xintercept=0, lty="dashed",alpha=0.5,color="grey")+ geom_hline(yintercept= -log10(0.05), lty="dashed",alpha=0.5,color="grey")+
  theme(legend.position="none")+  ggtitle("CRISPRi")+
  scale_shape_manual(values=c(1,19))+
  scale_alpha_manual(values=c(0.7,1))+
  geom_text_repel(data=as.data.frame(merge200kb[which(merge200kb$padj200kb<0.01),]), size=4,color="grey30",aes(label=gene_guide), min.segment.length = unit(0, 'lines'))+geom_point()
dev.off()

merge200kbsig<-merge200kb[which(merge200kb$padj200kb < 0.05),]
head(merge200kbsig)
table(merge200kbsig$gRNA_name)
#      APOE       BIN1 rs12703526 rs12721109   rs184017  rs6979218 rs72838287 
#         1          1          1          3          2          1          1 
#    SPDYE3      STAG3     ZNF594        ZYX 
#         1          2          1          1 

merge200kbsigUP<-merge200kbsig[which(merge200kbsig$log2FoldChange > 0),]
table(merge200kbsigUP$gRNA_name)
#rs12703526 rs12721109      STAG3        ZYX 
#         1          1          1          1 
merge200kbsigDOWN<-merge200kbsig[which(merge200kbsig$log2FoldChange < 0),]
table(merge200kbsigDOWN$gRNA_name)
#    APOE       BIN1 rs12721109   rs184017  rs6979218 rs72838287     SPDYE3 
#         1          1          2          2          1          1          1 
#     STAG3     ZNF594 
#         1          1 

write.csv(merge200kbsig, "DESeq2Results_newPadj_200kbCutoff_sig.csv")
