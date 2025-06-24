## Script to look which gRNAs targetinggene promoters successfully knocked down their target gene, generates the figures S9j,k,l

## Load libraries
library(data.table)
library(edgeR)
library(ggplot2)
library(tidyverse)

## Read in CRISPR data
crispr <- fread("/work/users/m/a/marielle/work/AD3D/CRISPRi/DESeq2_all_results_200kbCutoff_NewPadj.csv")
counts <- fread("/proj/phanstiel_lab/Share/Share_Cochran/CombinedCounts.csv") |> 
  as.data.frame()
rownames(counts) <- counts$V1
counts <- counts[grep("Q", colnames(counts))]
counts <- cpm(counts) |> 
  as.data.frame()
meta <- fread("/work/users/m/a/marielle/work/AD3D/CRISPRi/samplesheetWithIPSCmicrogliaScores.csv")
## Get the control samples only 
controls <- meta[meta$Type=="Control",]$Sample
controls <- counts[colnames(counts) %in% controls]
controls$mean <- apply(controls, 1, mean)

## Which promoters are successfully downregulated?
proms <- crispr[crispr$Type=="Promoter",]$gRNA_name |> unique() ## 11 promoter gRNAs 

## For each of these promoters, loop through and just get the gene that it corresponds to 
data <- data.frame()
for (i in 1:length(proms)){
  sub <- crispr[crispr$gRNA_name==proms[i],]
  sub[sub$Gene==proms[i],]
}  

## Find a good way to visualize this: 
promKD <- 
  crispr[crispr$gRNA_name==crispr$Gene,] |> 
  select(gRNA_name, padj200kb, log2FoldChange) |>
  unique() |> 
  mutate(log10p = -log10(padj200kb), sign = sign(log2FoldChange), 
         signedP = sign * log10p) |> 
  arrange(signedP)
promKD$sig <- FALSE
promKD[promKD$padj200kb<0.05,]$sig <- TRUE
promKD$gRNA_name <- factor(promKD$gRNA_name, levels = promKD$gRNA_name)
promKD$sig <- factor(promKD$sig, levels = c(TRUE, FALSE))

## Visualize the pvalues for these gRNAs and their targets:
ggplot(promKD, aes(x = gRNA_name, y = signedP, color = sig, fill = sig)) + 
  geom_point() + 
  geom_col(width = 0.05) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45)) + 
  xlab("") + 
  ylab("Signed -log10(padj)") + 
  scale_fill_manual(values = c("#528B8B", "gray")) + 
  scale_color_manual(values = c("#528B8B", "gray")) + 
  geom_hline(yintercept = -1.3, lty = 3) + 
  geom_hline(yintercept = 1.3, lty = 3)

## No detected gene (3) #ADAMTS4,EXOC3L2,SPDYE3: aren't even expressed in iTF microglia 
## significant downreg (3) #APOE,BIN1,ZNF594
## not sig (4) #EPHA1,FCGR3B,HLA-DRB1,STAG3
## upregulated (1) #ZYX

## grab the baseMeans from the control samples for each of these genes 
ref <- data.frame(gene = c("APOE", "BIN1", "ZNF594", "EPHA1", "FCGR3B", "HLA-DRB1", "STAG3"), 
                  class = c("sig", "sig", "sig", "nonsig", "nonsig", "nonsig", "nonsig"))

data <- data.frame()
for (i in 1:nrow(ref)){
  sub <- controls[rownames(controls)==ref[i,1],]
  sub$class <- ref[i,2]
  data <- rbind(data, sub)
}

data$gene <- rownames(data)
data <- 
  data |> 
  pivot_longer(cols = starts_with("Q"), values_to = "CPM")
data$class <- factor(data$class, levels = c("sig", "nonsig"))
ggplot(data, aes(x = gene, y = log(CPM), fill = class)) + 
  geom_boxplot() + 
  geom_point(shape = 21) + 
  facet_wrap(~class, scales = "free_x") + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90)) + 
  xlab("") + 
  scale_fill_manual(values = c("#528B8B", "gray"))


pval <- wilcox.test(data[data$class=="sig",]$CPM, data[data$class=="nonsig",]$CPM)$p.value
ggplot(data, aes(x = log(CPM), fill = class)) + 
  geom_density(alpha = 0.5) + theme_minimal() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#528B8B", "gray")) + 
  annotate("text", x = 3, y = 0.3, label = paste0("p = ", round(pval, 10)), size = 3) + 
  ylab("Density")