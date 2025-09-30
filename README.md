# AD_MPRA

Created by Marielle Bond (marielle_bond@med.unc.edu) 7/30/2024

## DOWNLOADING DATA

- Download files files from _GEO_ (GSE273887) to your local/virtual environment
   
## REPOSITORY ORGANIZATION

- `scripts/`
  - `processing/` — scripts relating to processing MPRA and other related data
  - `figure1/` — scripts that perfom analysis related to Figure 1
  - `figure2/` — scripts that perfom analysis related to Figure 2
  - `figure3/` — scripts that perfom analysis related to Figure 3
  - `figure4/` — scripts that perfom analysis related to Figure 4
  - `figure5/` — scripts that perfom analysis related to Figure 5
  - `figure6/` — scripts that perfom analysis related to Figure 6
  - `supplement/` — scripts that perfom analysis related to the supplementary figures

## DESCRIPTION OF SCRIPTS

`scripts/processing/`
- MPRA processing:
   - **MPRA_active_lmer.R:** a script that should be run on _AD_MPRA_aggregated_counts.txt_ from GEO to generate MPRA-active elements using the _lme4_ R package
   - **MPRA_allelic.R:** a script that should be run on _AD_MPRA_aggregated_counts.txt_ from GEO to generate MPRA-allelic variants using the _mpralm_ R package
   - **MPRA_active_lmer.R:** a script that should be run on the outputs of _MPRA_allelic.R_ and _MPRA_active_lmer.R_ to generate the 47 emVars described in this study
 
- Other processing:
   - **epiRMextractParallel.R:** a script that was written to extract the signal from ATAC-seq _bigwig_ files from the **Roadmap Epigenomics** project. This script was specifically written to extract variant coordinates from _bigwig_ files in parallel from each cell type using a _SLURM_ job scheduler
   - **prepmotifBreakR.R:** a script that runs the _motifbreakR_ function from the _motifbreakR_ R package on each variant tested in this MPRA. This script was specifically written to run this function in parallel on each variant using a _SLURM_ job scheduler 

`scripts/figure1/`
- Description: These scripts all relate to characterizing MPRA-active and MPRA-allelic variants, and intersecting them to identify 47 emVars. 
   - **compareMPRAactive.R:** this script uses the output from _MPRA_active_lmer.R_ and _AD_MPRA_aggregated_counts.txt_ from _GEO_ GSE273887 to look at the difference in MPRA activity between MPRA-active elements, MPRA-inactive elements, and both positive and negative controls. It generates the violin plots shown in Figure 1b and Supplementary Figure S2a
   - **compareMPRAactiveMotifEnrich.R:** this script uses the output from _MPRA_active_lmer.R_ to prepare MPRA-active and MPRA-inactive elements to run the _homer_ function _findMotifs.pl_ from the command line. The script then visualizes the output as both a density plot shown in Figure 1c and a bar plot shown in Supplementary Figure S2c.
   - **compareMPRAactiveMotifs_conditional.R:**  this script is similar to _compareMPRAactiveMotifEnrich.R_, as it takes the context-specific MPRA-active elements from _MPRA_active_lmer.R_, prepares these two sets of variants (control- and inflammatory-specific), and prepares them to run the _homer_ function _findMotifs.pl_ from the command line. This script then generates the bar plots shown in Figure 1d.
   - **figure1ExampleLoci.R:** this script uses the GWAS summary statistics from Jansen et al. 2019 and the MPRA barcode counts from this study (_AD_MPRA_aggregated_counts.txt_ from _GEO_ GSE273887) to visualize the GWAS/MPRA data as shown in Figures 1e/h.  
   - **compareMPRAallelic_volcano.R:** this script uses the outputs from _MPRA_active_lmer.R_ and _MPRA_allelic.R_ to generate a volcano plot shown in Figure 1f that shows the log2(risk/protective) of variants, highlighting MPRA-allelic variants, and emVars (which are MPRA-allelic variants that are also MPRA-active).
   - **variantsPerLocus.R:** this script uses the outputs from _MPRA_active_lmer.R_ and _MPRA_allelic.R_ to visualize the distribution of all tested MPRA-elements across the GWAS loci from Jansen et al. 2019. The result is visualized as a barplot in Figure 1g. 

`scripts/figure2/`
- Description: These scripts all describe the multi-omic data from Reed et al 2022. Cell Reports that was used in this study. 
   - **diffOmics.R:** this script uses the supplementary material from Reed et al 2022 Cell Reports and generates the pie charts shown in Figure 2a that show the number of differential loops, ATAC peaks, H3K27ac peaks, and genes. 
   - **exampleDiffLoop.R:** this script uses the processed data from Reed et al 2022 Cell Reports GEO superseries GSE201376 to visualize the locus shown in Figure 2b using the _plotgardener_ package.
   - **ldscEnrichment.R:**  this script prepares the ATAC/H3K27ac data from Reed et al 2022 for stratified LDSC regression analysis (see Gazal et al. Nat Genetics 2017). The resulting figures are shown in Figure 2c and S3a. 
   - **rna_vs_scMicroglia.R:** this script uses the scRNA-seq cluster marker genes from Sun et al 2023 and intersects with the RNA-seq from Reed et al 2022 to generate the boxplots shown in Figure 2d. 
   - **motifEnrichmentPlusExpression.R:** this script uses ATAC-seq and RNA-seq data from Reed et al 2022 to perform TF motif enrichment with the _homer_ function _findMotifsGenome.pl_. For selected TF motifs, it shows the expression of the genes related to that motif. These results are shown in Figure 2e and Supplementary Figure S3c.
   - **diffABCpairs_andGO.R:** this script uses the ABC pairs for resting and proinflammatory macrophages (this study, Supplementary Table S5) and identifies pairs that are specific to either condition. For the genes at ABC pairs unique to each condition, it uses the _gost_ function from _gprofiler2_ to identify enriched GO terms. This script generates the venn diagram from Figure 2g and the bar plots for GO terms in Figure 2h. 
   - **ABCgeneCorrelation.R:** this script uses the ABC pairs for resting and proinflammatory macrophages (this study, Supplementary Table S5) to identify the correlation between genes at ABC pairs and changes in ABC scores. It generates the density plot from Figure 2g. 

`scripts/figure3/`
- Description: The scripts in this directory are used in analyes pertaining to the characterization of MPRA variants in the epigenomic context of microglia and macrophages. 
   - **dipScoreandATAC.R:** this script calculates the ATAC-seq signal and H3K27ac dip scores for MPRA-tested variants in our macrophage data. It generates Figure 3a and Supplemental Figure S4a-d.
   - **epiRMheatmap.R:** this script uses the output from the processing script _epiRMextractParallel.R_ to visualize the ATAC-signal in our macrophage data and the Roadmap Epigenomics ATAC-seq data from other cell types. It generates the figure shown in Figure 3b. 
   - **motifBreakRenrichment.R:**  this script uses the output from _prepmotifBreakR.R_ to identify emVars that often disrupt TF motifs, it generates Figure 3d. 
   - **plotBrokenMotif.R:** this script uses _plotgardener_ to visualize one of the motifs that was prioritized from motifBreakR shown in Figure 3e. 
   - **brainEnhancerEnrichment.R:** this script intersects emVars from our study and frVars from Cooper et al 2022 (Science) with the microglia, astrocyte, neuron, and oligodendocyte data from Nott et al 2019 H3K27ac data. It generates the heatmaps shown in Figure 3f and Supplemental Figure S3d-e.
   - **emVarsOlapQTL:** this script intersects our emVarsa nd the frVars from Cooper et al 2022 (Science) with the microglia eQTLs from Kosoy et al 2023. It generates Figure 3g. 

`scripts/figure4/`
- Description: Thes scripts are used to map emVars to target genes and visualize examples.
   - **eQTLabcTargetGenes.R:** this script identifies the possible target genes of emVars from this study by intersecting with our macrophage ABC data (Figure 2, this study) and microglia eQTLs from Kosoy et al 2023. It characterizes these variant-gene pairs and generates the figures shown in Figure 4a-d
   - **adTargetEnrich.R:** this script performs gene enrichment analysis for emVar target genes and various gene sets from previous studies (Hassleman et al 2019, Sun et al 2023, Gazestani et al 2023, Li et al 2023.) The results are shown in Figure 4e.
   - **figure4Examples.R:**  this script uses _plotgardener_ to visualize 2 variants that were mapped to a target gene using either ABC or eQTL data, shown in Figure 4f-g. 

`scripts/figure5/`
- Description: These scripts are used to describe our case study of rs9270887 at the HLA locus
   - **emsaQuantification.R:** this script uses quantifications of EMSA data and performs statistical analyses, shown in Figure 5d. 

`scripts/figure6/`
- Description: These are the scripts relating to our CRISPRi screen and use our RNA-seq data. 
   - **crisprDEGanalysis.R:** this script performs differental gene expression analysis using _DESeq2_ on our CRISPRi data. It also calculates iPSC and Miroglia scores for each sample, which is used in the _DESeq2_ design for differential analysis. It generates Figure 6b-c and Supplemental Figure S9a-c. 
   - **crisprABC.R:** this script intersects the CRISPRi RNA-seq data with ABC data (Supplemental Table 5, this study) to look at the relationship between RNA-seq log2(fold-change) in variant-ABC predicted pairs vs random ABC pairs). The resulting barplot is shown in Figure 6d. 
   - **crisprVariantGeneVisualization.R:**  this script generates the visualizations shown in Figure 6e-f and Supplemental Figure 9j-i. It uses the _DESeq2_ results to plot the log2(fold-change) of genes within 200kb of the gRNA, shows normalized counts from _DESeq2_ for the gene in samples that received the negative control gRNAs, promoter-targeting gRNAs, and variant-targeting gRNAs. 

`scripts/supplement/`
- Description: These are scripts that generate miscellaneous supplemental figure panels that were not included in all of the scripts described above.
   - **barcodeDistribution.R:** this script generates a histogram showing the distribution of barcodes per variant in our MPRA library, shown in Supplemental Figure S1c. 
   - **mpraReproducibility.R:** this script calculates the reproducibility between MPRA replicates in both our control and proinflammatory states, shown in Supplemental Figure S1d-e. 
   - **distanceToGWASlead.R:**  this script calculates the distance between all MPRA elements (stratified by active/inactive/allelic/emVar) to the lead SNP from their locus, shown in Supplemental Figure S2d.
   - **compareToGWAS.R:** this script compares MPRA active vs inactive and MPRA allelic vs nonallelic variants according to their GWAS beta and GWAS p-values, shown in Supplemental Figure S2e-h. 
   - **diffATACPeakEnrich.R:** this script uses differential ATAC-seq peaks from Reed et al 2022 (Cell Reports) and prepares them for TF motif enrichment with _homer_, then visualizes the output as shown in Supplemental Figure S3b. 
   - **compareToHEKmpra.R:**  this script compares the MPRA variants in our study to the MPRA variants tested in Cooper et al 2022 (Science), and generates the figures shown in Supplementary Figure S5a-e.
   - **compareEpigeneticProfiles.R:** this script intersects MPRA variants from our study and the MPRA variatns from Cooper et al 2022 (Science), and compares to the epigenetic data from macrophages (Reed et al 2022), microglia (Yang et al 2023), HEK293T (Dong et al 2024), and primary macrophages (Alasoo et al 2018). It generates the panels shown in Supplemtary Figure 5f-j. 
   - **extendedTargetGeneMapping.R:** this script builds on the gene mapping performed by _eQTLabcTargetGenes.R_ for Figure 4 by adding on two additional datasets: the microglia ABC pairs from Kosoy et al 2023 and the Macrophage response eQTLs from Alasoo et al 2018. The results are shown in Supplementary Figure S6d-e
   - **crisprPromoterExpression.R:** this script analyzes the knockdown of target genes and the promoter gRNA that we designed for our CRISPRi experiments, the results are shown in Supplementary Figure S9j-l. 

## REFERENCE

Manuscript in revision.
Medrxiv : https://www.medrxiv.org/content/10.1101/2024.09.13.24313654v1
