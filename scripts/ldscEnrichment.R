### Script to perform LDSC regression analysis, scripts adapted from Hyejung Won
### Generates figures 2c, S3b

## List the paths to where you would like to run the ldsc analysis 
## For me, this is a directory for atac and k27 and within that, a different directory for each condition
atac <- list.files("/work/users/m/a/marielle/work/AD3D/data/ldsc/lima/atac", full.name=TRUE)
k27 <- list.files("/work/users/m/a/marielle/work/AD3D/data/ldsc/lima/k27", full.name=TRUE)

assays <- list(atac, k27)
names(assays) <- c("atac", "k27")

gwas <- "ad.jansen" ## name that will be used in the output files 

for (x in 1:length(assays)){
  sub <- assays[[x]]
  for (i in 1:length(sub)){
    bashout <- sub[[i]]
    ## for the first set of peaks to run, this script will prepare an sbatch command to be launched within a slurm system 
    
    ## header for the shell script to run 
    command1 <- paste0("echo '#!/bin/tcsh'> ",bashout, "/",  gwas, ".tcsh")
    ## adds the actual command that runs the ldscs.py script 
    ## specially uses gwas summary statisics, 1000G variant freq and baseline LD
    command2 <- paste0("echo '/proj/hyejunglab/program/ldsc/ldsc.py ",
                       "--h2 /proj/hyejunglab/crossdisorder/LDSC/sumstat/", gwas,".sumstats.gz ",
                       "--out ", bashout, "/", gwas, "_results ", # your result file location
                       "--frqfile-chr /proj/hyejunglab/program/ldsc/LDSC/1000G_Phase3_frq/1000G.EUR.QC. ",
                       "--overlap-annot --ref-ld-chr ", bashout, "/,/proj/hyejunglab/program/ldsc/LDSC/baselineLD_v1.1/baselineLD. ",
                       "--w-ld-chr /proj/hyejunglab/program/ldsc/LDSC/weights_hm3_no_hla/weights.",
                       "' >> " , bashout, "/", gwas,".tcsh")
    launch <- paste0("sbatch -n 1 --mem=40g -o ", bashout, "/",  gwas, ".out ", bashout, "/",  gwas,".tcsh")
    
    system(command1) #launch the command to generate the script
    system(command2) #launch the command to generate the script 
  }
}

### After this, you need to launch those shell scripts from the command line 
### Then within this directory you'll need to run the following per each chromosome: 
### python /proj/hyejunglab/program/ldsc/ldsc.py --l2 --bfile /proj/hyejunglab/program/ldsc/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr --ld-wind-cm 1 --annot dir/$chr.annot.gz --out dir/$chr --print-snps /proj/hyejunglab/program/ldsc/LDSC/list.txt
## You'll end up with a file that for my data is called ad.jansen_results.results, whicih is what you'll read in next

## Read in the output data and visualize
assays <- c("atac", "k27")
files <- 
  list.files(paste0("/work/users/m/a/marielle/work/AD3D/data/ldsc/lima/", assays), full.names = TRUE) |> 
  list.files("ad.jansen_results.results", full.names = TRUE)
atac <- files[grep("atac", files)]
k27 <- files[grep("k27", files)]

assayList <- list(atac, k27)
names(assayList) <- assays
cats <- c("resting", "activated")

data <- data.frame()
for (x in 1:length(assayList)){
  temp <- assayList[[x]]
  for (i in 1:length(temp)){
    df <- read.table(temp[[i]], header=T)
    df$group <- cats[i]
    df$assay <- assays[x]
    data <- rbind(data, df)
  }
}
data$group <- factor(data$group, levels = c("resting", "activated"))

## this data has lots of information, but specifically we want to just pull the rows for L2_0 
sub <- dplyr::filter(data, Category=="L2_0")

## Visualize! 
## Enrichment as seen in figure 2d
enrich <- ggplot(sub[sub$assay=="atac",], aes(x = group, y = Enrichment, color = group)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.05,
                position=position_dodge(.9)) + 
  ylim(0, 35) + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("AD Heritability Enrichment") + 
  scale_color_manual(values = c("#F0B0AA", "#EC6F61")) + 
  annotate("text", label = "p = 0.00033", x = 1, y = 30, size = 3) + 
  annotate("text", label = "p = 0.00382", x = 2, y = 35, size = 3)

## figure s3b
enrich <- ggplot(sub[sub$assay=="k27",], aes(x = group, y = Enrichment, color = group)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.05,
                position=position_dodge(.9)) + 
  ylim(0, 40) + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  xlab("") + 
  ylab("AD Heritability Enrichment") + 
  scale_color_manual(values = c("#F0B0AA", "#EC6F61")) + 
  annotate("text", label = "p = 0.00033", x = 1, y = 30, size = 3) + 
  annotate("text", label = "p = 0.00382", x = 2, y = 35, size = 3)
