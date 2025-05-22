### Script to identify emvars from all tested variants
### An emvar must be MPRA-active and MPRA-allelic 

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
emvar <- adjusted[adjusted$sig==TRUE & adjusted$active == TRUE,]$rsid


