# AD_MPRA

Created by Marielle Bond (marielle_bond@med.unc.edu) 7/30/2024

## BEFORE RUNNING SCRIPTS

1. Download files files from _GEO_ (GSE273887)
2. Git clone or download our files and scripts to your local/virtual environment.
   
## RUNNING SCRIPTS

1. Generate MPRA results:
    - run _MPRA_allelic.R_ on _AD_MPRA_aggregated_counts.txt_ from GEO to generate MPRA-allelic variants using _mpralm_
    - run _MPRA_active.R_ on _AD_MPRA_aggregated_counts.txt_ from GEO if you want to generate MPRA-active elements
    - run _MPRA_emVar.R_ on the outputs of the above two scripts to generate emVars
2. Perform downstream analysis:
    - nearly all of the results use the outputs from the _MPRA_allelic.R_, _MPRA_active.R_ and _MPRA_emVar.R_, scripts generated above
    - the remainder of the scripts in the _scripts_ directory are used to perform downstream analysis and generate the figures shown in the manuscript
    - if a script was used to generate any figures, it will state it in the header of the R script
    - if a script requires the downloading of external data to perform the analysis, it is stated within that script as well as the methods section of the manuscript pertianing to that analysis
  
   

## REFERENCE

Manuscript in revision.
Medrxiv : https://www.medrxiv.org/content/10.1101/2024.09.13.24313654v1
