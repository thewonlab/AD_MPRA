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

## REFERENCE

Manuscript in preparation.
