### Script to identify MPRA-allelic variants 

### Load libraries
library(mpra)

### Read in barcode counts
### This is available as processed data at GEO GSE273887 (AD_MPRA_aggregated_counts.txt)
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch_N9N6_combined_batchvar.txt", header=T)

### Read in variant-barcode map
### This is available as processed data at GEO GSE273887 (AD_MPRA_Barcode_Map.txt)
variant <- read.table("/work/users/m/a/marielle/proj/AD/MPRA/bcmap/fastq/variant_statistics.txt", header=T, sep="\t")

### Number of replicates
repnum <- 9 ## for resting
repnum2 <- 6 ## for proinflammatory

# dna matrix - all DNA counts for each variant for each allele
dnamat <- cbind(batchvar[seq(1,nrow(batchvar)-1,2),seq(3,4*repnum, by = 4)], batchvar[seq(2,nrow(batchvar),2),seq(3,4*repnum, by = 4)],
                batchvar[seq(1,nrow(batchvar)-1,2),seq(repnum*4+3, (repnum*4+3)*3, by = 4)[1:repnum2]], batchvar[seq(2,nrow(batchvar),2),seq(repnum*4+3, (repnum*4+3)*3, by = 4)[1:repnum2]])
rownames(dnamat) <- unique(batchvar$SNP)
colnames(dnamat) <- c(paste(rep("alt",repnum),paste0("s",1:repnum),sep="_"), paste(rep("ref",repnum),paste0("s",1:repnum),sep="_"),
                     paste(rep("alt",repnum2),paste0("s",1:repnum2),sep="_"), paste(rep("ref",repnum2),paste0("s",1:repnum2),sep="_"))

# rna matrix - all RNA counts for each variant for each allele
rnamat <- cbind(batchvar[seq(1,nrow(batchvar)-1,2),seq(2,4*repnum, by = 4)], batchvar[seq(2,nrow(batchvar),2),seq(2,4*repnum, by = 4)],
                batchvar[seq(1,nrow(batchvar)-1,2),seq(repnum*4+2, (repnum*4+2)*3, by = 4)[1:repnum2]], batchvar[seq(2,nrow(batchvar),2),seq(repnum*4+2, (repnum*4+2)*3, by = 4)[1:repnum2]])
rownames(rnamat) <- unique(batchvar$SNP)
colnames(rnamat) <- c(paste(rep("alt",repnum),paste0("s",1:repnum),sep="_"), paste(rep("ref",repnum),paste0("s",1:repnum),sep="_"),
                     paste(rep("alt",repnum2),paste0("s",1:repnum2),sep="_"), paste(rep("ref",repnum2),paste0("s",1:repnum2),sep="_"))

# variant IDs (rsIDs)
varid <- rownames(dnamat) # list of variant IDs (chr_position)

# variant sequences
variant <- variant[,1:2]
variant$rsid <- paste0(unlist(lapply(strsplit(variant$name, split="_"),'[[',1)), "_", unlist(lapply(strsplit(variant$name, split="_"),'[[',2)))
varseq <- variant[match(varid, variant$rsid), "variant"]

## Combine into mpraset
mpraset <- MPRASet(DNA=dnamat, RNA=rnamat, eid=varid, eseq=varseq, barcode=NULL)

## Make design matrix based off of replicate numbers: 
samplenum <- length(colnames(mpraset))
samples <- colnames(mpraset)
intcpt <- rep(1, samplenum)
alt <- unlist(strsplit(samples, "_"))[seq(1, 2*samplenum, by = 2)]
alt[alt=="alt"] <- TRUE
alt[alt=="ref"] <- FALSE
control <- c(rep(TRUE, 18), rep(FALSE, 12)) #when batches are using diff number of rep

## Number of samples: 
samples <- as.numeric(unlist(strsplit(samples, "_s"))[seq(2,(length(samples)*2), by=2)]) 

## Design matrix: 
designmat <- 
  data.frame(intcpt = intcpt, 
             alt = as.logical(alt),
             rep = samples,
             control = control)

designmat <- model.matrix(~alt+control+alt:control)

# analyze
mpraresult <- mpralm(object=mpraset, 
                    design=designmat, 
                    plot=F,
                    aggregate="none",
                    normalize=T, 
                    block=samples,
                    model_type="corr_groups")

allelic <- topTable(mpraresult, coef = 2, number = Inf)
conditional <- topTable(mpraresult, coef = 4, number = Inf)

allelic[allelic$adj.P.Val<0.05,] |> nrow() ## Report the number of MPRA-allelic variants
conditional[conditional$adj.P.Val<0.05,] |> nrow() ## Report the number of MPRA-allelic Conditional variants 

## Downstream processing of actual MPRA variants: 
allelic$sign <- sign(allelic$logFC)
allelic$sig <- ""
allelic[allelic$adj.P.Val<0.05,]$sig <- TRUE
allelic[allelic$adj.P.Val>0.05,]$sig <- FALSE

write.table(allelic, file = "/proj/hyejunglab/MPRA/RNAseq/AD/Marielle/mpralm_allelic_results.txt", sep = "\t", quote = F, row.names = FALSE)
write.table(conditional, file = "/proj/hyejunglab/MPRA/RNAseq/AD/Marielle/mpralm_conditional_results.txt", sep = "\t", quote = F, row.names = FALSE)
save(mpraresult, file = "/proj/hyejunglab/MPRA/RNAseq/AD/Marielle/mpraresult.rda")


# Adjust LFC according to risk/protection (from the beta value from the summary statistics) 
beta <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/betaValues.txt", header=T) ## this is an excerpt from the Jansen et al 2019 summary statistics
allelic$id <- rownames(allelic)
adjusted <- merge(allelic, beta, by = "id")
#adjusting for direction of risk/protective
adjusted <-
  adjusted %>%
  mutate(corrected_logFC = ifelse(beta < 0, -1 * (logFC), logFC))


