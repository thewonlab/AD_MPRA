## Script to identify MPRA-active elemenets with a linear mixed model 

## Load libraries
library(reshape2)
library(lme4)
library(lmerTest)

### Read in barcode counts
### This is available as processed data at GEO GSE273887 (AD_MPRA_aggregated_counts.txt)
batchvar <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batchvar_N9N6_combined_withControls.txt", header=T, sep = "\t")

## Remove NAs and Infs
batchvar <- na.omit(batchvar)
repnum <- 15
remove <- which(is.infinite(rowSums(batchvar[,2:(repnum*4)])))
if(length(remove)>0){
  batchvar <- batchvar[-remove,]
}

## Set the repnum: 
repnum <- 9 # resting replicate number aka con
repnum2 <- 6 #LPS+INFg replicate number aka treat

## Since we have different numbers of repnums for each condition, split now into control and treated
con <- cbind(batchvar[,1], batchvar[grep("con", colnames(batchvar))], batchvar[,62:63])
trt <- cbind(batchvar[,1], batchvar[grep("treat", colnames(batchvar))], batchvar[,62:63])

## Select JUST the rna/dna ratios for both control and treated samples:
con <- con[seq(5, (ncol(con)), by = 4)]
trt <- trt[seq(5, (ncol(trt)), by = 4)]

### RUN THIS FOR JUST CONTROL
ratio <- con
ratio <- cbind(batchvar$variant, ratio, batchvar$SNP, batchvar$allele)
colnames(ratio) <- c("variant", paste0("S", 1:repnum), "SNP", "allele")
repnum <- 9

## RUN THIS FOR JUST TREATED
ratio <- trt
ratio <- cbind(batchvar$variant, ratio, batchvar$SNP, batchvar$allele)
colnames(ratio) <- c("variant", paste0("S", 1:repnum2), "SNP", "allele")
repnum <- 6

### RUN THIS TO COMBINE ALL SAMPLES 
## combine for one large dataframe with 15 samples, combined both treatments into one dataframe
ratio <- cbind(con, trt)
ratio <- cbind(batchvar$variant, ratio, batchvar$SNP, batchvar$allele)
repnum <- 15
colnames(ratio) <- c("variant", paste0("S", 1:repnum), "SNP", "allele")

## Separate negative controls from AD vars
neg <- ratio[(grep("Scrambled", ratio$variant)),]
tested <- ratio[(grep("chr", ratio$variant)),]

## Remove outliers from negative controls: 
threshold <- 2
neg_mean = c()
neg_median = c()
outlier <- c()
for(i in 1:nrow(neg)){
  neg_i = as.numeric(neg[i,2:(repnum+1)])
  neg_scale = scale(neg_i)
  neg_i = neg_i[abs(neg_scale)<threshold]
  if(length(neg_i)<repnum){
    outlier[i] <- TRUE
  } else {
    outlier[i] <- FALSE
  }
  neg_mean = c(neg_mean, mean(neg_i))
  neg_median = c(neg_median, median(neg_i))
}
neg$outlier <- outlier
## Remove outliers
neg <- neg[neg$outlier==FALSE,]
neg <- neg[,1:(ncol(neg)-1)]
neg$class <- "neg"
neg$median <- apply(neg[,2:(repnum+1)], 1, median)
neg$mean <- apply(neg[,2:(repnum+1)], 1, mean)
## Also get the median across 9 samples: 
negMedian <- apply(neg[,2:(repnum+1)], 2, median)
negMean <- apply(neg[,2:(repnum+1)], 2, mean)
## Combine neg and tested into one large dataframe: 
tested$class <- "tested"
#ratio <- rbind(neg, tested)

## Visualize distribution of negative controls vs tested variants: 
tested$median <- apply(tested[,2:(repnum+1)], 1, median)
tested$mean <- apply(tested[,2:(repnum+1)], 1, mean)

## Now for every single variant tested, we have 15 different values, one for each replicate 
## We need to figure out which variants are statistically higher than the negative controls 
## For every single variant, test it's 15 values vs the collapsed 15 negative controls in order to control for differences in sequenicng depth  
variants <- unique(tested$SNP)
res <- data.frame()
for (i in 1:nrow(tested)){
  row <- tested[i,]
  
  glm_input = row[,c(2:16)]
  glm_input = data.frame(rbind(glm_input, negMean))
  glm_input$var = c("var", "neg")
  
  glm_input = melt(glm_input)
  glm_input$condition = c(rep("ctrl", 18), rep("treated", 12))
  glm_input$pairs = rep(paste0("pair", 1:15), each=2)
  
  lme_result = summary(lmer(value ~ var*condition + (1 | pairs), data = glm_input))
  glm_result = lme_result$coefficients

  df <- 
    data.frame(variant = row$variant, 
               pvalue_activity = glm_result[2,5], 
               beta_activity = glm_result[2,1], # var"var" = in var column, "var"/neg
               pvalue_condition = glm_result[3,5],
               beta_condition = glm_result[3,1], # condition"treated" = in condition column, "treated"/ctrl
               pvalue_interaction = glm_result[4,5],
               beta_interaction = glm_result[4,1]) # treated/ctrl
  
  res <- rbind(res, df)
  print(i)
}

res_var = c()
for(i in 1:length(variants)){
  var2test = res[grep(variants[i], res$variant),]
  var2test = var2test[which(min(var2test$pvalue_activity)==var2test$pvalue_activity),]
  res_var = rbind(res_var, var2test)
}

### MPRA-ACTIVE ELEMENT DEFINITION 
res_var$fdr_activity = p.adjust(res_var$pvalue_activity, "BH")
dim(res_var[res_var$fdr_activity < 0.05 & res_var$beta_activity > 0, ]) # 556 MPRA-active elements

res_act = res_var[res_var$fdr_activity < 0.05 & res_var$beta_activity > 0, ]
res_act$var = unlist(lapply(strsplit(res_act$variant, split=":"), '[[', 1))

res$var = unlist(lapply(strsplit(res$variant, split=":"), '[[', 1))
### MPRA-ACTIVE CONDITIONAL VARIANTS 
res_cond = res[res$var %in% res_act$var, ] #needs to be active in at least one allele
res_cond$fdr_interaction = p.adjust(res_cond$pvalue_interaction, "BH")
res_cond = res_cond[res_cond$fdr_interaction<0.05, ]
dim(res_cond)
length(unique(res_cond[res_cond$beta_interaction>0, "var"])) # treatment-specific: 36
length(unique(res_cond[res_cond$beta_interaction<0, "var"])) # ctrl(resting)-specific: 173



#### DO NOT COPY TO GITHUB FROM HERE DOWN

load("/proj/hyejunglab/MPRA/RNAseq/AD/Hyejung/activity_lme_interaction.rda")


## remove any variants that are active in both because they were listed twice (different alleles) 
treatedOnly <- treatedOnly[!treatedOnly %in% restingOnly]
restingOnly <- restingOnly[!restingOnly %in% treatedOnly]

contextspecific <- unique(c(treatedOnly, restingOnly))
active <- res_act$var
shared <- active[!active %in% contextspecific]

# chr11_59848232_ad:T:C_C = a variant with treatment-specific activity
res_cond[res_cond$variant=="chr11_59848232_ad:T:C_C",]
#                     variant pvalue_activity beta_activity pvalue_condition beta_condition pvalue_interaction beta_interaction fdr_interaction
#173 chr11_59848232_ad:T:C_C    2.601147e-13      1.286767     0.0008537093      0.3444379       2.287565e-13         2.055041  2.536909e-10

# chr7_99950006_ad:T:C_C = a variant with condition-specific activity
res_cond[res_cond$variant=="chr7_99950006_ad:T:C_C",]
#                    variant pvalue_activity beta_activity pvalue_condition beta_condition pvalue_interaction beta_interaction fdr_interaction
#6839 chr7_99950006_ad:T:C_C    1.028692e-13      1.022792      0.008291092      0.3444379       7.887759e-10       -0.7981789  2.915841e-07

save(res_var, res, res_cond, file="/proj/hyejunglab/MPRA/RNAseq/AD/Hyejung/activity_lme_interaction.rda")


active <- unlist(lapply(strsplit(active, split="_ad"), '[[', 1))

# Load output of lmer 
load("/proj/hyejunglab/MPRA/RNAseq/AD/Hyejung/activity_lme_interaction.rda")
## Full allele-levels results with raw-pvalues 
res$var <- unlist(lapply(strsplit(res$variant, split="_ad"), '[[', 1))
res$allele <- unlist(lapply(strsplit(res$variant, split=":"), '[[', 3))
res$allele <- unlist(lapply(strsplit(res$allele, split="_"), '[[', 2))
res$chr <- unlist(lapply(strsplit(res$var, split="_"), '[[', 1))
res$pos <- unlist(lapply(strsplit(res$var, split="_"), '[[', 2))
res <- res[,c("chr", "pos", "allele", "pvalue_activity", "beta_activity", "pvalue_condition", "beta_condition", "pvalue_interaction", "beta_interaction")]

## variant level MPRA-active analysis 
res_var$var <- unlist(lapply(strsplit(res_var$variant, split="_ad"), '[[', 1))
res_var$allele <- unlist(lapply(strsplit(res_var$variant, split=":"), '[[', 3))
res_var$allele <- unlist(lapply(strsplit(res_var$allele, split="_"), '[[', 2))
res_var$chr <- unlist(lapply(strsplit(res_var$var, split="_"), '[[', 1))
res_var$pos <- unlist(lapply(strsplit(res_var$var, split="_"), '[[', 2))
res_var <- res_var[,c("chr", "pos", "allele", "pvalue_activity", "beta_activity", "fdr_activity", "pvalue_condition", "beta_condition", "pvalue_interaction", "beta_interaction")]
res_var$active <- FALSE
res_var[res_var$fdr_activity<0.05 & res_var$beta_activity>0,]$active <- TRUE

## condition level analysis 
res_cond <- res_cond[res_cond$fdr_interaction<0.05, ]
res_cond$var <- unlist(lapply(strsplit(res_cond$variant, split="_ad"), '[[', 1))
res_cond$allele <- unlist(lapply(strsplit(res_cond$variant, split=":"), '[[', 3))
res_cond$allele <- unlist(lapply(strsplit(res_cond$allele, split="_"), '[[', 2))
res_cond$chr <- unlist(lapply(strsplit(res_cond$var, split="_"), '[[', 1))
res_cond$pos <- unlist(lapply(strsplit(res_cond$var, split="_"), '[[', 2))
res_cond$condition <- FALSE
res_cond[res_cond$beta_interaction>0 & res_cond$fdr_interaction<0.05,]$condition <- "treated"
res_cond[res_cond$beta_interaction<0 & res_cond$fdr_interaction<0.05,]$condition <- "resting"

write.table(res, file = "/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_res.txt", row.names = FALSE, sep = "\t", quote = F)
write.table(res_var, file = "/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_resvar.txt", row.names = FALSE, sep = "\t", quote = F)
write.table(res_cond, file = "/work/users/m/a/marielle/work/AD3D/tables/MPRA_lmer_rescond.txt", row.names = FALSE, sep = "\t", quote = F)

conOnly <- res_cond[res_cond$condition=="resting",]
conOnly <- paste0(conOnly$chr, "_", conOnly$pos)

trtOnly <- res_cond[res_cond$condition=="treated",]
trtOnly <- paste0(trtOnly$chr, "_", trtOnly$pos)

## Prepare for homer: 
## Regular MPRA data: 
mpra <- read.table("/work/users/m/a/marielle/work/AD3D/data/mpra/batch3_N9N6_combined_allelic.txt")
mpra$class <- ""
mpra[mpra$id %in% conOnly,]$class <- "control"
mpra[mpra$id %in% trtOnly,]$class <- "treated"
mpra[mpra$class=="",]$class <- "background"
mpra$id <- paste0("variant", 1:nrow(mpra))

## Prepare each set for Homer: 
cats <- unique(mpra$class)

homer <- list()
for (i in 1:length(cats)){
  data <- mpra[mpra$class==cats[i],]
  data <- data[,c(1,10,11,12)]
  data$start <- data$pos-75
  data$end <- data$pos+75
  data$score <- 1 #for bed file
  #data <- data[,c(1,2,5,6,4)] #homer file
  data <- data[,c(2,5,6,4,7,1)] #bed file
  #colnames(data) <- c("id", "chrom", "start", "end", "strand") #homer file
  colnames(data) <- c("chrom", "start", "end", "strand", "score", "name") #bed file
  homer[[i]] <- data
}
names(homer) <- cats

## Write outputs: 
#write.table(homer$shared, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM/shared/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer$control, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM_FDR05/control/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer$treated, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM_FDR05/treated/input/input.txt", row.names = F, sep = "\t", quote = F)
write.table(homer$background, "/work/users/m/a/marielle/work/AD3D/data/mpra/homer/diffActiveVars_GLM_FDR05/background.txt", row.names = F, sep = "\t", quote = F)




