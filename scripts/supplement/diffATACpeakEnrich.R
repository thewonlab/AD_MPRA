## Regular ATAC motif enrichment on just differential ATAC peaks
## Generates figure s3c-d

## Read in homer output 
treated <- read.delim("/work/users/m/a/marielle/worpk/AD3D/LIMA/atac/homer/gained/output/knownResults.txt")
control <- read.delim("/work/users/m/a/marielle/work/AD3D/LIMA/atac/homer/lost/output/knownResults.txt")

motifs <- list(control, treated)
names <- c("control", "treated")

## Clean up each motif file (new colnames, remove % signs)
newcolnames <- c("MotifName", "Consensus", "P.value", "Log P.value", "q value", 
                 "Number target sequences", "Percent_target_sequences", 
                 "Number background sequences", "Percent_background_sequences")

abridge <- function(longname){
  temp <- strsplit(longname, '[(]')
  shortname <- unlist(temp)[1]
}

for (i in 1:length(motifs)){
  
  #Subset by the motif list
  motif <- motifs[[i]]
  #Assign new colnames
  colnames(motif) <- newcolnames
  #Assign abreviated names
  motif$MotifName <- sapply(1:nrow(motif), function(x) abridge((as.character(motif[,1]))[x]))
  #Remove % signs from percent target sequences and percent background sequences
  motif$Percent_target_sequences <- as.numeric(sub("%", "", motif$Percent_target_sequences))
  motif$Percent_background_sequences <- as.numeric(sub("%", "", motif$Percent_background_sequences))
  #Calculate log2 enrichment and log10pval
  motif$log2enrichment <- log2(motif$Percent_target_sequences/motif$Percent_background_sequences)
  motif$log10pval <- -log10(exp(motif$`Log P.value`))
  motif$class <- names[i]
  #Reassign to the list 
  motifs[[i]] <- motif
}

## Gained Motifs
data <- motifs[[1]] 
data <- data[1:30,]
max <- max(data[!data$log10pval=="Inf",]$log10pval)
data$MotifName <- factor(data$MotifName, levels = rev(data$MotifName))
atacmotifsCon <- 
  ggplot(data, aes(x = MotifName, y = log10pval)) + 
  geom_bar(stat = "identity", fill = "#FFA6A6", alpha = 0.75) + 
  theme_classic() + coord_flip() + 
  labs(title = "Resting") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 8)) + 
  xlab("") + ylab("-log10(p-val)") + 
  geom_text(label = data$MotifName, hjust = "left", y = 0.2, color = "white", size = 2)


## Lost Motifs 
data <- motifs[[2]] 
data <- data[1:30,]
max <- max(data[!data$log10pval=="Inf",]$log10pval)
data[data$log10pval=="Inf",]$log10pval <- max
data$MotifName <- factor(data$MotifName, levels = rev(data$MotifName))
data$MotifName <- factor(data$MotifName, levels = rev(data$MotifName))
atacmotifsTrt <- 
  ggplot(data, aes(x = MotifName, y = log10pval)) + 
  geom_bar(stat = "identity", fill = "#E04F4F", alpha = 0.75) + 
  theme_classic() + coord_flip() + 
  labs(title = "Activated") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.title = element_text(size = 10), 
        axis.text.x = element_text(size = 6), 
        axis.title.x = element_text(size = 8)) + 
  xlab("") + ylab("-log10(p-val)") + 
  geom_text(label = data$MotifName, hjust = "left", y = 0.2, color = "white", size = 2)

