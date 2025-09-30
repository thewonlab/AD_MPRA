## Script to generate the dot/line plots and boxplots for the CRISPR analysis as shown in Figure 6e-g and Sup 9d-i

## Load libraries
library(readxl)
library(data.table)
library(cowplot)
library(edgeR)
library(ggrepel)
library(ggpubr)

## Read in data
crispr <- fread("/work/users/m/a/marielle/work/AD3D/CRISPRi/DESeq2_all_results_200kbCutoff_NewPadj.csv")
counts <- fread("/proj/phanstiel_lab/Share/Share_Cochran/CombinedCounts.csv") |> 
  as.data.frame()
rownames(counts) <- counts$V1
counts <- counts[grep("Q", colnames(counts))]
counts <- cpm(counts) |> 
  as.data.frame()
meta <- fread("/work/users/m/a/marielle/work/AD3D/CRISPRi/samplesheetWithIPSCmicrogliaScores.csv")

## Visualize the genes that are affected by each guideRNA 
guides <- unique(crispr$gRNA_name)
guides <- guides[grep("rs", guides)]

data <- data.frame()
plots <- list()
for (i in 1:length(guides)){
  g <- guides[i]
  sub <- crispr[crispr$gRNA_name==g,]
  sub$sig <- FALSE
  sub[sub$padj200kb<0.05,]$sig <- TRUE
  sub <- 
    sub |> 
    select(gRNA_name, log2FoldChange, padj200kb, Gene, sig) |> 
    unique()
  sub$signedP <- -log10(sub$padj200kb) * sign(sub$log2FoldChange)
  ## add data to df 
  data <- rbind(data, sub)
  ## get a range for the y axis centered around zero: 
  max <- max(abs(sub$signedP))
  ## order by the abs signedP value: 
  sub <- sub[order(sub$signedP),]
  sub$Gene <- factor(sub$Gene, levels = sub$Gene)
  ## look at gene expression of nearby genes
  plot <- 
    ggplot(sub, aes(x = Gene, y = signedP, color = sig, fill = sig)) + 
    geom_col(width = 0.005) + 
    geom_point() + 
    theme_minimal() + 
    scale_color_manual(values = c("gray", "#528B8B")) + 
    scale_fill_manual(values = c("gray", "#528B8B")) + 
    geom_hline(yintercept = 0, lty = 3) + 
    theme(legend.position = "none", 
          axis.text.x = element_text(angle=90, size = 6),
          axis.text.y = element_text(size = 7), 
          axis.title.y = element_text(size = 8),
          panel.grid.major.y = element_line(linetype = "dashed"), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_blank()) + 
    xlab("") + ylab("Signed -log10(padj)") + 
    ylim(-max, max)
  plots[[i]] <- plot
  
}

## variants to include in the main figure: 
main <- c("rs12721109", "rs72838287", "rs6979218", "rs184017")
sigplots <- plots[guides %in% main]
plot_grid(sigplots[[1]], sigplots[[4]], sigplots[[2]], sigplots[[3]], nrow = 1)

## nonsig for the supplement
nonsigplots <- plots[!guides %in% main]
plot_grid(nonsigplots[[1]], nonsigplots[[2]], nonsigplots[[3]], nonsigplots[[4]], nonsigplots[[5]],
            nonsigplots[[6]], nonsigplots[[7]], nonsigplots[[8]], nonsigplots[[9]], nonsigplots[[10]], 
            nrow = 2)

## Boxplots for the main figure: 
dfMain <- 
  list(data.frame(sample = c("QJ66", "QJ77", "QJ68", "QJ70", "QJ71", "QJ72", "QJ73", "QJ74", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("BIN1 Promoter", 4), rep("rs72838287", 4), rep("Control", 8)), 
                  gene = "BIN1", 
                  class = "main"), 
       data.frame(sample = c("QJ40", "QJ41", "QJ42", "QJ43", "QJ44", "QJ45", "QJ46", "QJ47", "QJ48", "QJ49", "QJ50", "QJ51", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("APOE Promoter", 8), rep("rs12721109", 4), rep("Control", 8)), 
                  gene = "APOE", 
                  class = "main"), 
       data.frame(sample = c("QJ40", "QJ41", "QJ42", "QJ43", "QJ44", "QJ45", "QJ46", "QJ47", "QJ48", "QJ49", "QJ50", "QJ51", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("APOE Promoter", 8), rep("rs12721109", 4), rep("Control", 8)), 
                  gene = "APOC1", 
                  class = "main"), 
       data.frame(sample = c("QJ56", "QJ57", "QJ58", "QJ59", "QJ60", "QJ97", "QJ98", "QJ99", "QJ100", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("STAG3 Promoter", 5), rep("rs6979218", 4), rep("Control", 8)), 
                  gene = "GPC2", 
                  class = "main"))

## for the supplement: 
dfSup <- 
  list(data.frame(sample = c("QJ75", "QJ76", "QJ77", "QJ82", "QJ83", "QJ84", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("EPHA1 Promoter", 3), rep("rs12703526", 3), rep("Control", 8)), 
                  gene = "EPHA1", 
                  class = "sup"), 
       data.frame(sample = c("QJ56", "QJ57", "QJ58", "QJ59", "QJ60", "QJ61", "QJ62", "QJ63", "QJ64", "QJ65", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("STAG3 Promoter", 5), rep("rs13230744", 5), rep("Control", 8)), 
                  gene = "STAG3", 
                  class = "sup"), 
       data.frame(sample = c("QJ15", "QJ16", "QJ17", "QJ18", "QJ19", "QJ20", "QJ21", "QJ22", "QJ23", "QJ24", "QJ25", "QJ26", "QJ27", "QJ28", "QJ29", "QJ30", "QJ31", "QJ32", "QJ33", "QJ34", "QJ35", "QJ36", "QJ37", "QJ38", "QJ39", 
                             "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("HLA-DRB1 Promoter", 5), rep("rs9271176", 4), rep("rs9270887", 5), rep("rs28732235", 3), rep("rs9271607", 4), rep("rs72847948", 4), rep("Control", 8)), 
                  gene = "HLA-DRB1", 
                  class = "sup"), 
       data.frame(sample = c("QJ40", "QJ41", "QJ42", "QJ43", "QJ44", "QJ45", "QJ46", "QJ47", "QJ52", "QJ53", "QJ54", "QJ55", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("APOE Promoter", 8), rep("rs184017", 4), rep("Control", 8)), 
                  gene = "TOMM40", 
                  class = "sup"), 
       data.frame(sample = c("QJ85", "QJ86", "QJ87", "QJ88", "QJ89", "QJ90", "QJ91", "QJ92", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("ZNF594 Promoter", 4), rep("rs77972827", 4), rep("Control", 8)), 
                  gene = "ZNF594", 
                  class = "sup"), 
       data.frame(sample = c("QJ78", "QJ79", "QJ80", "QJ81", "QJ82", "QJ83", "QJ84", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("ZYX Promoter", 4), rep("rs12703526", 3), rep("Control", 8)), 
                  gene = "ZYX", 
                  class = "sup"),
       data.frame(sample = c("QJ7", "QJ8", "QJ9", "QJ10", "QJ11", "QJ12", "QJ13", "QJ14", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("FCGR3B Promoter", 4), rep("rs4575098", 4), rep("Control", 8)), 
                  gene = "FCGR3B", 
                  class = "sup"), 
       data.frame(sample = c("QJ101", "QJ102", "QJ103", "QJ104", "QJ105", "QJ106", "QJ107", "QJ108", "QJ109", "QJ110", "QJ111", "QJ112", "QJ115", "QJ116", "QJ117", "QJ118"), 
                  type = c(rep("EXOC3L2 Promoter", 4), rep("rs6509201", 4), rep("Control", 8)), 
                  gene = "ERCC2", 
                  class = "sup")
  )

## combine into one df
dfs <- c(dfMain, dfSup)

plots <- list()
for (i in 1:length(dfs)){
  df <- dfs[[i]]
  sub <- counts[colnames(counts) %in% df$sample]
  sub <- 
    sub[rownames(sub) %in% df$gene,] |> 
    pivot_longer(cols = everything(), names_to = "sample", values_to = "CPM") |> 
    left_join(df, by = "sample")
  n <- length(unique(sub$type))
  cats <- unique(sub$type)
  ## figure out the appropriate order of the types
  sub$type <- factor(sub$type, levels = c(cats[grep("Control", cats)], cats[grep("Promoter", cats)], cats[grep("rs", cats)]))
  
  if (n==3){
    comp <- list(c("Control", cats[grep("Promoter", cats)]), 
                 c("Control", cats[grep("rs", cats)]))
    plots[[i]] <- 
      ggplot(sub, aes(x = type, y = CPM, fill = type)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 21) + 
      theme_minimal() + 
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=90, size = 7),
            axis.text.y = element_text(size = 7), 
            axis.title.y = element_text(size = 8)) + 
      xlab("") + ylab("Counts Per Million (CPM)") + 
      scale_fill_manual(values = c("#4E4E4E", "#A9DEE9", "#538B8B")) + 
      stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 2)
  } else if (n>3){
    comp <- list(c("Control", cats[grep("Promoter", cats)]), 
                 c("Control", cats[grep("rs", cats)][1]), 
                 c("Control", cats[grep("rs", cats)][2]), 
                 c("Control", cats[grep("rs", cats)][3]), 
                 c("Control", cats[grep("rs", cats)][4]), 
                 c("Control", cats[grep("rs", cats)][5]))
    plots[[i]] <- 
      ggplot(sub, aes(x = type, y = CPM, fill = type)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 21) + 
      theme_minimal() + 
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=90, size = 7),
            axis.text.y = element_text(size = 7), 
            axis.title.y = element_text(size = 8)) + 
      xlab("") + ylab("Counts Per Million (CPM)") + 
      scale_fill_manual(values = c("#4E4E4E", "#A9DEE9", rep("#538B8B", length(cats)-2))) + 
      stat_compare_means(comparisons = comp, method = "t.test", label = "p.format", hide.ns = FALSE, size = 2)
  }
}

genes <- c()
for (i in 1:length(dfs)){
  sub <- dfs[[i]]
  genes[i] <- unique(sub$gene)
}

## genes to include in the main figure: 
main <- c("BIN1", "APOE", "APOC1", "GPC2")
sigplots <- plots[genes %in% main]
plot_grid(sigplots[[1]], sigplots[[4]], sigplots[[2]], sigplots[[3]], nrow = 1)

## nonsig for the supplement
nonsigplots <- plots[!genes %in% main]
plot_grid(nonsigplots[[1]], nonsigplots[[2]], nonsigplots[[3]], nonsigplots[[4]], nonsigplots[[5]],
          nonsigplots[[6]], nonsigplots[[7]], nonsigplots[[8]], 
          nrow = 2)
