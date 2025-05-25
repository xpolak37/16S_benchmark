suppressMessages(suppressWarnings({
  library(data.table)
  library(ccrepe)
  library(igraph)
  library(readr)
  library(cowplot)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(phyloseq)
  library(MicrobiotaProcess)
  library(ggpubr)
  library(ggrepel)
  library(ggplotify)
  #library(radEmu)
  library(vegan)
  library(reshape2)
  library(pheatmap)
  library(mgcv)
  library(robustlmm)
  library(lmerTest)
  library(emmeans)
  library(magrittr)
  library(openxlsx)
  library(caret)
  library(MicrobiomeStat)
  library(glmnet)
  library(pROC)
  library(purrr)
  #library(umap)
  library(Maaslin2)
  library(ggvenn)
  library(ranger)
  library(doParallel)
  library(gbm)
  library(tidyr)
  library(kableExtra)
  library(tidyverse)
  library(picante)
  library(tidyr)
  library(mice)
}))

pca_plot_custom <- function(asv_table,taxa_table,metadata, 
                            measure="robust.aitchison",
                            variable = "Taxonomic_classifier",
                            size=3,
                            legend=TRUE){
  
  pca_asv_df <- asv_table %>% dplyr::select(-Color)
  pca_taxa_df <- taxa_table
  pca_metadata <- metadata
  
  cols <- colnames(pca_asv_df)
  rows <- rownames(pca_metadata)
  pca_asv_df <- pca_asv_df[c(cols[-2], cols[2])]
  pca_metadata <-  pca_metadata[c(rows[-1],rows[1]),]
  rownames(pca_metadata) <- NULL

  pc <- construct_phyloseq(pca_asv_df,pca_taxa_df,pca_metadata)
  
  ord_b <- ordinate(pc, method = "PCoA", distance = "bray")
  
  data_for_pca_bray <- as.data.frame(t(pc@otu_table)) %>%
    dplyr::mutate(SettingID= factor(pca_metadata$SettingID,levels=unique(pca_metadata$SettingID)),
                  Tool=factor(pca_metadata$process_tool,levels=unique(pca_metadata$process_tool)),
                  Replicate=factor(pca_metadata$Sample,levels=unique(pca_metadata$Sample)),
                  Taxonomic_classifier=factor(pca_metadata$taxonomy,levels=unique(pca_metadata$taxonomy)),
                  Preprocessing=factor(pca_metadata$prep_number,levels=unique(pca_metadata$prep_number)))
  
  imp_bray <- ord_b$values$Relative_eig
  pca_bray <- ord_b$vectors
  
  if (variable=="Taxonomic_classifier"){
    colors <- c(
      "#1E75B0", "#FA7D0E", "#7D7D7D", "#D22627","#1EFC1E"
    )
  }
  else if (variable=="Replicate"){
    colors =c("#9165B9", "#FAF33E","#3D52D5","#1EFC1E")
  } else if (variable=="Preprocessing"){
    colors=c(
      "#1E75B0", "#ABC3E3", "#FA7D0E","#FAB776",
      "#2B9D2B", "#D90368", "#D22627","#FA9593",
      "#9165B9", "#FAF33E", "#89544A","#C09991",
      "#DF75BE", "#F2B2CE","#1EFC1E")
  } else if (variable=="Tool"){
    colors=c("#89544A","#C09991","#1E75B0", "#D90368","#1EFC1E")
  }

  p_pcoa <- ggplot(data_for_pca_bray, 
                   aes(x=pca_bray[,1],y=pca_bray[,2], 
                       col= !!sym(variable))) +
    geom_point(show.legend =TRUE,size=size) +
    xlab(paste("PCo1 ", "(",round(imp_bray[1]*100,2),"%", ")", sep=""))+
    ylab(paste("PCo2 ", "(",round(imp_bray[2]*100,2),"%", ")", sep=""))+
    theme_bw() +
    scale_color_manual(values=colors) + 
    theme(
    axis.text.x = element_text(size = 17,color="black"),
    axis.text.y = element_text(size = 17,color="black"),
    axis.title.x = element_text(size = 17,color="black"),
    axis.title.y= element_text(size = 17,color="black"))
  
  if (!legend){
    p_pcoa <- p_pcoa + 
      theme(legend.position="none")
  }
  return(p_pcoa)
}


construct_phyloseq <- function(asv_table, taxa_table, metadata){
  # This function constructs phyloseq object and is mainly used in functions
  # with methods that require phyloseq object 
  
  otu_mat <- asv_table
  tax_mat <- taxa_table
  samples_df <- metadata
  
  # removing rownames
  rownames(otu_mat) <- NULL
  rownames(tax_mat) <- NULL
  rownames(samples_df) <- NULL
  
  # set SeqID as rownames
  otu_mat <- otu_mat %>%
    tibble::column_to_rownames("SeqID") 
  tax_mat <- tax_mat %>% 
    tibble::column_to_rownames("SeqID")
  
  # if more than one samples is provided
  if (nrow(samples_df)>1){
    samples_df <- as.data.frame(samples_df) %>% 
      tibble::column_to_rownames("SampleID") 
  } else{
    # if only one sample - rownames have to be set manually
    rownames(samples_df) <- samples_df$SampleID
    samples_df <- samples_df[, -which(names(samples_df) == "SampleID")]
  }
  
  # construct matrix
  otu_mat <- as.matrix(otu_mat)
  tax_mat <- as.matrix(tax_mat)
  
  # preprocessing for phyloseq object
  OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(tax_mat) 
  samples = sample_data(samples_df)
  
  # resulting object
  object <- phyloseq(OTU, TAX, samples)
  return(object)
}

bray_custom <- function(asv_table,metadata){
  boxplot_df <- metadata %>%
    dplyr::mutate(prep_and_tool= paste(prep_number,process_tool,sep="_"))
  
  boxplot_df[boxplot_df$SampleID=="Mock_reference","prep_and_tool"] <- "Mock_reference"
  mock_b_dist <- as.data.frame(as.matrix(
    vegdist(t(asv_table %>% dplyr::select(-SeqID,-Color)), method="bray")))
  samples <- rownames(mock_b_dist)
  mock_b_dist <- mock_b_dist %>% dplyr::select(Mock_reference)
  
  rownames(mock_b_dist) <- samples
  colnames(mock_b_dist) <- "Bray_Curtis"

  boxplot_df <- merge(mock_b_dist %>% rownames_to_column("SampleID"),
                      boxplot_df,by="SampleID") 
  
  which_row <- which(boxplot_df$SampleID == "Mock_reference")
  rows <- rownames(boxplot_df)
  boxplot_df <- boxplot_df[c(rows[which_row],rows[-which_row]),]
  rownames(boxplot_df) <- NULL
  
  p1 <- ggplot(data=boxplot_df,
              aes(x=prep_and_tool, y=Bray_Curtis, fill=prep_and_tool)) +
    geom_boxplot() + 
    scale_fill_manual(values=c(
    "#1E75B0", "#ABC3E3", "#FA7D0E","#FAB776",
    "#2B9D2B", "#D90368", "#D22627","#FA9593",
    "#9165B9", "#FAF33E", "#89544A","#C09991",
    "#DF75BE", "#F2B2CE", "#7D7D7D","#C3C3C3",
    "#B8B921", "#D7D78A", "#17BACB","#9BD6E1",
    "#FFD131", "#F4AC32", "#DBB3B1","#A5BE00",
    "#3D52D5", "#FFEEDD","#1EFC1E","#C1ADD1",
    "#95DB87","#000000","#89544A")) + 
    scale_x_discrete(guide = guide_axis(angle = 90))  +
    theme_bw() + theme(legend.position = "none") + 
    ylab(label = "Bray-Curtis distance") +  
    xlab(label = "Pipeline") 
  
  boxplot_df_grouped <- boxplot_df %>%
    dplyr::group_by(prep_and_tool,taxonomy) %>%
    summarise(
      Bray_Curtis = mean(Bray_Curtis),
      taxonomy=unique(taxonomy)
    ) 
  
  boxplot_df_grouped$prep_and_tool <- factor(boxplot_df_grouped$prep_and_tool,
                                             levels=unique(boxplot_df$prep_and_tool))
  
  boxplot_df_grouped$taxonomy <-  factor(boxplot_df_grouped$taxonomy,
                                         levels=unique(boxplot_df$taxonomy))
  
  p2 <- ggplot(data=boxplot_df_grouped,
               aes(x=prep_and_tool, y=Bray_Curtis, color=taxonomy)) +
    geom_jitter(height = 0,width = 0.2,size=3) +
    scale_color_manual(values=c(
      "#1EFC1E","#1E75B0", "#FA7D0E", "#7D7D7D", "#D22627"
    )) + 
    theme_bw() +
    ylab(label = "Bray-Curtis") +  
    xlab(label = "") +
    scale_x_discrete(guide = guide_axis(angle = 90)) + 
    theme(
      axis.text.x = element_text(size = 8,color="black"),
      axis.text.y = element_text(size = 17,color="black"),
      axis.title.y= element_text(size = 17,color="black"))
  
  return(p2)
}

precision_recall <- function(asv_table, metadata){
  asv_table %<>% dplyr::select(-Color)
  
  zymo_taxa <- asv_table[,c("SeqID","Mock_reference")]
  zymo_taxa <- zymo_taxa[zymo_taxa$Mock_reference!=0,]
  
  precs <- c(1)
  recs <- c(1)
  samples <- c("Mock_reference")
  # calculating precision
  for (col in 3:ncol(asv_table)){
    sample <- colnames(asv_table)[col]
    taxa_sample <- asv_table[,c(1,col)]
    taxa_sample <- taxa_sample[taxa_sample[,2]!=0,]
    tp <- sum(taxa_sample$SeqID %in% zymo_taxa$SeqID)
    fp <- sum(!(taxa_sample$SeqID %in% zymo_taxa$SeqID))
    prec = tp/(tp + fp)
    precs <- c(precs,prec)

    tp <- sum(taxa_sample$SeqID %in% zymo_taxa$SeqID)
    fn <- sum(!(zymo_taxa$SeqID %in% taxa_sample$SeqID))
    rec = tp/(tp+fn)
    recs <- c(recs,rec)
    samples <- c(samples,sample)
  }
  
  new_df <- data.frame(SampleID=samples,
                       precision=precs,
                       recall=recs,
                       F1=(2*precs*recs)/(precs+recs))
  metadata <- merge(metadata,new_df,by="SampleID",all=TRUE)
  return(metadata)
}

performance_plot <- function(metadata){
  boxplot_df <- metadata %>%
    dplyr::mutate(prep_and_tool= paste(prep_number,process_tool,sep="_"))
  boxplot_df[boxplot_df$SampleID=="Mock_reference","prep_and_tool"] <- "Reference"
  
  boxplot_df$SettingID <- gsub("Mock_reference","Reference",boxplot_df$SettingID) 
  boxplot_df$SettingID <- factor(boxplot_df$SettingID,levels=(c("Reference",
                                                                         paste0("DA",1:48),
                                                                         paste0("VU",1:32),
                                                                         paste0("DE",1:32))))
  
  boxplot_df_grouped <- boxplot_df %>%
    dplyr::group_by(SettingID) %>%
    summarise(
      precision = mean(precision),
      recall = mean(recall),
      F1 = mean(F1)
    ) 
  
  boxplot_df_grouped$SettingID = factor(boxplot_df_grouped$SettingID,levels=(c("Reference",
                                                                               paste0("DA",1:48),
                                                                               paste0("VU",1:32),
                                                                               paste0("DE",1:32))))
  
  long_df <- melt(boxplot_df_grouped)
  wide_width <- 0.9
  thin_width <- 0.5
  
  recall <- "#3fbdc3"
  precision <- "#f07b70"
  # Create the plot
  p <- ggplot(long_df, aes(x = SettingID, y = value, fill = variable)) +
    geom_bar(data = subset(long_df, variable == "F1"),
             aes(y = value),
             width=wide_width,
             stat = "identity",
             fill = "gray") + 
    geom_bar(data = subset(long_df, variable == "recall"),
             width=thin_width,
             aes(y = value),
             stat = "identity",
             fill="#3fbdc3") +
    geom_bar(data = subset(long_df, variable == "precision"),
             aes(y = value),
             stat = "identity",
             width = thin_width,
             fill="#f07b70") +
    theme_bw() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    xlab("") + 
    ylab("Value") +  
    labs(fill = "Variable Type") + 
    theme(legend.position="bottom",
          axis.text.x = element_text(size = 8,color="black"),
          axis.text.y = element_text(size = 12,color="black"),
          axis.title.x = element_text(size = 12,color="black"),
          axis.title.y= element_text(size = 12,color="black"))
    
  
  return(p)
}
