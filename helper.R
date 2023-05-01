library(ggplot2)
library(dplyr)
library(tidyverse)

## Helper functions for Sample Data (Sample Tab)
create_summary <- function(df){
  summary <- tibble(ColumnName=colnames(df), Type="", Values="")
  for(i in 1:length(colnames(df))){
    summary[i,2] <- class(df[[i]])
    if (summary[i,2] == "character"){
      summary[i,3] <- paste(df[[i]], collapse = ", ")
    }
    else if (summary[i,2] == "integer" || summary[i,2] == "numeric"){
      t <- c(mean(df[[i]], na.rm = TRUE), sd(df[[i]], na.rm = TRUE))
      summary[i,3] <- paste(t, collapse = " +/- ")
    }
    else if (summary[i,2] == "factor"){
      summary[i,3] <- paste(levels(df[[i]]), collapse = ", ")
    }
    else {
      summary[i,3] <- NULL
    }
  }
  return(summary)
}

HD_control <- function(data, col){
  p <- ggplot(data) +
    geom_violin(aes(x=Condition,y=!!sym(col),fill=Condition))
  return(p)
}

HD_plot <- function(data, col){
  p <- filter(data, Condition == "HD") %>%
    drop_na() %>%
    ggplot() +
    geom_density(aes(x=!!sym(col), fill=Condition))
  return(p)
}

filter_data <- function(data, var_filter, nonzero_filter){
  length_samples <- length(colnames(data))
  data$vars <- apply(data[-1], 1, var)
  data$zero_sum <- apply(data[-1] == 0, 1, sum) 
  data$filter <- as.logical(ifelse((data$vars > quantile(data$vars, var_filter)) == TRUE & (data$zero_sum < (length_samples - nonzero_filter)) == TRUE, TRUE, FALSE))
  return(data)
}

## Helper Function for Counts Data (Counts Tab)
filter_summary <- function(fltr_data){
  num_samples <- length(colnames(fltr_data)) - 2
  num_genes <- length(fltr_data$gene)
  num_pass <- length(dplyr::filter(fltr_data, filter ==TRUE)$filter)
  num_fail <- length(dplyr::filter(fltr_data, filter ==FALSE)$filter)
  fltr_summary <- tibble(Description=c("Number of Samples", "Number of Genes", "Number that Pass Filter", "Number that Fail Filter"),
                         Statistic= c(num_samples, num_genes, num_pass, num_fail), 
                         Percentage=c(NA, NA, paste0(round(num_pass/num_genes*100, 2), "%"), paste0(round(num_fail/num_genes*100, 2), "%")))
  return(fltr_summary)
}

plot_var <- function(fltr_data){
  fltr_data$med_counts <- apply(fltr_data[2:70], 1, median)
  p_var <- ggplot(fltr_data, mapping=aes(x=log10(vars), y=med_counts)) +
    geom_point(mapping=aes(color=filter)) + 
    scale_color_manual(values = c("#652CBA", "#FFCF56")) +
    labs(title="log(10) Variance of Genes that Pass Selected Filter",
         x="log10(variance)", 
         y="Median Counts") 
  return(p_var)
}

plot_nonzero <- function(fltr_data){
  fltr_data$med_counts <- apply(fltr_data[2:70], 1, median)
  p_zero <- ggplot(fltr_data, mapping=aes(x=zero_sum, y=med_counts)) +
    geom_point(mapping=aes(color=filter)) + 
    scale_color_manual(values = c("#652CBA", "#FFCF56")) +
    labs(title="Number of Zero Counts for Samples where Genes Pass Selected Filter",
         x="Number of Zero Counts", 
         y="Median Counts") 
  return(p_zero)
}

plot_heatmap <- function(fltr_data){
  ht_df <- filter(fltr_data, filter==TRUE)[2:70] %>%
    as.matrix()
  rownames(ht_df) <- filter(fltr_data, filter==TRUE)$gene
  ht_map <- heatmap(ht_df, col = brewer.pal(n = 11, name = "RdBu"))
  return(ht_map)
}

plot_pca <- function(norm_counts, pca_1, pca_2){
  fltr_counts <- filter(norm_counts, filter==TRUE)
  pca_vals <- prcomp(t(fltr_counts[2:70]))
  plot_pca_vals <- tibble(PC1 = pca_vals$x[,pca_1], PC2 = pca_vals$x [,pca_2], Condition = str_extract(colnames(fltr_counts[2:70]), "[A-Z]*"))
  exp_var <- pca_vals$sdev^2/sum(pca_vals$sdev^2)
  
  pca <- ggplot(plot_pca_vals, aes(x=PC1,y=PC2)) + 
    geom_point(aes(color=Condition)) +
    labs(title="DESeq2 Normalized PCA",
         x=paste0("PC", pca_1,": ", round(exp_var[pca_1]*100,2),"% variance"),
         y=paste0("PC", pca_2,": ", round(exp_var[pca_2]*100,2),"% variance"))
  return(pca)
}

## Helper Function for Differential Expression (Differential Expression Tab)
plot_deseq_res <- function(deseq2_results, x_name, y_name, slider, color1, color2) {
  p_df <- drop_na(deseq2_results) 
  p <- ggplot(p_df, mapping=aes(x=!!sym(x_name), y=-log10(as.numeric(!!sym((y_name)))))) +
    geom_point(size=1, aes(color=(!!sym(y_name) < (10^slider)))) + 
    scale_colour_manual(name = paste0(y_name, ' < 1e', slider), values = setNames(c(color1,color2),c(F, T))) 
  return(p)
}

draw_table_deseq2 <- function(deseq2_results, slider_padj) {
  fltr_df <- filter(deseq2_results, padj < 10^slider_padj) %>%
    mutate(across(baseMean:stat, ~ round(., 4)))
  fltr_df$padj <- format(fltr_df$padj)
  fltr_df$pvalue <- format(fltr_df$pvalue) 
  return(fltr_df)
}

## Helper functions for FGSEA (FGSEA Tab)
pathway_plot <- function(fgsea_results, n_paths){
  fgsea_results <- arrange(fgsea_results, desc(NES))
  top_pos <- fgsea_results %>% 
    top_n(n = n_paths, wt = NES) %>% 
    dplyr::mutate(sign = "pos")
  top_neg <- fgsea_results %>% 
    top_n(n = -n_paths, wt = NES) %>% 
    dplyr::mutate(sign = "neg")
  top_nes <- rbind(top_pos, top_neg) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway = gsub("_"," ", pathway), 
                  pathway = stringr::str_wrap(pathway, width = 80))
  p <- top_nes %>% 
    ggplot2::ggplot() + geom_col(aes(NES, pathway, fill=sign)) + 
    scale_y_discrete(limits = pull(top_nes, pathway)) +
    scale_fill_manual(values=c("red","blue")) +
    ggtitle("FGSEA results for Hallmark MSigDB gene sets") +
    xlab("Normalized Enrichment Score (NES)") + 
    theme(axis.text.y=element_text(size = 6), legend.position="bottom",
          axis.title.y=element_blank(), axis.title.x=element_text(size = 8),
          title=element_text(size = 8.5))
  return(p)
}

fgsea_table <- function(fgsea_results, padj_slider, NES_choice){
  if(NES_choice == "Positive"){
    fgsea_results <- filter(fgsea_results, NES > 0)
  }
  if(NES_choice == "Negative"){
    fgsea_results <- filter(fgsea_results, NES < 0)
  }
  fgsea_results <- filter(fgsea_results, padj < 10^(padj_slider)) %>%
    dplyr::mutate(pathway = gsub("_"," ", pathway)) %>%
    dplyr::mutate(leadingEdge = gsub(" ", ", ", leadingEdge)) %>%
    mutate(across(log2err:NES, ~ round(., 4)))
  
  fgsea_results$padj <- format(fgsea_results$padj)
  fgsea_results$pval <- format(fgsea_results$pval) 
  #fgsea_results$leadingEdge <- abbreviate(fgsea_results$leadingEdge)
  return(fgsea_results)
}

scatter_fgsea <- function(fgsea_results, padj_slider){
  p_df <- fgsea_results
  p_df$label <- fgsea_results$padj < 10^padj_slider
  p <- ggplot(p_df, mapping=aes(x=NES, y=-log10(padj))) + 
    geom_point(mapping=aes(color=padj<10^(padj_slider))) + 
    scale_color_manual(values = c("grey", "#652CBA")) + 
    labs(title="FGSEA Results that Pass Filter")
  return(p)
}