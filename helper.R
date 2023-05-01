library(ggplot2)
library(dplyr)
library(tidyverse)

########### Helper functions for Sample Data (Sample Tab) ###########
#' Create summary table for sample statistics tab 
#'
#' @param df Data frame loaded by load_data_sample()
#'
#' @return Tibble summarizing each column in sample statisitics data frame
#' with a column for ColumnName, column type, and unique/mean values in each column
#' 
#' @details This function initializes a tibble with each row being a column name 
#' in the sample statistics dataframe. If/else statements are used to evaluate
#' each column type to append to the returned tibble. The unique values if the 
#' column is a character or factor is filled in for value or for a numeric/integer
#' column the mean and sd are calculated and appended to the value column for that row
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

#' Violin Plot HD and control samples against user-specified column 
#'
#' @param data Data frame loaded by load_data_sample()
#' @param col User specified column to plot each sample type against
#'
#' @return ggplot of the violon plot of each sample plotted against user specified
#' column
#' 
#' @details This function takes a dataframe of the sample statistics and plots the 
#' user-specified column as the y-axis and Condition as the x-axis. Plot is filled 
#' by condition
HD_control <- function(data, col){
  p <- ggplot(data) +
    geom_violin(aes(x=Condition,y=!!sym(col),fill=Condition))
  return(p)
}

#' Density Plot HD samples against user-specified column 
#'
#' @param data Data frame loaded by load_data_sample()
#' @param col User specified column to plot HD sample type against
#'
#' @return ggplot of the denisty plot of HD sample plotted against user specified
#' column
#' 
#' @details This function takes a dataframe of the sample statistics and plots the 
#' user-specified column as the y-axis and HD sample as the x-axis. Data frame is 
#' filtered to only contain HD samples and NA values are dropped. Plot is filled
#' by condition
HD_plot <- function(data, col){
  p <- filter(data, Condition == "HD") %>%
    drop_na() %>%
    ggplot() +
    geom_density(aes(x=!!sym(col), fill=Condition))
  return(p)
}

########### Helper Function for Counts Data (Counts Tab) ###########
#' Dataframe filtered by variance and number of non-zeros as specified by user  
#'
#' @param data Data frame loaded by load_data_counts()
#' @param var_filter User specified quantile of variances to include 
#' @param nonzero_filter User specified number of nonzeros to include
#'
#' @return dataframe that adds columns with variance calculation and number of zeros 
#' and a boolean filter column that labels which columns pass both filters
#' 
#' @details This function takes a dataframe of counts and adds a column called var 
#' with the variance calculation for that row and a column called zero_sum with the
#' sum of zeros for that row. A third column is added called filter that is a boolean
#' that computes whether the row passes both filters as TRUE and if not as FALSE. 
filter_data <- function(data, var_filter, nonzero_filter){
  length_samples <- length(colnames(data))
  data$vars <- apply(data[-1], 1, var)
  data$zero_sum <- apply(data[-1] == 0, 1, sum) 
  data$filter <- as.logical(ifelse((data$vars > quantile(data$vars, var_filter)) == TRUE & (data$zero_sum < (length_samples - nonzero_filter)) == TRUE, TRUE, FALSE))
  return(data)
}

#' Summary Dataframe that describes how many genes passed or failed the filters  
#'
#' @param fltr_data Data frame from filter_data() of count data
#'
#' @return Tibble that summarizes the user-specified filters applied
#' 
#' @details This function takes a dataframe as returned by filter_data() and computes
#' the number of samples, number of genes, number that pass the filters and number that 
#' fails the filters. Additionally, the percentages for number of genes that pass/fail 
#' filters are computed. Those that pass/fail filters are found by using the filter
#' column to find. 
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

#' Diagnostic variance scatter plot 
#'
#' @param fltr_data Data frame from filter_data() of count data
#'
#' @return ggplot of the scatter plot with log10(variance) vs. Median Counts
#' 
#' @details This function takes a dataframe as returned by filter_data() and 
#' computes the median counts for each row. A scatter plot of the log10(variance)
#' on the x-axis and median counts on the y-axis is plotted with points passing
#' the filter being yellow and points not passing the filter being purple. 
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

#' Diagnostic nonzero scatter plot 
#'
#' @param fltr_data Data frame from filter_data() of count data
#'
#' @return ggplot of the scatter plot with Number of Zero counts vs. Median Counts
#' 
#' @details This function takes a dataframe as returned by filter_data() and 
#' computes the median counts for each row. A scatter plot of the Number of 
#' Zero counts on the x-axis and median counts on the y-axis is plotted with 
#' points passing the filter being yellow and points not passing the filter
#'  being purple. 
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

#' Heatmap of genes passing the filter
#'
#' @param fltr_data Data frame from filter_data() of count data
#'
#' @return Heatmap of genes passing filter
#' 
#' @details This function takes a dataframe as returned by filter_data() and 
#' uses only genes that pass both variance and non-zero filters. The filtered 
#' counts dataframe is passed as a matrix with genes as the rownames to the 
#' heatmap function. 
plot_heatmap <- function(fltr_data){
  ht_df <- filter(fltr_data, filter==TRUE)[2:70] %>%
    as.matrix()
  rownames(ht_df) <- filter(fltr_data, filter==TRUE)$gene
  ht_map <- heatmap(ht_df, col = brewer.pal(n = 11, name = "RdBu"))
  return(ht_map)
}

#' PCA plot of genes that pass the filter
#'
#' @param norm_counts Data frame from filter_data() of count data
#' @param pca_1 Principal Component Vector 1 as specified by user
#' @param pca_2 Principal Component Vector 2 as specified by user
#'
#' @return ggplot of the PCA plot of user-specified principal components
#' 
#' @details This function takes a dataframe as returned by filter_data() and 
#' uses only genes that pass both variance and non-zero filters. The transposed filtered 
#' counts dataframe is then passed into prcomp() to run PCA. The explained variance
#' is computed from the pca_vals for each principal component. ggplot is used to 
#' plot pca_1 and pca_2 as specified by the user. The explained variance is 
#' attached to each axis of the plot
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

########### Helper Function for Differential Expression (Differential Expression Tab) ########### 
#' User specified Plot of DeSeq2 Results
#'
#' @param deseq2_results Data frame from load_deseq2() of deseq2 results
#' @param x_name User specified x-axis
#' @param y_name User specified y-axis
#' @param slider User specified padj threshold
#' @param color1 User specified color for genes that fail filter
#' @param color2 User specified color for genes that pass filter
#'
#' @return ggplot of the deseq2_results with user specified columns
#' 
#' @details This function takes deseq2 results and drops any NA values. x_name 
#' is plotted as the x-axis and y_name is plotted as the y-axis with genes passing 
#' padj slider filter colored with color2 and genes failing padj slider filter
#' with color1 using ggplot
plot_deseq_res <- function(deseq2_results, x_name, y_name, slider, color1, color2) {
  p_df <- drop_na(deseq2_results) 
  p <- ggplot(p_df, mapping=aes(x=!!sym(x_name), y=-log10(as.numeric(!!sym((y_name)))))) +
    geom_point(size=1, aes(color=(!!sym(y_name) < (10^slider)))) + 
    scale_colour_manual(name = paste0(y_name, ' < 1e', slider), values = setNames(c(color1,color2),c(F, T))) 
  return(p)
}

#' Dataframe with the deseq2 results passing the filter 
#'
#' @param deseq2_results Data frame from load_deseq2() of deseq2 results
#' @param slider_padj User specified padj threshold filter
#'
#' @return Dataframe with genes passing padj threshold 
#' 
#' @details This function takes a dataframe of deseq2 results and filters
#' by user-specified padj threshold. The numeric values are formatted and
#' padj and pvalue are formatted by scientific notation. 
draw_table_deseq2 <- function(deseq2_results, slider_padj) {
  fltr_df <- filter(deseq2_results, padj < 10^slider_padj) %>%
    mutate(across(baseMean:stat, ~ round(., 4)))
  fltr_df$padj <- format(fltr_df$padj)
  fltr_df$pvalue <- format(fltr_df$pvalue) 
  return(fltr_df)
}

########### Helper functions for FGSEA (FGSEA Tab) ########### 
#' Plot of top pathways from FGSEA 
#'
#' @param fgsea_results Data frame from load_data_fgsea() of FGSEA results
#' @param n_paths User specified number of top pathways to include
#'
#' @return ggplot of n top pathways by NES value
#' 
#' @details This function takes FGSEA results and arranges the columns by NES. 
#' The top n positive and negative pathways are saved to use to make a ggplot of
#' top pathways. The pathway names are made human readable by replacing _ with 
#' spaces. ggplot is used to plot NES values by pathway with red bars being 
#' negative NES values and blue bars being positive NES values
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
                  pathway = stringr::str_wrap(pathway, width = 50))
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

#' Dataframe with FGSEA results passing filters 
#'
#' @param fgsea_results Data frame from load_data_fgsea() of FGSEA results
#' @param padj_slider User specified padj threshold filter
#' @param NES_choice positive, negative, or all NES values
#'
#' @return Dataframe with genes passing padj threshold and NES values
#' 
#' @details This function takes a dataframe of FGSEA results and filters
#' by user-specified padj threshold and user-specified NES value choice. The NES
#' values are filtered if user-specifies positive or negative values. The numeric
#' values are rounded and padj and pvalue are formatted by scientific notation. 
#' The leadingEdge column is reformatted to a list for easy download. 
fgsea_table <- function(fgsea_results, padj_slider, NES_choice){
  if(NES_choice == "Positive"){
    fgsea_results <- filter(fgsea_results, NES > 0)
  }
  if(NES_choice == "Negative"){
    fgsea_results <- filter(fgsea_results, NES < 0)
  }
  fgsea_results <- filter(fgsea_results, padj < 10^(padj_slider)) %>%
    dplyr::mutate(pathway = gsub("_"," ", pathway)) %>%
    mutate(across(log2err:NES, ~ round(., 4)))
  fgsea_results$leadingEdge <- apply(fgsea_results, 1, function(x) strsplit(x[8], " ")[[1]])
  fgsea_results$padj <- format(fgsea_results$padj)
  fgsea_results$pval <- format(fgsea_results$pval) 
  return(fgsea_results)
}

#' Scatter Plot of FGSEA results as filtered by padj threshold
#'
#' @param fgsea_results Data frame from load_data_fgsea() of FGSEA results
#' @param padj_slider User specified padj threshold filter
#'
#' @return scatter plot of genes that pass padj threshold 
#' 
#' @details This function takes FGSEA results and plots the NES value vs. -log10(padj)
#' value. Genes that do not pass the padj_slider threshold are colored gray and 
#' genes that do pass the padj_slider are colored purple. 
scatter_fgsea <- function(fgsea_results, padj_slider){
  p_df <- fgsea_results
  p <- ggplot(p_df, mapping=aes(x=NES, y=-log10(padj))) + 
    geom_point(mapping=aes(color=padj<10^padj_slider)) + 
    scale_color_manual(name = paste0('padj <1e^', padj_slider), values = c("grey", "#652CBA")) + 
    labs(title="FGSEA Results that Pass Filter")
  return(p)
}