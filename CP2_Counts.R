## Author: Mahek Pannu
## mkpannu@bu.edu
## BU BF591
## Final Project

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(tidyverse)
library(plotly)
library(colourpicker) 
library('RColorBrewer')


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Final Project")),
  sidebarLayout(
    sidebarPanel(
      HTML("<h2>Counts Matrix Input</h2>"),
      fileInput("counts_input", paste0("Load normalized counts input"),  placeholder = "GSE64810_norm_counts.csv", accept='.csv'),
      sliderInput("var_filter", "Select to include genes with at least X percentile of variance:", value =0, min=0, max=1.00),
      sliderInput("nonzero_filter", "Select to include genes with at least X samples that are non-zero:", value=4, min=0, max=70),
      actionButton("render_counts", "Submit", class = "btn-success"), 
      width=3
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Filter Summary", dataTableOutput("filtersummary")),
                  tabPanel("Diagnostic Scatter Plots", plotlyOutput("var_scatter"), plotlyOutput("nonzero_scatter")),
                  tabPanel("Heat Map", plotOutput("heatmap")), 
                  tabPanel("PCA", 
                           sidebarLayout(
                             sidebarPanel(
                               sliderInput("pca_1", "Select First Principal Component:", value =1, min=1, max=69),
                               sliderInput("pca_2", "Select First Principal Component:", value=2, min=1, max=69),
                               actionButton("plot_pca", "Plot PCA", class = "btn-success"), 
                               width=3
                             ),
                             mainPanel(
                               plotlyOutput("pca_plot"))
                             )
                           )
                           
      ),
      width=9
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  load_data_counts <- eventReactive(input$render_counts, {
    isolate({
      # Change this when submitting
      norm_counts <- read.csv(input$counts_input$datapath)
      #norm_counts <- read.csv("data/GSE64810_norm_counts.csv")
      return(norm_counts)
    })
  })
  
  filter_data <- function(data, var_filter, nonzero_filter){
    length_samples <- length(colnames(data))
    data$vars <- apply(data[-1], 1, var)
    data$zero_sum <- apply(data[-1] == 0, 1, sum) 
    data$filter <- as.logical(ifelse((data$vars > quantile(data$vars, var_filter)) == TRUE & (data$zero_sum < (length_samples - nonzero_filter)) == TRUE, TRUE, FALSE))
    return(data)
  }
  
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
  
  output$filtersummary <- DT::renderDataTable(
      DT::datatable(filter_summary(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter)), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
  )
  
  
  output$var_scatter <- renderPlotly({
    input$render_counts
    isolate({
      plot_var(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
    })
  })
  
  output$nonzero_scatter <- renderPlotly({
    input$render_counts
    isolate({
      plot_nonzero(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
    })
  })
  
  output$heatmap <- renderPlot({
    input$render_counts
    isolate({
      plot_heatmap(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
    })
  })
  
  output$pca_plot <- renderPlotly({
    input$plot_pca
    isolate({
      plot_pca(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter), input$pca_1, input$pca_2)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
