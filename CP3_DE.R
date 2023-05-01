## Author: Mahek Pannu
## mkpannu@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(tidyverse)
library(plotly)
library(colourpicker) # you might need to install this.


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Assignment 7")),
  sidebarLayout(
    sidebarPanel(
      fileInput("deseq2_results", paste0("Load differential expression results"),  placeholder = "GSE64810_DESeq2_results.csv", accept='.csv'),
      HTML("<p>A volcano plot can be generated with <b>\"log2 fold-change\"</b> on the x-axis and <b>\"p-adjusted\"</b> on the y-axis.</p>"),
      radioButtons("DE_x", label = "Choose the column for the x-axis", choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected="log2FoldChange"),
      radioButtons("DE_y", label = "Choose the column for the y-axis", choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected="padj"),
      colourInput("DE_base", "Base point color", "#22577A"),
      colourInput("DE_highlight", "Highlight point color", "#FFCF56"), 
      sliderInput("DE_slider_padj", "Select the magnitude of the p adjusted coloring:", value =-10, min=-35, max=0),
      actionButton("render_deseq2", "Plot", class = "btn-success")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("DEseq2 Plot", plotlyOutput("deseq2_plot")),
                  tabPanel("DEseq2 Results", dataTableOutput("deseq2_table"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  load_deseq2 <- reactive({
    input$deseq2_results
    isolate({
      # Change this when submitting
      #deseq2_results <- read.csv(input$sample_stats$datapath)
      deseq2_results <- read.csv("data/GSE64810_DESeq2_results.csv")
      return(deseq2_results)
    })
  })
 
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
  
  output$deseq2_plot <- renderPlotly({
    input$render_deseq2
    isolate({
      plot_deseq_res(load_deseq2(), input$DE_x, input$DE_y, input$DE_slider_padj, input$DE_base, input$DE_highlight)
    })
  })
  
  output$deseq2_table <- DT::renderDataTable(
    DT::datatable(draw_table_deseq2(load_deseq2(), input$DE_slider_padj),
                  class = "display",
                  options = list(paging = TRUE, 
                                 fixedColumns = TRUE, 
                                 ordering = TRUE, 
                                 dom = 'Brtip'
                  )
    )
  )
}

# Run the application
shinyApp(ui = ui, server = server)
