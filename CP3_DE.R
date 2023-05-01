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
source("helper.R")


ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Final Project")),
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

server <- function(input, output, session) {
  load_deseq2 <- reactive({
    input$deseq2_results
    isolate({
      # Change this when submitting
      deseq2_results <- read.csv(input$deseq2_results$datapath)
      #deseq2_results <- read.csv("data/GSE64810_DESeq2_results.csv")
      return(deseq2_results)
    })
  })
  
  observeEvent(input$render_deseq2, {
    output$deseq2_plot <- renderPlotly({
      input$render_deseq2
      isolate({
        plot_deseq_res(load_deseq2(), input$DE_x, input$DE_y, input$DE_slider_padj, input$DE_base, input$DE_highlight)
      })
    })
  })
  
  observeEvent(input$render_deseq2, {
    output$deseq2_table <- DT::renderDataTable(
      DT::datatable(isolate(draw_table_deseq2(load_deseq2(), input$DE_slider_padj)),
                    class = "display",
                    options = list(paging = TRUE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
