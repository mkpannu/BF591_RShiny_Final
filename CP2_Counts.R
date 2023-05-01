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
source('helper.R')


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
  
  observeEvent(input$render_counts, {
    output$filtersummary <- DT::renderDataTable(
      DT::datatable(filter_summary(isolate(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  observeEvent(input$render_counts, {
    output$var_scatter <- renderPlotly({
      isolate({
        plot_var(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
      })
    })
  })
  
  observeEvent(input$render_counts, {
    output$nonzero_scatter <- renderPlotly({
      isolate({
        plot_nonzero(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
      })
    })
  })
  
  observeEvent(input$render_counts, {
    output$heatmap <- renderPlot({
      isolate({
        showModal(modalDialog("Producing Heatmap", footer=NULL))
        plot_heatmap(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
        removeModal()
      })
    })
  })
  
  observeEvent(input$plot_pca, {
    output$pca_plot <- renderPlotly({
      input$plot_pca
      isolate({
        plot_pca(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter), input$pca_1, input$pca_2)
      })
    })
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
