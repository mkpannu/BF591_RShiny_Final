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

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Final Project")),
  sidebarLayout(
    sidebarPanel(
      HTML("<h2>Sample Statistics Input</h2>"),
      fileInput("sample_stats", paste0("Load sample statistics"),  placeholder = "sample_statistics.csv", accept='.csv'),
      actionButton("render_sample", "Submit", class = "btn-success"), 
      width=3
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Summary", dataTableOutput("summarystats")),
                  tabPanel("Table", dataTableOutput("sample_table")),
                  tabPanel("Plot", 
                           sidebarLayout(
                             sidebarPanel(
                               HTML("<h4>Plotting continuous variables for both control and HD samples</h4>"), 
                               radioButtons("samp_var", label = "Choose column to plot", choiceNames = c("Post-Mortem Interval", "Age Of Death", "RNA Integrity Number"), choiceValues = c("PMI", "Age.Of.Death", "RIN"), selected="Age.Of.Death"),
                               actionButton("sample_HDC_plot", "Plot", class = "btn-success"), 
                               HTML("<h4>Plotting continuous variables HD samples only</h4>"), 
                               radioButtons("HD_var", label = "Choose column to plot", choiceNames = c("Post-Mortem Interval", "Age Of Death", "RNA Integrity Number", "Age of Onset", "Duration", "CAG in HTT Gene", "H.V Striatal Score", "H.V Corticol Score"), choiceValues = c("PMI", "Age.Of.Death", "RIN", "Age.of.Onset", "Duration", "CAG", "H.V.Striatal.Score", "H.V.Cortical.Score"), selected="Age.Of.Death"),
                               actionButton("sample_HD_plot", "Plot", class = "btn-success"), 
                               width=3
                             ),
                             mainPanel(
                               HTML("<h4>Violon Plot of continuous variables for both control and HD samples</h4>"), 
                               plotlyOutput("HD_control_plot"), 
                               HTML("<h4>Density Plot of continuous variables for HD samples only</h4>"), 
                               plotlyOutput("HD_plot"), 
                               width=9
                             )
                           )
                  )
      ),
      width=9
    )
  )
)

server <- function(input, output, session) {
  load_data_sample <- reactive({
    input$render_sample
    isolate({
      # Change this when submitting
      sample_stats <- read.csv(input$sample_stats$datapath)
      #sample_stats <- read.csv("data/sample_statistics.csv")
      sample_stats$Condition <- as.factor(sample_stats$Condition)
      return(sample_stats)
    })
  })
  
  observeEvent(input$render_sample, {
    output$summarystats <- DT::renderDataTable(
      DT::datatable(create_summary(load_data_sample()), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  observeEvent(input$render_sample, {
    output$sample_table <- DT::renderDataTable(
      DT::datatable(load_data_sample(), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  observeEvent(input$sample_HDC_plot, {
    output$HD_control_plot <- renderPlotly({
      input$sample_HDC_plot
      isolate({
        HD_control(load_data_sample(), input$samp_var)
      })
    })
  })
  
  observeEvent(input$sample_HD_plot, {
    output$HD_plot <- renderPlotly({
      input$sample_HD_plot
      isolate({
        HD_plot(load_data_sample(), input$HD_var)
      })
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
