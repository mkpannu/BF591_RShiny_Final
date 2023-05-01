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
      HTML("<h2>FGSEA Results Input</h2>"),
      fileInput("fgsea_csv", paste0("Load FGSEA results"),  placeholder = "fgsea_results.csv", accept='.csv'),
      actionButton("render_fgsea", "Submit", class = "btn-success"), 
      width=3
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Top Pathways", 
                           sidebarLayout(
                             sidebarPanel(
                               sliderInput("top_pathways", "Select top pathways to display:", value=10, min=0, max=25),
                               actionButton("render_paths", "Submit", class = "btn-success")
                             ), 
                             mainPanel(
                               plotOutput("plot_pathways")
                             ))),  
                  tabPanel("FGSEA Table", 
                           sidebarLayout(
                             sidebarPanel(
                               sliderInput("fgsea_padj_slider", "Filter by padj:", value=-10, min=-35, max=0),
                               radioButtons("NES_choice", label = "Select NES pathway:", choices = c("All", "Positive", "Negative"), choiceValues = c("PMI", ">", "<"), selected="All"),
                               actionButton("render_fgsea_table", "Submit", class = "btn-success"),
                               downloadButton('download',"Download the data")
                             ),
                             mainPanel(
                               dataTableOutput("fgsea_table")
                             ))),
                  tabPanel("Scatter Plot", 
                           sidebarLayout(
                             sidebarPanel(
                               sliderInput("scatter_fgsea_slider", "Filter by padj:", value=-4, min=-35, max=0),
                               actionButton("render_fgsea_scatter", "Submit", class = "btn-success")
                             ), 
                             mainPanel(
                               plotlyOutput("scatter_NES")
                             )
                           ))
      ),
      width=9
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  load_data_fgsea <- reactive({
    input$render_fgsea
    isolate({
      # Change this when submitting
      fgsea_results <- read.csv(input$fgsea_csv$datapath)
      #fgsea_results <- read.csv("data/fgsea_results.csv")
      return(fgsea_results)
    })
  })
  
  observeEvent(input$render_paths, {
    output$plot_pathways <- renderPlot({
      input$render_paths
      isolate({
        pathway_plot(load_data_fgsea(), input$top_pathways)
      })
    })
  })
  
  observeEvent(input$render_fgsea_table, {
    output$fgsea_table <- DT::renderDataTable(
      DT::datatable(isolate(fgsea_table(load_data_fgsea(), input$fgsea_padj_slider, input$NES_choice)), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  observeEvent(input$render_fgsea_table, {
    output$download <- downloadHandler(
      filename = function(){"fgsea_filtered_results.csv"}, 
      content = function(fname){
        write.csv(isolate(fgsea_table(load_data_fgsea(), input$fgsea_padj_slider, input$NES_choice), fname))
      }
    )
  })
  
  observeEvent(input$render_fgsea_scatter, {
    output$scatter_NES <- renderPlotly({
      isolate({
        scatter_fgsea(load_data_fgsea(), input$scatter_fgsea_slider)
      })
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
