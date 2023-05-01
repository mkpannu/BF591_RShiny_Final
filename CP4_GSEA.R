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
    input$render_sample
    isolate({
      # Change this when submitting
      #fgsea_results <- read.csv(input$fgsea_csv$datapath)
      fgsea_results <- read.csv("data/fgsea_results.csv")
      return(fgsea_results)
    })
  })
  
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
      dplyr::mutate(leadingEdge = gsub(" ", ", ", leadingEdge))
      
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
  
  output$plot_pathways <- renderPlot({
    input$render_paths
    isolate({
      pathway_plot(load_data_fgsea(), input$top_pathways)
    })
  })
  
  output$fgsea_table <- DT::renderDataTable(
    DT::datatable(fgsea_table(load_data_fgsea(), input$fgsea_padj_slider, input$NES_choice), 
                  class = "display",
                  options = list(paging = FALSE, 
                                 fixedColumns = TRUE, 
                                 ordering = TRUE, 
                                 dom = 'Brtip'
                  )
    )
  )
  
  output$download <- downloadHandler(
    filename = function(){"fgsea_filtered_results.csv"}, 
    content = function(fname){
      write.csv(fgsea_table(load_data_fgsea(), input$fgsea_padj_slider, input$NES_choice), fname)
    }
  )
  
  output$scatter_NES <- renderPlotly({
    input$render_fgsea_scatter
    isolate({
      scatter_fgsea(load_data_fgsea(), input$scatter_fgsea_slider)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
