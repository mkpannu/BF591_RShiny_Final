## Author: Mahek Pannu
## mkpannu@bu.edu
## BU BF591
## Final Project

# Load necessary libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(tidyverse)
library(plotly)
library(colourpicker) 
library(data.table)
library(RColorBrewer)
source("helper.R")

# Declare UI for the app
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Final Project")),
  titlePanel(h3("Mahek Pannu")),
  tabsetPanel(type = "tabs", 
              ## Sample Exploration Tab
              tabPanel("Sample", sidebarLayout(
                sidebarPanel(
                  HTML("<h2>Sample Statistics Input</h2>"),
                  fileInput("sample_stats", paste0("Load sample statistics as a csv"),  placeholder = "sample_statistics.csv", accept='.csv'),
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
              )), 
              ## Counts Exploration Tab
              tabPanel("Counts", sidebarLayout(
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
                                           sliderInput("pca_2", "Select Second Principal Component:", value=2, min=1, max=69),
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
              )),
              ## Differential Expression Tab
              tabPanel("Differential Expression", sidebarLayout(
                sidebarPanel(
                  fileInput("deseq2_results", paste0("Load differential expression results"),  placeholder = "GSE64810_DESeq2_results.csv", accept='.csv'),
                  HTML("<p>A volcano plot can be generated with <b>\"log2 fold-change\"</b> on the x-axis and <b>\"p-adjusted\"</b> on the y-axis.</p>"),
                  radioButtons("DE_x", label = "Choose the column for the x-axis", choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected="log2FoldChange"),
                  radioButtons("DE_y", label = "Choose the column for the y-axis", choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"), selected="padj"),
                  colourInput("DE_base", "Base point color", "#22577A"),
                  colourInput("DE_highlight", "Highlight point color", "#FFCF56"), 
                  sliderInput("DE_slider_padj", "Select the magnitude of the p adjusted coloring:", value =-10, min=-35, max=0),
                  actionButton("render_deseq2", "Submit", class = "btn-success"), 
                  width=3
                ),
                mainPanel(
                  tabsetPanel(type = "tabs",
                              tabPanel("DEseq2 Plot", plotlyOutput("deseq2_plot")),
                              tabPanel("DEseq2 Results", dataTableOutput("deseq2_table"))
                  ),
                  width=9
                )
              )), 
              ## FGSEA Exploration Tab
              tabPanel("FGSEA", sidebarLayout(
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
              ))
  )
)

server <- function(input, output, session) {
  ## Declares max file size (specifically needed for counts csv)
  options(shiny.maxRequestSize=30*1024^2)
  
  ## Sample Exploration Features 
  ### Load sample statistics data as a reactive
  ### Return sample statistics as a dataframe
  load_data_sample <- reactive({
    input$render_sample
    isolate({
      sample_stats <- read.csv(input$sample_stats$datapath)
      sample_stats$Condition <- as.factor(sample_stats$Condition)
      return(sample_stats)
    })
  })
  
  ### Renders summary data table when submit button pressed
  ### Takes reactive sample dataframe 
  ### Uses helper function create_summary()
  observeEvent(input$render_sample, {
    output$summarystats <- DT::renderDataTable(
      DT::datatable(create_summary(isolate(load_data_sample())), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  ### Renders full sample data table when submit button pressed
  ### Takes reactive sample dataframe 
  observeEvent(input$render_sample, {
    output$sample_table <- DT::renderDataTable(
      DT::datatable(isolate(load_data_sample()), 
                    class = "display",
                    options = list(paging = FALSE, 
                                   fixedColumns = TRUE, 
                                   ordering = TRUE, 
                                   dom = 'Brtip'
                    )
      )
    )
  })
  
  ### Renders Violin plot of continous variables when submit button for specific plot is pressed
  ### Takes reactive sample dataframe and variable to plot control and HD samples from user specified radio button
  ### Uses helper function HD_control()
  observeEvent(input$sample_HDC_plot, {
    output$HD_control_plot <- renderPlotly({
      input$sample_HDC_plot
      isolate({
        HD_control(load_data_sample(), input$samp_var)
      })
    })
  })
  
  ### Renders Density plot of continous variables for HD samples only when submit button for specific plot is pressed
  ### Takes reactive sample dataframe and variable to plot HD samples from user specified radio button
  ### Uses helper function HD_plot()
  observeEvent(input$sample_HD_plot, {
    output$HD_plot <- renderPlotly({
      input$sample_HD_plot
      isolate({
        HD_plot(load_data_sample(), input$HD_var)
      })
    })
  })
  
  ## Counts Exploration Tab Features
  ### Load normalized count data as a reactive
  ### Return normalized count data as a dataframe
  load_data_counts <- eventReactive(input$render_counts, {
    isolate({
      norm_counts <- read.csv(input$counts_input$datapath)
      return(norm_counts)
    })
  })
  
  ### Renders summary of filter results when submit button pressed
  ### Takes reactive normalized count dataframe and inputs of variance and nonzero from user specified radio input
  ### Uses helper function filter_data() and filter_summary()
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
  
  ### Renders variance scatter plot of continuous variables for normalized counts only when general submit button is pressed
  ### Variance Plot: x-axis: log10(variance) and y-axis: median counts
  ### Takes reactive normalized count dataframe and variance and nonzero filter from user specified radio button
  ### Uses helper function filter_data() and plot_var()
  observeEvent(input$render_counts, {
    output$var_scatter <- renderPlotly({
      isolate({
        plot_var(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
      })
    })
  })
  
  ### Renders number of nonzero scatter plot of continuous variables for normalized counts only when general submit button is pressed
  ### Nonzero Plot: x-axis: Number of Zero counts and y-axis: median counts
  ### Takes reactive normalized count dataframe and variance and nonzero filter from user specified radio button
  ### Uses helper function filter_data() and plot_nonzero()
  observeEvent(input$render_counts, {
    output$nonzero_scatter <- renderPlotly({
      isolate({
        plot_nonzero(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
      })
    })
  })
  
  ### Renders heatmap of filtered genes when general submit button is pressed
  ### Takes reactive normalized count dataframe and variance and nonzero filter from user specified radio button
  ### Uses helper function filter_data() and plot_heatmap()
  observeEvent(input$render_counts, {
    output$heatmap <- renderPlot({
      isolate({
        showModal(modalDialog("Producing Heatmap", footer=NULL))
        plot_heatmap(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter))
        removeModal()
      })
    })
  })
  
  ### Renders pca plot of filtered genes when plot specific submit button is pressed
  ### Takes reactive normalized count dataframe, variance and nonzero filter, two pca vectors from user specified radio button
  ### Uses helper function filter_data() and plot_pca()
  observeEvent(input$plot_pca, {
    output$pca_plot <- renderPlotly({
      input$plot_pca
      isolate({
        plot_pca(filter_data(load_data_counts(), input$var_filter, input$nonzero_filter), input$pca_1, input$pca_2)
      })
    })
  })
  
  ## Differential Exploration Tab Features
  ### Load deseq2 results data as a reactive
  ### Return deseq2 results as a dataframe
  load_deseq2 <- reactive({
    input$deseq2_results
    isolate({
      deseq2_results <- read.csv(input$deseq2_results$datapath)
      return(deseq2_results)
    })
  })
  
  ### Renders customizable plot when general submit button is pressed
  ### Takes reactive deseq2 results dataframe and user specified x-axis and y-axis columns to plot, and two user specified colors
  ### Highlight color are genes that pass the filter vs. base color are genes that fail the filter
  ### Uses helper function plot_deseq_res()
  observeEvent(input$render_deseq2, {
    output$deseq2_plot <- renderPlotly({
      input$render_deseq2
      isolate({
        plot_deseq_res(load_deseq2(), input$DE_x, input$DE_y, input$DE_slider_padj, input$DE_base, input$DE_highlight)
      })
    })
  })
  
  ### Renders deseq2 results when submit button pressed filtered by user specified padj 
  ### Takes reactive deseq2 results dataframe and inputs of specified padj from user specified radio input
  ### Uses helper function draw_table_deseq2()
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
  
  ## FGSEA Exploration Tab Features
  ### Load FGSEA results data as a reactive
  ### Return FGSEA results data as a dataframe
  load_data_fgsea <- reactive({
    input$render_fgsea
    isolate({
      # Change this when submitting
      fgsea_results <- read.csv(input$fgsea_csv$datapath)
      #fgsea_results <- read.csv("data/fgsea_results.csv")
      return(fgsea_results)
    })
  })
  
  ### Renders Top pathways based on NES score when plot specific submit button is pressed
  ### Takes reactive FGSEA results dataframe and user specified number of pathways from slider input
  ### Uses helper function pathway_plot()
  observeEvent(input$render_paths, {
    output$plot_pathways <- renderPlot({
      input$render_paths
      isolate({
        pathway_plot(load_data_fgsea(), input$top_pathways)
      })
    })
  })
  
  ### Renders Top pathways based on NES score when plot specific submit button is pressed
  ### Takes reactive FGSEA results dataframe and user specified number of pathways from slider input
  ### Uses helper function pathway_plot()
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
  
  ### Renders FGSEA results when table specific submit button pressed
  ### Takes reactive FGSEA results dataframe and inputs of padj filter and NES pathways as specified by user
  ### Uses helper function fgsea_table() and fwrite() from data.table library
  observeEvent(input$render_fgsea_table, {
    output$download <- downloadHandler(
      filename = function(){"fgsea_filtered_results.csv"}, 
      content = function(fname){
        fwrite(fgsea_table(load_data_fgsea(), input$fgsea_padj_slider, input$NES_choice), file=fname, sep=",", sep2=c("", " ", ""))
      }
    )
  })
  
  ### Renders scatter plot when plot specific submit button is pressed
  ### X-axis: NES and y-axis: -log10(padj)
  ### Takes reactive FGSEA results dataframe and user specified padj filter from slider input
  ### Uses helper function scatter_fgsea()
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
