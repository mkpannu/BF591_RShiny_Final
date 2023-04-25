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
library(colourpicker) 


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
  ),
  titlePanel(h1("BF591 Final Project")),
  tabsetPanel(type = "tabs", 
              tabPanel("Sample"), 
              tabPanel("Counts"),
              tabPanel("DE"), 
              tabPanel("Other")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
}

# Run the application
shinyApp(ui = ui, server = server)
