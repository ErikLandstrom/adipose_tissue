# Load libraries ----------------------------------------------------------

library(shiny)
library(tidyverse)
library(ggbeeswarm)


# Source helper functions -------------------------------------------------

source("impute_NAs_with_10_percent_of_mean.R")
source("plot_beeswarm.R")


# UI ----------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Differential expression analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file",
                label = "Choose a input file."),
      
      uiOutput("VarsInput")
      ,
                
      
      helpText("Missing values were imputed with 10% of the mean intensity for the purpose of vizualization.")),
    
    mainPanel(
      plotOutput("plot")
    )
  )
)
server <- function(input, output) {
  
  in_file <- output$file
  
  file <- impute_NAs_with_10_percent_of_mean(infile$datapath)
  
  output$VarsInput <- renderUI({
    selectizeInput("name",
                   label = "Select a protein to be plotted.",
                   choices = file$name,
                   options = list(maxOptions = 10))
  })
  
  output$plot <- renderPlot({
    gene_name <- input$name
    
    plot_beeswarm(file, gene_name)
  
  })
  
}


# Run app -----------------------------------------------------------------

shinyApp(ui, server)
