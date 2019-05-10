# Load libraries ----------------------------------------------------------

library(shiny)
library(tidyverse)
library(ggbeeswarm)

# Load data ---------------------------------------------------------------

data <- read_tsv("MIDY_adipose_no742_filtered_prots_190121_long.txt")

# Source helper functions -------------------------------------------------

source("impute_NAs_with_10_percent_of_mean.R")
source("plot_beeswarm.R")


# Impute ------------------------------------------------------------------

data <- impute_NAs_with_10_percent_of_mean(data)


# UI ----------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Differential expression analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("data2",
                label = "Choose a input file."),
      
      selectizeInput("name",
                label = "Select a protein to be plotted.",
                choices = data$name,
                options = list(maxOptions = 10)),
                
      
      helpText("Missing values were imputed with 10% of the mean intensity for the purpose of vizualization.")),
    
    mainPanel(
      plotOutput("plot")
    )
  )
)
server <- function(input, output) {
  
  output$plot <- renderPlot({
    gene_name <- input$name
    
    plot_beeswarm(data, gene_name)
  
  })
  
}


# Run app -----------------------------------------------------------------

shinyApp(ui, server)
