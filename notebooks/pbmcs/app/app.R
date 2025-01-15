library(shiny)
library(reactable)
library(dplyr)

source('../functions.R')
data <- readRDS('../results/nested_table.rds') %>%
  dplyr::select(-path) %>%
  dplyr::rename(details = nested_table)

ui <- fluidPage(
  reactableOutput("table")
)

server <- function(input, output) {
  output$table <- renderReactable({
    reactable(
      dplyr::select(data, -c(celltype_sanitized, details)),
      details = make_get_details(data))
  })
}

shinyApp(ui, server)
