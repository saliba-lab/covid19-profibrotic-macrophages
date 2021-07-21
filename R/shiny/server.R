################################################################################
# Shiny dashboard - Server calculations

# Load heler functions
source("utils.R")

# Main function
server <- function(input, output, session) {
  
  # Close session when app closed
  shiny::onSessionEnded(
    session = session,
    fun     = stopApp
  )
  
  rv <- shiny::reactiveValues()
  
  # ----------------------------------------------------------------------------
  # BAL
  
  # Load data on click
  shiny::observeEvent(input$load_bal, {
    rv$bal <- readRDS("../../data/BAL.Rds")
  })
  
  # Show object dimensions
  output$bal_shape <- shiny::renderPrint({
    shiny::req(input$load_bal)
    glue::glue("{dim(rv$bal)[1]} genes across {dim(rv$bal)[2]} cells")
  })
  
  # Plot 
  output$bal_embedding_metadata <- shiny::renderPlot({
    shiny::req(rv$bal)
    plot.embedding(rv$bal)
  })
  
  # ----------------------------------------------------------------------------
  # BAL macrophages
  
  # ----------------------------------------------------------------------------
  # Monocytes
  
}
# server



# end of document
################################################################################