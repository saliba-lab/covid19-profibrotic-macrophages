################################################################################
# Shiny dashboard - Server calculations

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
  
  # Data upload ================================================================
  # Load data on click
  shiny::observeEvent(input$load_bal, {
    rv$bal <- readRDS("../../data/BAL.Rds")
  })
  # Show object dimensions
  output$bal_shape <- shiny::renderPrint({
    shiny::req(input$load_bal)
    glue::glue("{dim(rv$bal)[1]} genes across {dim(rv$bal)[2]} cells")
  })
  
  # Embedding metadata =========================================================
  # Plot 
  output$bal_embedding_metadata <- shiny::renderPlot({
    shiny::req(rv$bal)
    plot.embedding(
      rv$bal, coldata = input$bal_embedding_metadata_coldata,
      pt.size = 0.1
      )
  })
  # Select coldata
  output$bal_embedding_metadata_coldata <- shiny::renderUI({
    shiny::req(rv$bal)
    shiny::selectInput(
      inputId  = "bal_embedding_metadata_coldata",
      label    = "Select cell-specific metadata",
      choices  = names(rv$bal@meta.data), 
      selected = tail(names(rv$bal@meta.data), 1),
      multiple = FALSE
    )
  })
  
  # Embedding expression =======================================================
  # Plot 
  output$bal_embedding_expression <- shiny::renderPlot({
    shiny::req(rv$bal)
    plot.embedding(
      rv$bal, coldata = input$bal_embedding_expression_coldata, 
      pt.size = 0.1
      )
  })
  # Select coldata
  shiny::observe({
    shiny::req(rv$bal)
    choices <- rownames(rv$bal@assays$RNA@data)
    choices <- convertFeatures(choices, rv$bal@misc$features, 1, 2)
    shiny::updateSelectizeInput(
      session = session,
      inputId = "bal_embedding_expression_coldata",
      choices = choices,
      server  = TRUE
    )
  })
  
  # ----------------------------------------------------------------------------
  # BAL macrophages
  
  # Load monocyte data on click
  datasetServer(id = "balmac", datapath = "../../data/BAL.Rds", rv)
  
  # Create scatterplot
  scatterServer(id = "balmac_meta", key = "balmac", rv, "Metadata")
  scatterServer(id = "balmac_expr", key = "balmac", rv, "Gene Expression")
  
  # ----------------------------------------------------------------------------
  # Monocytes
  
  # Load monocyte data on click
  datasetServer(id = "Monocytes", datapath = "../../data/Monocytes.Rds", rv)
  
  # Create scatterplot
  scatterServer(id = "mono_meta", key = "Monocytes", rv, "Metadata")
  scatterServer(id = "mono_expr", key = "Monocytes", rv, "Gene Expression")
}
# server



# end of document
################################################################################