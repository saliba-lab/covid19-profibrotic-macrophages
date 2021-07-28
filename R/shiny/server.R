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
  
  # Load data on click
  datasetServer(id = "bal", datapath = "../../data/BAL.Rds", rv)
  
  # Create scatterplot
  scatterServer(id = "bal_meta", key = "bal", rv, "Metadata")
  scatterServer(id = "bal_expr", key = "bal", rv, "Gene Expression")
  
  # ----------------------------------------------------------------------------
  # BAL macrophages
  
  # Load data on click
  datasetServer(id = "balmac", datapath = "../../data/BAL-macrophages.Rds", rv)
  
  # Create scatterplot
  scatterServer(id = "balmac_meta", key = "balmac", rv, "Metadata")
  scatterServer(id = "balmac_expr", key = "balmac", rv, "Gene Expression")
  
  # ----------------------------------------------------------------------------
  # Monocytes
  
  # Load data on click
  datasetServer(id = "Monocytes", datapath = "../../data/Monocytes.Rds", rv)
  
  # Create scatterplot
  scatterServer(id = "mono_meta", key = "Monocytes", rv, "Metadata")
  scatterServer(id = "mono_expr", key = "Monocytes", rv, "Gene Expression")
}
# server



# end of document
################################################################################
