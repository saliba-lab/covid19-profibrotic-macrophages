################################################################################
# Shiny dashboard - User Interface

# Global options
options(
  shiny.maxRequestSize = 30*10^9,
  spinner.color        = "#0dc5c1"
)

# Main function
ui <- shiny::fluidPage(
  theme = shinythemes::shinytheme(theme = "flatly")
  ,
  # Main page organization
  shiny::navbarPage(
    title = paste(
      "SARS-CoV-2 infection triggers profibrotic",
      "macrophage responses and lung fibrosis"
    )
    ,
    # --------------------------------------------------------------------------
    # BAL
    shiny::tabPanel(
      title = "BAL"
      ,
      # Top row
      shiny::fluidRow(
        # Column
        shiny::column(
          width = 1, offset = 0.5
          ,
          shiny::actionButton(
            inputId = "load_bal",
            label   = "Load BAL data", 
            icon = shiny::icon("database")
          )
        )
        ,
        # Column
        shiny::column(
          width = 11, offset = 0.5
          ,
          shinycssloaders::withSpinner(
            shiny::verbatimTextOutput("bal_shape"), 
            type = 5, size = 0.5, proxy.height = "50px"
          )
        )
      )
      ,
      # Navigation list
      shiny::navlistPanel(
        well = FALSE, widths = c(2, 10)
        ,
        "Navigation"
        ,
        shiny::tabPanel(
          title = "Overview"
          ,
          shiny::column(
            width = 6, offset = 0
            ,
            shiny::plotOutput("bal_embedding_metadata")
          )
          ,
          shiny::column(
            width = 6, offset = 0
          )
        )
        ,
        shiny::tabPanel(
          title = "Differential expression"
        )
      )
    )
    ,
    # --------------------------------------------------------------------------
    # BAL macrophages
    shiny::tabPanel(
      title = "BAL macrophages"
    )
    ,
    # --------------------------------------------------------------------------
    # Stimulated monocytes
    shiny::tabPanel(
      title = "Stimulated monocytes"
    )
  )
  # navbarPage
)
# ui


# end of document
################################################################################