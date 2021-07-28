################################################################################
# Shiny dashboard - User Interface

# Global options
options(
  shiny.maxRequestSize = 30*10^9,
  spinner.color        = "#0dc5c1"
)

# Load heler functions
source("utils.R")
source("scatterplot.R")

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
    selected = "Stimulated monocytes"
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
          # Metadata on embedding
          shiny::column(
            width = 6, offset = 0
            ,
            # Plot
            shiny::plotOutput("bal_embedding_metadata")
            ,
            # Left side
            shiny::column(
              width = 6, offset = 0
              ,
              shiny::uiOutput("bal_embedding_metadata_coldata")
            )
            ,
            # Right side
            shiny::column(
              width = 6, offset = 0
            )
          )
          ,
          # Gene expression on embedding
          shiny::column(
            width = 6, offset = 0
            ,
            # Plot
            shiny::plotOutput("bal_embedding_expression")
            ,
            # Left side
            shiny::column(
              width = 6, offset = 0
              ,
              shiny::selectizeInput(
                inputId  = "bal_embedding_expression_coldata",
                label    = "Select gene",
                choices  = NULL,
                multiple = FALSE
              )
            )
            ,
            # Right side
            shiny::column(
              width = 6, offset = 0
            )
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
      ,
      datasetInput("balmac")
      ,
      # Navigation list
      shiny::navlistPanel(
        well = FALSE, widths = c(2, 10)
        ,
        "Navigation"
        ,
        shiny::tabPanel(
          title = "Overview",
          scatterUI("balmac_meta"),
          scatterUI("balmac_expr")
        )
      )
    ) # tabPanel
    ,
    # --------------------------------------------------------------------------
    # Stimulated monocytes
    shiny::tabPanel(
      title = "Stimulated monocytes"
      ,
      datasetInput("Monocytes")
      ,
      # Navigation list
      shiny::navlistPanel(
        well = FALSE, widths = c(2, 10)
        ,
        "Navigation"
        ,
        shiny::tabPanel(
          title = "Overview",
          scatterUI("mono_meta"),
          scatterUI("mono_expr")
        )
      )
    ) # tabPanel
    
  )
  # navbarPage
)
# ui


# end of document
################################################################################