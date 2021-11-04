################################################################################
# Shiny dashboard - User Interface

# Global options
options(
  shiny.maxRequestSize = 30*10^9,
  spinner.color        = "#0dc5c1"
)

# Load helper functions
source("utils.R")
source("scatterplot.R")

# User interface function
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
      datasetInput("bal")
      ,
      shiny::navlistPanel(
        well = FALSE, widths = c(2, 10)
        ,
        "Navigation"
        ,
        shiny::tabPanel(
          title = "Overview",
          scatterUI("bal_meta"),
          scatterUI("bal_expr")
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
    )
    ,
    # --------------------------------------------------------------------------
    # Stimulated monocytes
    shiny::tabPanel(
      title = "Stimulated monocytes"
      ,
      datasetInput("Monocytes")
      ,
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
    )
    ,
    # --------------------------------------------------------------------------
    # App Info
    shiny::tabPanel(
      title = "App Info"
      ,
      shiny::h3("General"),
      shiny::p("This app contains three scRNA-seq datasets. To start browsing, 
               open any tab and select the 'Load data' button.")
      ,
      shiny::h3("Zooming"),
      shiny::p("Brush over points to zoom. Double click to return.")
    )
    
  )
  # navbarPage
)
# ui


# end of document
################################################################################