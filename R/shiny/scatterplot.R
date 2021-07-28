################################################################################
# Shiny dashboard - Scatterplot Module

#' 2D scatterplot - UI
#' 
#' @id Shiny module key
#' @datapath Path to data
#' @rv Reactive values to store data
scatterUI <- function(id) {
  shiny::column(
    width = 6, offset = 0
    ,
    # Plot
    shiny::plotOutput(
      outputId = shiny::NS(id, "plot"), 
      brush    = shiny::brushOpts(id = "plot.brush"),
      dblclick = "plot.dblclick"
      )
    ,
    # Left side
    shiny::column(
      width = 6, offset = 0
      ,
      shiny::selectInput(
        inputId = shiny::NS(id, "coldata"), 
        label   = "Select cell-specific metadata", choices = NULL
        )
    )
    ,
    # Right side
    shiny::column(
      width = 6, offset = 0
      ,
      shiny::sliderInput(
        inputId = shiny::NS(id, "pointsize"), label = "Select point size",
        min = 0.1, max = 10, step = 0.1, value = 0.5
      )
    )
  )
}

#' 2D scatterplot - Server 
#' 
#' @id Shiny module key
#' @key Data storage key (identical to datasetInput() id)
#' @rv Reactive values list containing data (access by key)
#' @type 'Metadata' or 'Gene expression'
scatterServer <- function(id, key, rv, type) {
  shiny::moduleServer(id, function(input, output, session) {
    # Plot
    output$plot <- shiny::renderPlot({
      shiny::req(rv[[key]])
      plot.embedding(
        rv[[key]], coldata = input$coldata,
        pt.size = input$pointsize
      )
    })
    # Select coldata
    shiny::observeEvent(rv[[key]], {
      # Decide between metadata and expression
      if (type == "Metadata") {
        choices <- names(rv[[key]]@meta.data)
      } else {
        choices <- convertFeatures(
          rownames(rv[[key]]), convtab = ds@misc$features, 1, 2
          )
      }
      shiny::updateSelectizeInput(
        session, "coldata", 
        choices  = choices
        )
    })
  })
}

# end of document
################################################################################