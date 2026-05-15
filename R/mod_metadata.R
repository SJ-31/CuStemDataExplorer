#' metadata UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_metadata_ui <- function(id) {
  ns <- NS(id)
  bslib::nav_panel(
    "Metadata",
    reactable::reactableOutput(ns("meta_tab"))
  )
}

#' metadata Server Functions
#'
#' @noRd
mod_metadata_server <- function(id, sheet_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    meta_tab <- reactable::reactable(
      sheet_data$meta,
      columnGroups = list(
        reactable::colGroup(
          name = "Dates",
          columns = c("Sample collection", "Received data")
        )
      ),
      columns = list(
        `Sample collection` = reactable::colDef(
          format = reactable::colFormat(date = TRUE)
        ),
        `Received data` = reactable::colDef(
          format = reactable::colFormat(date = TRUE)
        )
      )
    )
    output$meta_tab <- reactable::renderReactable(meta_tab)
  })
}

## To be copied in the UI
# mod_metadata_ui("metadata_1")

## To be copied in the server
# mod_metadata_server("metadata_1")
