#' samples UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_samples_ui <- function(id) {
  ns <- NS(id)
  tagList(
    bslib::nav_panel(
      "Catalog",
      reactable::reactableOutput("mtab")
    )
  )
}

#' samples Server Functions
#'
#' @noRd
mod_samples_server <- function(id, sheet_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$sample_tab <- reactable::renderReactable(reactable::reactable(
      sheet_data$all,
      pagination = FALSE,
      searchable = TRUE,
      showPagination = TRUE,
      resizable = TRUE,
      wrap = FALSE,
      columns = list(
        PBMC = bool_col,
        Raw = bool_col,
        Processed = bool_col,
        Tumor = bool_col,
        Path = reactable::colDef(filterable = FALSE)
      ),
      elementId = "main_tab",
      columnGroups = list(
        reactable::colGroup(
          name = "Available Data",
          columns = c("Processed", "Raw")
        ),
        reactable::colGroup(
          name = "Available Samples",
          columns = c("Tumor", "PBMC")
        )
      )
    ))
  })
}

## To be copied in the UI
# mod_samples_ui("samples_1")

## To be copied in the server
# mod_samples_server("samples_1")
