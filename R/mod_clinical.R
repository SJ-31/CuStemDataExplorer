#' clinical UI Function
#'
#' @description Tab for clinical data
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_clinical_ui <- function(id, sb) {
  ns <- NS(id)
  bslib::nav_panel(
    "Clinical",
    bslib::layout_sidebar(
      reactable::reactableOutput(ns("clin_tab")),
      sidebar = sb
    )
  )
}

#' clinical Server Functions
#'
#' @noRd
mod_clinical_server <- function(id, sheet_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$clin_tab <- reactable::renderReactable(reactable::reactable(
      sheet_data$clinical,
      searchable = TRUE,
      wrap = FALSE,
      columns = list(
        Note = reactable::colDef(
          cell = function(value, index, name) "",
          html = TRUE,
          details = reactable::JS(
            "function(rowInfo) {
        return `<br>${rowInfo.values['Note']}<br>`
}"
          ),
        ),
        `Has data?` = bool_col
      )
    ))
  })
}

## To be copied in the UI
# mod_clinical_ui("clinical_1")

## To be copied in the server
# mod_clinical_server("clinical_1")
