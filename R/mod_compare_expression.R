#' compare_expression UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_compare_expression_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' compare_expression Server Functions
#'
#' @noRd
mod_compare_expression_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_compare_expression_ui("compare_expression_1")

## To be copied in the server
# mod_compare_expression_server("compare_expression_1")
