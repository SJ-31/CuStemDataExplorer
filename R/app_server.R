#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  mod_clinical_server("clinical_1", SHEETS)
  mod_samples_server("samples_1", SHEETS)
  mod_metadata_server("metadata_1", SHEETS)
  mod_compare_expression_server("compare_expression_1")
}
