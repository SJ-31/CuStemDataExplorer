#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  sheet_data <- get_sheets()
  mod_clinical_server("clinical_1", sheet_data)
  mod_samples_server("samples_1", sheet_data)
  mod_metadata_server("metadata_1", sheet_data)
}
