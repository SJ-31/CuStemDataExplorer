#' samples UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_samples_ui <- function(id, sb) {
  ns <- NS(id)
  bslib::nav_panel(
    "Available modalities",
    bslib::layout_sidebar(
      reactable::reactableOutput(ns("sample_tab")),
      sidebar = sb
    )
  )
}

#' samples Server Functions
#'
#' @noRd
mod_samples_server <- function(id, sheet_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data <- sheet_data$all$data()
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
        Path = reactable::colDef(
          filterable = FALSE,
          html = TRUE,
          style = "font-family: monospace;",
          cell = \(v, i, n) {
            val <- data$Path[[i]]
            stringr::str_trunc(val, width = 10, side = "right")
          },
          details = \(i) {
            glue::glue("<pre>{data$Path[i]}</pre>")
          },
        ),
        `Processed files` = reactable::colDef(
          cell = function(value, index, name) {
            val <- data$`Processed files`[index]
            if (!is.na(val)) {
              "📁"
            } else {
              ""
            }
          },
          html = TRUE,
          details = \(index) {
            val <- data$`Processed files`[index]
            if (!is.na(val)) {
              paste0(
                "<pre>",
                stringr::str_replace_all(val, ";", "</br>"),
                "</pre>"
              )
            }
          }
        ),
        Warnings = reactable::colDef(
          cell = function(value, index, name) {
            if (!is.na(data$`Warnings`[index])) {
              "⁉"
            } else {
              ""
            }
          },
          html = TRUE,
          details = \(index) {
            val <- data$`Warnings`[index]
            if (!is.na(val)) {
              stringr::str_replace_all(val, ";", "</br>")
            }
          }
        )
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
