#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  # Globals set up
  set_logger()
  BFC <<- get_validate_cache()

  SHEETS <<- from_bfc("sheets")
  CACHE <<- cachem::cache_mem()

  sheet_sidebar <- get_sheet_sidebar(SHEETS)
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    bslib::page_navbar(
      mod_stats_ui("stats_1"),
      bslib::nav_menu(
        "Sample catalog",
        mod_clinical_ui("clinical_1", sheet_sidebar),
        mod_metadata_ui("metadata_1", sheet_sidebar),
        mod_samples_ui("samples_1", sheet_sidebar),
      ),
      bslib::nav_menu(
        "RNA-seq",
        mod_compare_expression_ui("compare_expression_1", "Expression: Bulk"),
        mod_compare_expression_ui(
          "compare_expression_2",
          "Expression: Single-cell (pseudobulk)"
        ),
        mod_de_ui("de_1", "DE Analysis: Bulk")
      ),
      theme = bslib::bs_theme(bootswatch = "flatly"),
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "StemDataExplorer"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
