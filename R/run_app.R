#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  # Globals set up
  set_logger()
  BFC <<- get_validate_cache()
  SHEETS <<- from_bfc("sheets")
  CACHE <<- cachem::cache_mem()

  if (Sys.getenv("SHINYLOADTEST") == "") {
    logger::log_info("Loading authentication module")
    app_ui <- shinymanager::secure_app(
      app_ui,
      theme = bslib::bs_theme(bootswatch = "flatly"),
      enable_admin = TRUE
    )
  } else {
    logger::log_warn("Running without authentication")
  }
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
