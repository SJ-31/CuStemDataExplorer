#' stats UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_stats_ui <- function(id) {
  ns <- NS(id)
  bslib::nav_panel(
    "Overview",
    shiny::h2("Statistics"),
    shiny::textOutput(ns("description")),
    shiny::h2("Cached data"),
    reactable::reactableOutput(ns("cache_tab"))
  )
}

#' stats Server Functions
#'
#' @noRd
mod_stats_server <- function(id, sheets) {
  box::use(reactable[colDef, colFormat, reactable], glue[glue])

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    n_samples <- sheets

    names <- purrr::map_chr(get_cache_spec(), \(x) x$name)
    cache_latest <- BiocFileCache::bfcinfo(BFC) |>
      dplyr::filter(rname %in% names) |>
      dplyr::select(rname, create_time) |>
      dplyr::rename(Resource = "rname", `Last updated` = "create_time")
    cache_tab <- reactable(
      cache_latest,
      columns = list(
        `Last updated` = colDef(
          format = colFormat(date = TRUE)
        )
      )
    )

    output$description <- shiny::renderText({
      shiny::h4(glue("Number of samples: {99}"))
      glue(
        "Number of organoids"
      )
    })
    output$cache_tab <- reactable::renderReactable(cache_tab)
  })
}

## To be copied in the UI
# mod_stats_ui("stats_1")

## To be copied in the server
# mod_stats_server("stats_1")
