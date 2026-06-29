#' bulk_de UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_de_ui <- function(id, label) {
  ns <- NS(id)
  bslib::nav_panel(
    label,
    bslib::layout_columns(
      bslib::card(
        bslib::layout_columns(
          shiny::h3("Contrast"),
          shiny::actionButton(
            ns("reset_palette"),
            "Reset palette",
            style = "color:white;margin: 5px;vertical-align:baseline"
          ),
          htmltools::tags$button(
            "📥",
            onclick = sprintf(
              "Reactable.downloadDataCSV(\"%s\")",
              ns("de_results")
            )
          ),
          col_widths = c(7, 4, 1)
        ),
        shiny::selectizeInput(
          ns("contrast"),
          choices = NULL,
          label = NULL,
          multiple = FALSE
        ),
        reactable::reactableOutput(ns("de_results"))
      ),
      bslib::navset_card_tab(
        bslib::nav_panel(
          "Volcano plot",
          shiny::plotOutput(ns("volcano")),
          value = "volcano"
        ),
        bslib::nav_panel(
          "Scatter plot for selected genes (maximum 20)",
          shiny::plotOutput(ns("scatter")),
          value = "scatter"
        ),
        id = ns("nav")
      ),
      col_widths = c(5, 7)
    )
  )
}

#' de Server Functions
#'
#' @noRd
mod_de_server <- function(id, cached, expression_key) {
  box::use(
    reactable[colDef, colFormat, reactable],
    glue[glue],
    DBI[dbGetQuery]
  )

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cfg <- get_golem_config("de_viewer")

    expr <- from_bfc(expression_key)
    keys <- from_bfc(cached)
    dbs <- list(
      ovr = from_bfc(keys$ovr),
      pairwise = from_bfc(keys$pairwise)
    )

    table_lookup <- "SELECT * FROM INFORMATION_SCHEMA.TABLES"

    all_contrasts <- c(
      dbGetQuery(dbs$ovr, table_lookup)$table_name,
      dbGetQuery(dbs$pairwise, table_lookup)$table_name
    ) |>
      sort()

    get_contrast <- function(contrast, indices = NULL, format = TRUE) {
      if (nchar(contrast) == 0) {
        return(data.frame())
      }
      if (stringr::str_ends(contrast, "vs. Rest")) {
        con <- dbs$ovr
      } else {
        con <- dbs$pairwise
      }
      query <- glue("SELECT * FROM '{contrast}' ORDER BY log2FoldChange")
      df <- dbGetQuery(con, query) |> dplyr::filter(!is.na(log2FoldChange))
      if (format) {
        df <- df |>
          dplyr::rename(lfc = "log2FoldChange") |>
          dplyr::relocate(stat, .after = dplyr::everything()) |>
          dplyr::relocate(padj, .before = pvalue)
      }
      if (!is.null(indices)) {
        df[indices, ]
      } else {
        df
      }
    }

    to_numeric <- c("stat", "pvalue", "lfcSE", "baseMean", "lfc")
    de_col_format <- lapply(
      rep(1, length(to_numeric)),
      \(x) {
        reactable::colDef(format = reactable::colFormat(digits = 3))
      }
    ) |>
      `names<-`(to_numeric)

    shiny::updateSelectizeInput(
      session,
      "contrast",
      choices = all_contrasts,
      selected = all_contrasts[1],
      server = TRUE
    )
    ## * Update palettes

    palette_key_group <- "bulk_de"
    n_ttypes <- length(unique(expr$meta$tumor_type))

    shiny::observeEvent(input$reset_palette, {
      if (input$nav == "volcano") {
        randomize_palette(
          "volcano",
          key_group = palette_key_group,
          min_length = 2
        )
      } else {
        randomize_palette(
          "de_scatter",
          key_group = palette_key_group,
          min_length = n_ttypes
        )
      }
    })

    ## * Plots & tables
    output$volcano <- shiny::renderPlot({
      input$reset_palette
      volcano_plot(
        get_contrast(
          input$contrast,
          format = FALSE,
        ),
        key_group = "bulk_de"
      )
    }) |>
      shiny::bindCache(input$contrast, input$reset_palette) |>
      shiny::bindEvent(input$contrast, input$reset_palette)

    output$de_results <- reactable::renderReactable(reactable(
      get_contrast(
        input$contrast
      ),
      columns = de_col_format,
      selection = "multiple",
      searchable = TRUE,
      onClick = "select"
    )) |>
      shiny::bindCache(input$contrast) |>
      shiny::bindEvent(input$contrast)

    chosen_genes <- reactive({
      index <- reactable::getReactableState(
        "de_results",
        name = "selected"
      )
      get_contrast(input$contrast, indices = index)$gene
    })
    output$scatter <- shiny::renderPlot({
      input$reset_palette
      de_scatter_plot(
        expr_tbs = expr,
        chosen_genes = chosen_genes(),
        key_group = palette_key_group,
        contrast = input$contrast,
        factor = "tumor_type"
      )
    }) |>
      shiny::bindCache(input$contrast, chosen_genes(), input$reset_palette) |>
      shiny::bindEvent(input$contrast, chosen_genes(), input$reset_palette)
  })
}

## To be copied in the UI
# mod_de_ui("de_1")

## To be copied in the server
# mod_de_server("de_1")
