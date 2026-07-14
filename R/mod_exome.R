#' exome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_exome_ui <- function(id, label) {
  ns <- NS(id)
  bslib::nav_panel(
    label,
    bslib::layout_sidebar(
      reactable::reactableOutput(ns("table")),
      reactable::reactableOutput(ns("vtable")),
      sidebar = bslib::sidebar(
        shiny::h3("Tumor type"),
        shiny::selectizeInput(
          ns("ttype"),
          label = NULL,
          choices = NULL,
          multiple = FALSE
        ),
        variant_table_download_button("Download", ns("table")),
        shiny::br(),
        shiny::h3("Selected gene:"),
        shiny::htmlOutput(ns("shown_gene")),
        shiny::conditionalPanel(
          condition = "output.shown_gene",
          variant_table_download_button(
            "Download",
            ns("vtable"),
            cols = GENE_VARIANT_COLS_DOWNLOAD
          ),
          ns = ns
        )
      ),
      fillable = FALSE
    )
  )
}

#' exome Server Functions
#'
#' @noRd
mod_exome_server <- function(id, cached) {
  box::use(DBI[dbGetQuery])
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    con <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = from_bfc(cached)$genes,
      read_only = TRUE
    )
    ttypes <- dbGetQuery(
      con,
      "SELECT * FROM INFORMATION_SCHEMA.TABLES"
    )$table_name

    shiny::updateSelectizeInput(
      inputId = "ttype",
      choices = ttypes,
      selected = ttypes[1],
      server = TRUE
    )

    make_table <- function(type) {
      data <- dbGetQuery(con, sprintf("SELECT * FROM %s", type))
      make_variant_table_at_level(data, level = "genes")
    }

    chosen_row <- reactive({
      req(input$ttype)
      idx <- reactable::getReactableState("table", name = "selected")
      req(idx)
      query <- sprintf(
        "SELECT * FROM %s OFFSET %s ROWS FETCH NEXT 1 ROWS ONLY",
        input$ttype,
        idx - 1
      )
      dbGetQuery(con, query)
    })

    make_vtable <- function(type, row) {
      data <- dbGetQuery(con, sprintf("SELECT * FROM %s", type))
      make_variant_table_at_level(
        data,
        level = "variants",
        chosen = row$Gene
      )
    }

    # Could also make a separate view of the variant-level details, rather than have it be a nested table
    # could make things faster

    output$table <- reactable::renderReactable({
      req(input$ttype)
      make_table(input$ttype)
    }) |>
      shiny::bindCache(input$ttype) |>
      shiny::bindEvent(input$ttype)

    output$shown_gene <- renderUI({
      shiny::h5(chosen_row()$Symbol)
    })

    output$vtable <- reactable::renderReactable({
      req(input$ttype)
      req(chosen_row())
      make_vtable(input$ttype, chosen_row())
    }) |>
      shiny::bindCache(input$ttype, chosen_row()) |>
      shiny::bindEvent(input$ttype, chosen_row())
  })
}

## To be copied in the UI
# mod_exome_ui("exome_1")

## To be copied in the server
# mod_exome_server("exome_1")
