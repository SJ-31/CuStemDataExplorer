#' compare_expression UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_compare_expression_ui <- function(id, label) {
  box::use(bslib[nav_panel], shiny[observeEvent])

  ns <- NS(id)
  bslib::nav_panel(
    label,
    bslib::layout_sidebar(
      bslib::navset_card_tab(
        nav_panel(
          "Heatmap",
          shiny::plotOutput(ns("expr_comparison")),
          value = "heatmap"
        ),
        nav_panel("PCA", shiny::plotOutput(ns("pca")), value = "pca"),
        id = ns("nav")
      ),
      sidebar = bslib::sidebar(
        shiny::conditionalPanel(
          condition = "input.nav != 'pca'",
          shiny::h4("Gene selection"),
          shiny::selectizeInput(
            ns("gene_selection"),
            label = NULL,
            choices = NULL,
            multiple = TRUE,
            options = list(maxItems = 15)
          ),
          ns = ns
        ),
        shiny::h4("Filters"),
        shiny::selectizeInput(
          ns("tumor_type"),
          "Tumor Type",
          choices = NULL,
          multiple = TRUE,
          options = list(items = NULL)
        ),
        shiny::selectizeInput(
          ns("cohort"),
          "Cohort",
          choices = NULL,
          multiple = TRUE,
          options = list(items = NULL)
        ),
        shiny::conditionalPanel(
          condition = "input.nav != 'pca'",
          shiny::h4("Aggregation"),
          shiny::selectInput(
            ns("group_by"),
            label = NULL,
            choices = c(
              "---" = "none",
              "Tumor type" = "tumor_type",
              "Cohort" = "cohort"
            ),
            selected = "---"
          ),
          ns = ns,
        ),
        shiny::conditionalPanel(
          condition = "input.nav == 'pca'",
          shiny::h4("Color by"),
          shiny::selectizeInput(
            ns("color_by"),
            label = NULL,
            choices = c(),
          ),
          ns = ns,
        ),
        shiny::h4("Reset palettes"),
        shiny::selectizeInput(
          ns("palette_choices"),
          label = NULL,
          choices = NULL,
          multiple = TRUE
        ),
        shiny::actionButton(ns("reset_palette"), "Reset")
      )
    )
  )
}

#' compare_expression Server Functions
#'
#' @noRd
mod_compare_expression_server <- function(id, cached) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cfg <- get_golem_config("expression_viewer")
    combined_expr <- from_bfc(cached)
    pca_key <- paste0(cached, "_pca")
    pca <- from_bfc(pca_key)
    gene_ids <- purrr::discard(combined_expr$expr$gene_id, is.na) |>
      `names<-`(NULL)

    shiny::updateSelectizeInput(
      session,
      "gene_selection",
      choices = gene_ids,
      server = TRUE
    )

    shiny::updateSelectizeInput(
      session,
      "color_by",
      choices = purrr::discard(
        colnames(combined_expr$meta),
        \(x) x %in% c("sample")
      )
    )

    shiny::updateSelectizeInput(
      session,
      "tumor_type",
      choices = unique(combined_expr$meta$tumor_type),
      server = TRUE
    )
    shiny::updateSelectizeInput(
      session,
      "cohort",
      choices = unique(combined_expr$meta$cohort),
      server = TRUE
    )

    # Randomize whenever `input$reset_palette` gets pressed
    shiny::observe({
      if (input$nav != "pca") {
        randomize_palette(
          input$palette_choices,
          key_group = cached
        )
      } else {
        randomize_palette(
          input$palette_choices,
          key_group = pca_key
        )
      }
    }) |>
      shiny::bindEvent(input$reset_palette)

    cur_palettes <- shiny::reactiveVal(rlang::hash(CACHE$get(cached)))
    cur_palettes_pca <- shiny::reactiveVal(rlang::hash(CACHE$get(pca_key)))

    shiny::observeEvent(input$reset_palette, {
      new <- rlang::hash(CACHE$get(cached))
      new_pca <- rlang::hash(CACHE$get(pca_key))
      cur_palettes_pca(new_pca)
      cur_palettes(new)
    })

    output$expr_comparison <- shiny::renderPlot(
      {
        input$reset_palette
        update_palette_input(
          session,
          input,
          "palette_choices",
          key_group = cached
        )
        do_expr_plot(
          combined_expr,
          genes = input$gene_selection,
          cfg = cfg$palette %||% list(),
          tumor_types = input$tumor_type,
          key_group = cached,
          cohorts = input$cohort,
          group_by = input$group_by
        )
      },
      res = 120
    ) |>
      shiny::bindCache(
        input$gene_selection,
        input$tumor_type,
        input$cohort,
        input$group_by,
        cur_palettes()
      ) |>
      shiny::bindEvent(
        input$gene_selection,
        input$tumor_type,
        input$cohort,
        input$group_by,
        cur_palettes()
      )

    output$pca <- shiny::renderPlot(
      {
        input$reset_palette
        update_palette_input(
          session,
          input,
          "palette_choices",
          key_group = pca_key
        )
        plot_pca(
          pca,
          combined_expr,
          color_by = input$color_by,
          tumor_types = input$tumor_type,
          cohorts = input$cohort,
          key_group = pca_key
        )
      },
      res = 120
    ) |>
      shiny::bindCache(
        input$tumor_type,
        input$cohort,
        input$color_by,
        cur_palettes_pca()
      ) |>
      shiny::bindEvent(
        input$tumor_type,
        input$cohort,
        input$color_by,
        cur_palettes_pca()
      )
  })
}

## To be copied in the UI
# mod_compare_expression_ui("compare_expression_1")

## To be copied in the server
# mod_compare_expression_server("compare_expression_1")
