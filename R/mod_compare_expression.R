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
        shiny::downloadButton(ns("download_plot"), "Download plot"),
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
        shiny::conditionalPanel(
          condition = "input.nav != 'pca'",
          shiny::selectizeInput(
            ns("palette_choices"),
            label = NULL,
            choices = NULL,
            multiple = TRUE
          ),
          ns = ns
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
          input$color_by,
          key_group = pca_key
        )
      }
    }) |>
      shiny::bindEvent(input$reset_palette)

    cur_palettes <- shiny::reactiveVal(rlang::hash(CACHE$get(cached)))

    shiny::observeEvent(input$reset_palette, {
      new <- rlang::hash(CACHE$get(cached))
      cur_palettes(new)
    })

    ## ** Heatmap and boxplots
    comparison_plot_reactive <- reactive({
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
    }) |>
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
    output$expr_comparison <- renderPlot(comparison_plot_reactive(), res = 120)

    ## ** PCA

    shiny::updateSelectizeInput(
      session,
      "color_by",
      choices = purrr::discard(
        colnames(combined_expr$meta),
        \(x) x %in% c("sample")
      )
    )

    pca_plot_reactive <- reactive({
      input$reset_palette
      plot <- plot_pca(
        pca,
        combined_expr,
        color_by = input$color_by,
        tumor_types = input$tumor_type,
        cohorts = input$cohort,
        key_group = pca_key
      )
      plot
    }) |>
      shiny::bindCache(
        input$tumor_type,
        input$cohort,
        input$color_by,
        input$reset_palette
      ) |>
      shiny::bindEvent(
        input$tumor_type,
        input$cohort,
        input$color_by,
        input$reset_palette
      )
    output$pca <- shiny::renderPlot(pca_plot_reactive(), res = 120)

    ## ** Save plot
    output$download_plot <- shiny::downloadHandler(
      filename = \() {
        paste0(input$nav, ".pdf")
      },
      content = \(file) {
        if (input$nav != "pca") {
          width <- 8 +
            as.integer(log2(length(unique(combined_expr$meta$sample))))
          height <- 7 + as.integer(log2(length(input$gene_selection)))
          plot <- comparison_plot_reactive()
        } else {
          width <- 9
          height <- 9
          plot <- pca_plot_reactive()
        }
        suppressMessages(ggplot2::ggsave(
          filename = file,
          plot = plot,
          width = ifelse(is.na(width), 8, width),
          height = ifelse(is.na(height), 7, height)
        ))
      }
    )

    ## ** Update selection
    observeEvent(input$nav, {
      update_palette_input(
        session,
        input,
        input_id = "palette_choices",
        key_group = cached
      )
    })
  })
}

## To be copied in the UI
# mod_compare_expression_ui("compare_expression_1")

## To be copied in the server
# mod_compare_expression_server("compare_expression_1")
