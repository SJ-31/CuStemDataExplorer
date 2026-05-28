rename_df <- function(df) {
  lookup <- c(
    Case = "case_name",
    Path = "path",
    `N biopsies` = "n_biopsies",
    PBMC = "has_pbmc",
    Tumor = "has_tumor",
    `Processed files` = "processed_files",
    Warnings = "warnings",
    Raw = "has_raw",
    Processed = "has_processed",
    Modality = "modality",
    Diagnosis = "diagnosis",
    `NGS platform` = "platform",
    `NGS provider` = "ngs_provider",
    `Sample collection` = "sample_collection_date",
    `Received data` = "received_date",
    `Tumor Type` = "tumor_type",
    Cohort = "cohort",
    Note = "note"
  )
  df |> dplyr::rename(dplyr::any_of(lookup))
}

read_from_other <- function(link, sheet_name, grouped_sample_df) {
  googlesheets4::read_sheet(link, sheet = sheet_name) |>
    dplyr::mutate(
      tumor_type = stringr::str_to_upper(tumor_type),
      case_name = dplyr::first(case_name),
    ) |>
    dplyr::left_join(
      grouped_sample_df,
      by = c(
        "tumor_type",
        "case_name",
        "cohort"
      )
    ) |>
    dplyr::mutate(`Has data?` = ifelse(is.na(Modality), "F", "T")) |>
    dplyr::select(-Modality)
}

#' get_sheets
#'
#' @description Retrieve data availability, clinical, and sample metadata
#'
#' @return A list of tibbles containing the sample data
#'
#' @noRd
get_sheets <- function(cfg) {
  data <- list()
  sheets <- cfg$data_sheets
  other_sheet <- sheets$other
  if (file.exists(sheets$samples)) {
    logger::log_info("Reading from on-disk sample sheet")
    sample_sheet <- readr::read_csv(sheets$samples)
  } else {
    sample_sheet <- googlesheets4::read_sheet(sheets$samples)
  }
  df <- sample_sheet |>
    dplyr::mutate(
      tumor_type = stringr::str_to_upper(tumor_type),
      modality = dplyr::replace_values(
        modality,
        "exome" ~ "Exome",
        "rna_seq" ~ "RNA seq",
        "scrna_seq" ~ "scRNA seq",
        "sc_atac_seq" ~ "scATAC seq",
        "tcr_seq" ~ "TCR seq"
      ),
    ) |>
    dplyr::relocate(path, .after = dplyr::everything())
  grouped <- df |>
    dplyr::group_by(cohort, case_name, tumor_type) |>
    dplyr::summarise(Modality = paste0(modality, collapse = "; "))
  data$all <- rename_df(df) |> crosstalk::SharedData$new(group = "tables")
  data$clinical <- read_from_other(other_sheet, "clinical", grouped) |>
    rename_df() |>
    crosstalk::SharedData$new(group = "tables")
  data$meta <- read_from_other(other_sheet, "metadata", grouped) |>
    rename_df() |>
    crosstalk::SharedData$new(group = "tables")
  data
}

get_sheet_sidebar <- function(sheets) {
  bslib::sidebar(
    shiny::h3("Filters"),
    crosstalk::filter_checkbox("modality", "Modality", sheets$all, ~Modality),
    crosstalk::filter_select("cohort", "Cohort", sheets$all, ~Cohort),
    crosstalk::filter_checkbox(
      "ttype",
      "Tumor Type",
      sheets$all,
      ~`Tumor Type`
    ),
  )
}
