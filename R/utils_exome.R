count_sbs <- function(vcf) {
  box::use(dplyr[replace_values, case_when, rename])
  subs <- mcols(vcf)[
    VariantAnnotation::isSubstitution(vcf) &
      VariantAnnotation::isSNV(vcf),
  ]
  refs <- as.character(subs$REF)
  split_alt <- S4Vectors::unstrsplit(subs$ALT) |> as.character()

  paste0(
    replace_values(refs, "A" ~ "T", "G" ~ "C"),
    ">",
    case_when(
      stringr::str_starts(refs, "A|G") ~
        replace_values(
          split_alt,
          "A" ~ "T",
          "C" ~ "G",
          "G" ~ "C",
          "T" ~ "A"
        ),
      .default = split_alt
    )
  ) |>
    table() |>
    as.data.frame() |>
    rename(substitution = "Var1", count = "Freq")
}

get_vep_anno <- function(vcf, info_field = "ANN") {
  box::use(VariantAnnotation[info])
  info <- VariantAnnotation::header(vcf) |> info()
  if (info_field %notin% rownames(info)) {
    stop("The provided info tag isn't present in the VCF file")
  }
  colnames <- info[rownames(info) == info_field, ]$Description |>
    stringr::str_extract("Format: (.*)", group = 1) |>
    stringr::str_split_1("\\|")
  ann <- as.matrix(info(vcf)[[info_field]])[, 1]
  tb <- read.table(
    textConnection(ann),
    sep = "|",
    col.names = colnames
  ) |>
    tibble::as_tibble()
  tb |> dplyr::select(-dplyr::any_of("TRANSCRIPTION_FACTORS"))
}
