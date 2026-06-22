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
  box::use(VariantAnnotation[info], dplyr[rename, any_of])
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
  tb |>
    dplyr::select(-any_of("TRANSCRIPTION_FACTORS")) |>
    rename(Somatic_VEP = "SOMATIC") |>
    dplyr::rename_with(
      stringr::str_to_title,
      any_of(c(
        "SYMBOL",
        "IMPACT",
        "EXON",
        "BIOTYPE",
        "INTRON",
        "CANONICAL",
        "FLAGS",
        "PHENO"
      ))
    ) |>
    rename(PubMed = "PUBMED", ClinSig = "CLIN_SIG")
}

#' Read a VCF file and format into a tibble for later visualization
#'
format_sample_vcf <- function(
  vcf_file,
  sample_name,
  filters = NULL,
  wanted_cols = c(
    "AN",
    "GERMQ",
    "MBQ",
    "MMQ",
    "MQ",
    "QSS",
    "SOMATIC",
    "SomaticEVS",
    "TLOD",
    "Consequence",
    "Impact",
    "Symbol",
    "Gene",
    "Feature_type",
    "Feature",
    "Biotype",
    "Exon",
    "Intron",
    "HGVSc",
    "HGVSp",
    "Existing_variation",
    "Canonical",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "HGVSg",
    "ClinSig",
    "Somatic_VEP",
    "Pheno",
    "PubMed"
  )
) {
  box::use(VariantAnnotation[geno, info], tibble[tibble])
  vcf <- VariantAnnotation::readVcf(vcf_file)
  # Apply filters based on vcf statistics
  from_geno <- tibble(
    AF = geno(vcf)$AF[, sample_name],
    AD = geno(vcf)$AD[, sample_name],
    GT = geno(vcf)$GT[, sample_name]
  ) |>
    tidyr::hoist(
      "AD",
      AD_REF = 1,
      .remove = TRUE
    ) |>
    dplyr::mutate(
      AD_MAX = purrr::map_dbl(AD, max),
      AD = purrr::map_chr(AD, \(x) paste0(as.character(x), collapse = ","))
    )

  if (!is.null(filters)) {
    if ("min_alt_depth" %in% names(filters)) {
      mask <- from_geno$AD_MAX >= filters$min_alt_depth
      vcf <- vcf[mask, ]
      from_geno <- from_geno[mask, ]
    }
    if ("max_ref_depth" %in% names(filters)) {
      mask <- from_geno$AD_REF < filters$max_ref_depth
      vcf <- vcf[mask, ]
      from_geno <- from_geno[mask, ]
    }
    if ("min_af" %in% names(filters)) {
      mask <- from_geno$AF >= filters$min_af
      vcf <- vcf[mask, ]
      from_geno <- from_geno[mask, ]
    }
  }

  anno <- get_vep_anno(vcf)

  # Deduplicate
  granges <- MatrixGenerics::rowRanges(vcf)
  dupes <- duplicated(tibble(
    REF = as.character(granges$REF),
    ALT = S4Vectors::unstrsplit(as(granges$ALT, "CharacterList"), ",")
  )) &
    duplicated(granges)

  anno <- anno[!dupes, ]
  vcf <- vcf[!dupes, ]
  from_geno <- from_geno[!dupes, ]

  info(vcf)$ANN <- NULL
  info(vcf) <- suppressMessages(cbind(info(vcf), anno))

  variant_counts <- lapply(info(vcf)$Consequence, \(x) str_split_1(x, "&")) |>
    unlist() |>
    table()

  tb <- dplyr::bind_cols(
    from_geno,
    dplyr::select(as.data.frame(info(vcf)), dplyr::any_of(wanted_cols))
  )

  list(tb = tb, sbs_counts = count_sbs(vcf), variant_counts = variant_counts)
}

combine_vcf_tbs <- function(tbs) {
  box::use(dplyr[across, any_of])

  vep_cols <- c(
    "Consequence",
    "Impact",
    "Symbol",
    "Gene",
    "Feature_type",
    "Feature",
    "Biotype",
    "Exon",
    "Intron",
    "HGVSc",
    "HGVSp",
    "Existing_variation",
    "Canonical",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "HGVSg",
    "ClinSig",
    "Somatic_VEP",
    "Pheno",
    "PubMed"
  )

  avg <- c("AF", "AD_REF", "AD_MAX", "AN", "GERMQ", "QSS", "SomaticEVS")
  uniq <- c("GT", "Sample")

  checkmate::assert_list(tbs, names = "unique")
  lapply(names(tbs), \(x) dplyr::mutate(tbs[[x]], Sample = x)) |>
    dplyr::bind_rows() |>
    dplyr::group_by(HGVSg, Consequence) |>
    dplyr::summarise(
      `Sample count` = length(unique(Sample)),
      across(any_of(vep_cols), dplyr::first),
      across(any_of(avg), \(x) mean(x, na.rm = TRUE)),
      across(any_of(uniq), \(x) list(unique(x))),
      SOMATIC = any(SOMATIC)
    )
}
