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

VALID_CONSEQUENCE <- c(
  "3_prime_UTR_variant",
  "5_prime_UTR_variant",
  "coding_sequence_variant",
  "downstream_gene_variant",
  "frameshift_variant",
  "inframe_deletion",
  "inframe_insertion",
  "intergenic_variant",
  "intron_variant",
  "mature_miRNA_variant",
  "missense_variant",
  "NMD_transcript_variant",
  "non_coding_transcript_exon_variant",
  "non_coding_transcript_variant",
  "protein_altering_variant",
  "regulatory_region_variant",
  "splice_acceptor_variant",
  "splice_donor_5th_base_variant",
  "splice_donor_region_variant",
  "splice_donor_variant",
  "splice_polypyrimidine_tract_variant",
  "splice_region_variant",
  "start_lost",
  "start_retained_variant",
  "stop_gained",
  "stop_lost",
  "stop_retained_variant",
  "synonymous_variant",
  "upstream_gene_variant"
)

VALID_CLINSIG <- c(
  "benign",
  "likely_benign",
  "pathogenic",
  "likely_pathogenic",
  "not_provided",
  "other",
  "",
  "risk_factor",
  "uncertain_significance"
)


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

filter_list_col <- function(tb, accepted, valid, column, sep = "&") {
  checkmate::assert_names(accepted, subset.of = valid)
  checkmate::assert_tibble(tb)
  if (class(tb[[column]]) != "list") {
    tb[[column]] <- lapply(
      tb[[column]],
      \(cons) stringr::str_split_1(cons, sep)
    )
  }
  mask <- purrr::map_lgl(tb[[column]], \(x) length(intersect(x, accepted)) >= 1)
  tb[mask, ]
}

combine_vcf_tbs <- function(tbs, filters = NULL) {
  box::use(dplyr[across, any_of, mutate, select])

  vep_cols <- c(
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
    "HGVSg",
    "ClinSig",
    "Somatic_VEP",
    "Pheno",
    "PubMed"
  )

  avg <- c("AF", "AD_REF", "AD_MAX", "AN", "GERMQ", "QSS", "SomaticEVS")
  uniq <- c("GT", "Sample")

  gene_cols <- c(
    "Symbol",
    "Gene",
    "Feature",
    "Feature_type",
    "Biotype",
    "Canonical"
  )

  cons_count <- c("synonymous_variant", "missense_variant")

  checkmate::assert_list(tbs, names = "unique")
  comb <- lapply(names(tbs), \(x) {
    tb <- dplyr::mutate(
      tbs[[x]],
      Sample = x,
      Consequence = lapply(Consequence, \(cons) str_split_1(cons, "&"))
    )
    if ("consequence" %in% names(filters)) {
      tb <- filter_list_col(
        tb,
        filters$consequence,
        VALID_CONSEQUENCE,
        "Consequence"
      )
    }
    if ("clinsig" %in% names(filters)) {
      tb <- filter_list_col(
        tb,
        filters$clinsig,
        VALID_CLINSIG,
        "ClinSig"
      )
    }

    ## for (cons in cons_count) {
    ##   tb[[paste0("is_", cons)]] <- purrr::pluck(tb$Consequence, \(x) x == cons)
    ## }

    tb
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(HGVSg) |>
    dplyr::summarise(
      `Sample count` = length(unique(Sample)),
      across(any_of(vep_cols), dplyr::first),
      across(any_of(avg), \(x) mean(x, na.rm = TRUE)),
      across(any_of(uniq), \(x) list(unique(x))),
      SOMATIC = any(SOMATIC),
      Consequence = list(unique(unlist(Consequence)))
    ) |>
    dplyr::relocate(
      dplyr::starts_with("HGVS"),
      .before = dplyr::everything()
    ) |>
    dplyr::mutate(
      HGVSc = stringr::str_extract(HGVSc, ".*c\\.(.*)", group = 1),
      HGVSp = stringr::str_extract(HGVSp, ".*p\\.(.*)", group = 1)
    )

  grouped <- comb |>
    group_by(across(all_of(gene_cols))) |>
    nest() |>
    mutate(
      `N variants` = nrow(data),
      `Consequence counts` = lapply(data, \(x) table(unlist(x$Consequence)))
    )

  list(genes = grouped, vars = comb)
}

#' Produce interactive reactable to explore variant data
#'
#' @param gene_tb tibble returned from `combine_vcf_tbs$grouped`
#'
make_variant_table <- function(gene_tb) {
  box::use(reactable[reactable, colGroup])
  reactable(
    dplyr::select(gene_tb, -data),
    details = function(index) {
      cur <- gene_tb$data[[index]]
      reactable(
        cur,
        outlined = TRUE,
        columnGroups = list(
          colGroup(
            name = "IDs",
            columns = c("HGVSg", "HGVSc", "HGVSp")
          ),
          colGroup(name = "Location", columns = c("Exon", "Intron")),
          colGroup(
            name = "Classification",
            columns = c("Consequence", "Impact", "ClinSig"),
          ),
          colGroup(
            name = "Cohort statistics",
            columns = c(
              "Sample count",
              "Sample",
              "AF",
              "AD_REF",
              "AD_MAX",
              "AN",
              "QSS",
              "GERMQ",
              "SOMATIC",
              "SomaticEVS",
              "GT"
            )
          ),
          colGroup(
            name = "External links",
            columns = c("Existing_variation", "PubMed")
          )
        )
      )
    }
  )
}
