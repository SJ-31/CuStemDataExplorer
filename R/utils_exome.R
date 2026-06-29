count_sbs <- function(vcf) {
  box::use(dplyr[replace_values, case_when, rename])
  subs <- S4Vectors::mcols(vcf)[
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
    "PubMed"
  )
) {
  box::use(
    VariantAnnotation[geno, info, `info<-`, header, samples],
    tibble[tibble]
  )
  vcf <- VariantAnnotation::readVcf(vcf_file)
  samples <- Biobase::samples(header(vcf))
  if (sample_name %notin% samples) {
    stop(sprintf(
      "Provided sample %s not present in VCF file. Available samples: %s",
      sample_name,
      paste0(samples, collapse = ", ")
    ))
  }
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
    if (
      "min_alt_depth" %in% names(filters) && !is.null(filters$min_alt_depth)
    ) {
      mask <- from_geno$AD_MAX >= filters$min_alt_depth
      vcf <- vcf[mask, ]
      from_geno <- from_geno[mask, ]
    }
    if (
      "max_ref_depth" %in% names(filters) && !is.null(filters$max_ref_depth)
    ) {
      mask <- from_geno$AD_REF < filters$max_ref_depth
      vcf <- vcf[mask, ]
      from_geno <- from_geno[mask, ]
    }
    if ("min_af" %in% names(filters) && !is.null(filters$min_af)) {
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

  variant_counts <- lapply(
    info(vcf)$Consequence,
    \(x) stringr::str_split_1(x, "&")
  ) |>
    unlist() |>
    table() |>
    as.data.frame() |>
    dplyr::rename(Consequence = "Var1", Count = "Freq")

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
    "Canonical",
    "HGVSg",
    "ClinSig"
  )

  avg <- c("AF", "AD_REF", "AD_MAX", "AN", "GERMQ", "QSS", "SomaticEVS")
  uniq <- c("GT", "Sample")
  list_uniq <- c("Consequence", "Existing_variation", "PubMed")

  gene_cols <- c(
    "Symbol",
    "Gene",
    "Feature",
    "Feature type",
    "Biotype",
    "Canonical"
  )

  checkmate::assert_list(tbs, names = "unique")
  comb <- lapply(names(tbs), \(x) {
    tb <- mutate(tbs[[x]], Sample = x) |>
      mutate(across(
        any_of(c("Consequence", "Existing_variation", "PubMed")),
        \(col) {
          lapply(
            col,
            \(e) {
              if (nchar(e) == 0) {
                character(0)
              } else {
                stringr::str_split_1(e, "&")
              }
            }
          )
        }
      ))
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
      across(any_of(list_uniq), \(x) list(unique(unlist(x))))
    ) |>
    dplyr::relocate(
      dplyr::starts_with("HGVS"),
      .before = dplyr::everything()
    ) |>
    mutate(
      HGVSc = stringr::str_extract(HGVSc, ".*c\\.(.*)", group = 1),
      HGVSp = stringr::str_extract(HGVSp, ".*p\\.(.*)", group = 1)
    ) |>
    dplyr::rename(
      `Existing variation` = "Existing_variation",
      `Feature type` = "Feature_type"
    )

  stats <- select(comb, any_of(c("HGVSg", avg, "SOMATIC", "GT"))) |>
    mutate(
      across(any_of(avg), as.character),
      SOMATIC = as.character(SOMATIC),
      GT = purrr::map_chr(GT, \(gt) {
        val <- purrr::discard(unique(gt), is.na) |> paste0(collapse = ", ")
        if (nchar(val) == 0) {
          "NaN"
        } else {
          val
        }
      })
    ) |>
    tidyr::pivot_longer(-HGVSg) |>
    mutate(name = paste0("<strong>", name, ":</strong>")) |>
    tidyr::unite(
      col = "Caller statistics",
      name,
      value,
      sep = "  ",
      remove = TRUE
    ) |>
    dplyr::group_by(HGVSg) |>
    dplyr::summarise(`Caller statistics` = list(`Caller statistics`))

  comb <- select(comb, -any_of(c(avg, "SOMATIC", "GT"))) |>
    dplyr::inner_join(stats, by = dplyr::join_by(HGVSg)) |>
    mutate(
      `Exon/Intron` = purrr::map2_chr(Exon, Intron, \(ex, int) {
        if ((is.na(int) && is.na(ex)) || (nchar(int) == 0 && nchar(ex) == 0)) {
          ""
        } else if (!is.na(ex) && nchar(ex) > 0) {
          glue::glue("{ex} (exon)")
        } else {
          glue::glue("{int} (intron)")
        }
      })
    ) |>
    select(-Exon, -Intron)

  grouped <- comb |>
    dplyr::group_by(across(dplyr::all_of(gene_cols))) |>
    tidyr::nest() |>
    mutate(
      N = purrr::map_dbl(data, nrow),
      `Consequence counts` = lapply(data, \(x) table(unlist(x$Consequence))),
      Canonical = dplyr::recode_values(Canonical, "YES" ~ "Y", default = "N")
    )

  list(genes = grouped, vars = comb)
}

GENE_VARIANT_COLS_DOWNLOAD <- c(
  "HGVSp",
  "HGVSg",
  "HGVSc",
  "Sample count",
  "Sample",
  "PubMed",
  "Consequence",
  "Existing variation",
  "Exon/Intron"
)

variant_table_download_button <- function(title, elementId, file = "data.csv") {
  cols <- c(
    "Symbol",
    "Gene",
    "Feature type",
    "Feature",
    "Biotype",
    "Canonical",
    "N"
  ) |>
    to_js_array()
  onclick <- sprintf(
    "Reactable.downloadDataCSV('%s', '%s', {
columnIds: %s
})",
    elementId, # Pass the outputId of the reactable
    file,
    cols
  )
  htmltools::tags$button(title, onclick = onclick)
}

#' Produce interactive reactable to explore variant data
#'
#' @param gene_tb tibble returned from `combine_vcf_tbs$grouped`
#'
make_variant_table <- function(gene_tb) {
  box::use(reactable[reactable, colGroup, colDef])
  reactable(
    tibble::rowid_to_column(dplyr::select(gene_tb, -data), "Download"),
    searchable = TRUE,
    columns = list(
      Download = colDef(
        align = "center",
        cell = \(v, i, n) {
          table_name <- sprintf("nested_%s", v)
          onclick <- sprintf(
            "Reactable.downloadDataCSV('%s', 'data.csv', {
columnIds: %s
})",
            table_name,
            to_js_array(GENE_VARIANT_COLS_DOWNLOAD)
          )
          htmltools::tags$button(
            "📥",
            onclick = onclick,
            style = "background-color:white;border-color:white;box-shadow:none"
          )
        }
      ),
      `Consequence counts` = colDef(
        cell = \(v, i, n) sum(v),
        html = TRUE,
        details = \(i) {
          val <- gene_tb$`Consequence counts`[[i]]
          paste0(
            "<strong>",
            names(val),
            ": </strong>",
            val,
            collapse = "<br>"
          ) |>
            html_pad_div()
        }
      )
    ),
    details = function(index) {
      cur <- gene_tb$data[[index]]
      reactable(
        cur,
        compact = TRUE,
        elementId = sprintf("nested_%s", index),
        searchable = TRUE,
        columns = list(
          `Sample count` = colDef(html = TRUE, details = \(i) {
            html_join_newlines(cur$Sample[[i]]) |>
              html_with_header("Sample names") |>
              html_pad_div()
          }),
          Sample = colDef(show = FALSE),
          `Existing variation` = colDef(
            html = TRUE,
            cell = \(v, i, n) length(v),
            details = \(i) reactable_display_list(i, cur, "Existing variation")
          ),
          PubMed = colDef(
            html = TRUE,
            cell = \(v, i, n) length(v),
            details = \(i) reactable_display_list(i, cur, "PubMed")
          ),
          `Caller statistics` = colDef(
            html = TRUE,
            cell = \(v, i, n) "",
            details = \(i) reactable_display_list(i, cur, "Caller statistics")
          ),
          Consequence = colDef(
            html = TRUE,
            cell = \(v, i, n) length(v),
            details = \(i) reactable_display_list(i, cur, "Consequence")
          )
        ),
        columnGroups = list(
          colGroup(
            name = "Change",
            columns = c("HGVSg", "HGVSc", "HGVSp")
          ),
          colGroup(
            name = "Classification",
            columns = c("Consequence", "Impact", "ClinSig"),
          ),
          colGroup(
            name = "External links",
            columns = c("Existing variation", "PubMed")
          )
        ),
        bordered = TRUE,
        theme = reactable::reactableTheme(
          backgroundColor = "#c3e7eb",
          headerStyle = list(backgroundColor = "#559f9f", color = "#feffff"),
          borderColor = "#feffff",
          groupHeaderStyle = list(
            backgroundColor = "#feffff",
            color = "#509d9d"
          ),
          searchInputStyle = list(width = "100%")
        )
      )
    },
    borderless = TRUE,
    theme = reactable::reactableTheme(
      borderColor = "#7eccd3"
    )
  )
}

#' Read a variant specification file
#'
#' @return
#' A list with names corresponding to tumor types. Each
#' value is a named list with four tibbles: genes, vars, sbs_counts, var_counts
#' The former two are aggregated and single variant statistics, respectively
read_variant_spec <- function(file, cfg) {
  box::use(dplyr[mutate])
  contents <- yaml::read_yaml(file)

  read <- function(spec) {
    to_combine <- list()
    from_cohort <- lapply(spec, \(lst) {
      cohort <- lst$cohort %||% "-"
      prefix <- lst$prefix %||% ""
      suffix <- lst$suffix %||% ""
      snames <- names(lst$samples)

      vc_filters <- cfg$filters$variant_calling

      res <- lapply(snames, \(n) {
        vcf_key <- glue::glue("{prefix}{n}{suffix}")
        tmp <- format_sample_vcf(
          lst$samples[[n]],
          sample_name = vcf_key,
          filters = vc_filters
        )
        if (!is.null(lst$cohort)) {
          combine_key <- glue::glue("{cohort}::{vcf_key}")
        } else {
          combine_key <- vcf_key
        }
        to_combine[[combine_key]] <<- tmp$tb

        tmp$tb <- NULL
        tmp$sbs_counts$sample <- n
        tmp$sbs_counts$cohort <- cohort
        tmp$variant_counts$cohort <- cohort
        tmp$variant_counts$sample <- n
        tmp
      }) |>
        `names<-`(snames)

      sbs <- dplyr::bind_rows(lapply(res, \(x) x$sbs_counts))
      var_counts <- dplyr::bind_rows(lapply(res, \(x) x$variant_counts))
      list(s = sbs, v = var_counts)
    })
    combined <- combine_vcf_tbs(to_combine, filters = cfg$filters)
    all_sbs <- lapply(from_cohort, \(x) x$s) |> dplyr::bind_rows()
    all_var_counts <- lapply(from_cohort, \(x) x$v) |> dplyr::bind_rows()

    list(
      genes = combined$genes,
      vars = combined$vars,
      var_counts = all_var_counts,
      sbs_counts = all_sbs
    )
  }

  lapply(contents, read) |> `names<-`(names(contents))
}

#' Read sample specification files for variant data
#'
#' @description
#' See the README for the description of the spec files
read_variant_spec_all <- function(cfg, allowed_exts = c("vcf", "vcf.gz")) {
  dir <- cfg$spec_directory
  checkmate::assert_directory_exists(dir)
  per_file <- lapply(
    list.files(dir, pattern = ".yml|yaml$", full.names = TRUE),
    \(f) read_variant_spec(f, cfg)
  )
  per_ttypes <- list(
    vars = list(),
    genes = list(),
    sbs_counts = list(),
    var_counts = list()
  )
  for (lst in per_file) {
    for (tb_name in names(per_ttypes)) {
      for (ttype in names(lst)) {
        prev <- per_ttypes[[tb_name]][[ttype]]
        cur <- lst[[ttype]][[tb_name]]
        per_ttypes[[tb_name]][[ttype]] <- dplyr::bind_rows(prev, cur)
      }
    }
  }

  per_ttypes
}
