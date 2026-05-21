#' Access files in the current app
#'
#' NOTE: If you manually change your package name in the DESCRIPTION,
#' don't forget to change it here too, and in the config file.
#' For a safer name change mechanism, use the `golem::set_golem_name()` function.
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#' within your package. The default, none, returns the root of the app.
#'
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "StemDataExplorer")
}

golem_config_env <- function() {
  env_lookup <- Sys.getenv("GOLEM_CONFIG_FILE")
  if (nchar(env_lookup) > 0) {
    if (file.exists(env_lookup)) {
      env_lookup
    } else {
      logger::log_warn(
        "Configuration file {env_lookup} specified by GOLEM_CONFIG_FILE doesn't exist! Using only defaults"
      )
      NULL
    }
  } else {
    NULL
  }
}

#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, R_CONFIG_ACTIVE.
#' If unset, "default".
#' @param use_parent Logical, scan the parent directory for config file.
#' @param file Location of the config file
#'
#' @noRd
get_golem_config <- function(
  value = NULL,
  config = Sys.getenv(
    "GOLEM_CONFIG_ACTIVE",
    Sys.getenv(
      "R_CONFIG_ACTIVE",
      "default"
    )
  ),
  use_parent = TRUE,
  # Modify this if your config file is somewhere else
  file = app_sys("golem-config.yml")
) {
  defaults <- config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
  from_env <- golem_config_env()
  if (!is.null(from_env)) {
    modified <- config::get(
      value = value,
      config = config,
      file = from_env,
      use_parent = TRUE
    )
    config::merge(defaults, modified)
  } else {
    defaults
  }
}
