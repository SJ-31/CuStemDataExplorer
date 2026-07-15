suppressMessages({
  library(shinymanager)
  library(here)
})

for_prod <- Sys.getenv("GOLEM_CONFIG_ACTIVE")

if (for_prod == "production") {
  credentials <- readr::read_csv(here("deploy", "users.csv "))
} else {
  credentials <- data.frame(
    user = c("admin", "default"),
    password = c("foobar", "foobar"),
    admin = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}

sqlite <- here("inst", "authen.sqlite")
create_db(credentials_data = credentials, sqlite_path = sqlite)
