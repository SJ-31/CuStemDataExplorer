html_join_newlines <- function(vec, wrap_pre = FALSE) {
  if (wrap_pre) {
    paste0("<pre>", vec, "</pre>", collapse = "<br>")
  } else {
    paste0(vec, collapse = "<br>")
  }
}

html_center_div <- function(content) {
  paste0("<div style='text-align: center'>", content, "</div>")
}
