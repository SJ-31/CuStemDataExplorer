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

html_with_header <- function(
  content,
  header,
  bolden = TRUE,
  wrapper = identity
) {
  if (bolden) {
    head <- wrapper(htmltools::strong(header))
  } else {
    head <- wrapper(header)
  }
  paste0(head, "<br>", content)
}

reactable_display_list <- function(index, table, col, center = TRUE) {
  val <- table[[col]][[index]]
  if (length(val) > 0) {
    joined <- html_join_newlines(val)
    if (center) {
      html_center_div(joined)
    } else {
      joined
    }
  }
}
