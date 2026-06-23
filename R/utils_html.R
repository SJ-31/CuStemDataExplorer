html_join_newlines <- function(vec, wrap_pre = FALSE) {
  if (wrap_pre) {
    paste0("<pre>", vec, "</pre>", collapse = "<br>")
  } else {
    paste0(vec, collapse = "<br>")
  }
}

html_center_div <- function(content) {
  ## paste0("<div style='text-align: center'>", content, "</div>")
  paste0("<div style='text-align-last: center'>", content, "</div>")
}

html_pad_div <- function(content, pad = "20%") {
  paste0(
    "<div style='padding-left:",
    pad,
    "'>",
    content,
    "</div>"
  )
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
      html_pad_div(joined)
    } else {
      joined
    }
  }
}
