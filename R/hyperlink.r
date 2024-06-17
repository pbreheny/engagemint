#' Create an html hyperlink

hyperlink <- function(link, text) {
  paste0("<a href='",
         link,
         "'target='_blank'>",
         text)
}
