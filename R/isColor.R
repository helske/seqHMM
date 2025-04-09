# Copied from a post by Josh O'Brien on Stack Overflow:
# https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
isColor <- function(x) {
  sapply(x, \(X) {
    tryCatch(is.matrix(grDevices::col2rgb(X)), error = function(e) FALSE)
  })
}
