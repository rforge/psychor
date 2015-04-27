scores <- function(x) UseMethod("scores")

## extract the scores for the first dimension
scores.gifi <- function(x) {
  x$scoremat[,,1]
}  