## uses the biplot function from stats

plot.vmu <- function(x, ...) {
  stats::biplot(x$conf.row, x$conf.col, ...)
}