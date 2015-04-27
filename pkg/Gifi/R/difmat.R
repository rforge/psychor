## utility function for monotone B-splines
difmat <- function (n) {
  m1 <- ifelse(outer(1:(n - 1),1:n,"-") == -1, 1, 0)
  m2 <- ifelse(outer(1:(n - 1),1:n,"-") == 0,-1, 0)
  return (m1 + m2)
}
