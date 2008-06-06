`aspectDeterminant` <-
function(r,extra) {
  list(f = -log(det(r)), g = -solve(r))
}

