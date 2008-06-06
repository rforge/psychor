`aspectAbs` <-
function(r, pow) {
  list(f = sum(abs(r)^pow)/2, g = pow*sign(r)*(abs(r)^(pow-1)))
}

