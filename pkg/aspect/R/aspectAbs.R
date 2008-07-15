`aspectAbs` <-
function(r, pow = 1) {
  rl <- r[lower.tri(r)]
  list(f = sum(abs(rl)^pow), g = pow*sign(rl)*(abs(rl)^(pow-1)))
}

