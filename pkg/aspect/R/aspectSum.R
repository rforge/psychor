#sum of the r^pow elements of the correlation matrix 

`aspectSum` <-
function(r, pow = 1)                 
{
  rl <- r[lower.tri(r)]
  list(f = sum(rl^pow), g = pow*rl^(pow-1))
}

