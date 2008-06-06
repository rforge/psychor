`aspectSum` <-
function(r,pow)                 
{
  m <- dim(r)[1]
  list(f = sum(r^pow), g = pow*r^(pow-1))
}

