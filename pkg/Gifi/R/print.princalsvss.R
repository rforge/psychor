print.princalsvss <- function(x, ...) {
  
  cat("Call: ")
  print(x$call)
  
  nvss1 <- which.max(x$VSS1)         ## VSS complexity 1
  vss1 <- x$VSS1[nvss1]
  nvss2 <- which.max(x$VSS2)         ## VSS complexity 2
  vss2 <- x$VSS2[nvss2]
  nmap <- which.min(x$MAP)           ## MAP
  map <- x$MAP[nmap]
  
  cat("\n")
  if(nvss1 == 1) {
    cat("VSS complexity 1 achieves a maximimum of", round(vss1, 2), "with", nvss1,"factor.")
  } else {
    cat("VSS complexity 1 achieves a maximimum of", round(vss1, 2), "with", nvss1,"factors.")
  }
  
  if(nvss2 == 1) {
    cat("\nVSS complexity 2 achieves a maximimum of", round(vss2, 2), "with", nvss2,"factor.\n")
  } else {
    cat("\nVSS complexity 2 achieves a maximimum of", round(vss2, 2), "with", nvss2,"factors.\n")
  }
  
  if(nmap == 1) {
    cat("\nThe MAP achieves a minimum of", round(map, 2), "with", nmap,"factor.\n")
  } else {
    cat("\nThe MAP achieves a minimum of", round(map, 2), "with", nmap,"factors.\n")
  }
  cat("\n")
}