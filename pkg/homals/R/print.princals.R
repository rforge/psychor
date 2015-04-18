print.princals <- function(x, ...)
{
    nvar <- dim(x$dframe)[2]
    cat("Call: ")
    print(x$call)
    cat("\nLoss:",x$loss)
    cat("\nNumber of iterations:" ,x$niter, "\n")
    cat("\nEigenvalues:")
    eigen.val <- round(x$eigenvalues, 5)
    cat("\n")
    names(eigen.val) <- paste("Comp.", 1:length(eigen.val), sep = "")
    print(eigen.val)
    cat("\n")
}