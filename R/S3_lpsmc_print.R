#' Print a lpsmc object.
#'
#' @description Print method for a \code{lpsmc} object.
#'
#' @param x An object of class \code{lpsmc}.
#' @param ... Further arguments to be passed to print routine.
#'
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @export

print.lpsmc <- function(x,...){

  tabcolnames <- c("Estimate","sd","CI90%.low","CI90.up%",
                   "CI95%.low","CI95%.up")
  K <- x$K

  # Table for the incidence part
  tabincidrow <- x$p
  tabincidout <- as.data.frame(matrix(0, nrow = tabincidrow, ncol = 6))
  rownames(tabincidout) <- colnames(x$X)
  colnames(tabincidout) <- tabcolnames
  tabincidout[,1] <- x$betahat
  tabincidout[,2] <- sqrt(diag(x$Covhat)[(K + 1):(K + tabincidrow)])
  tabincidout[,3] <- as.numeric(x$CI90[,1])[1:tabincidrow]
  tabincidout[,4] <- x$CI90[,2][1:tabincidrow]
  tabincidout[,5] <- x$CI95[,1][1:tabincidrow]
  tabincidout[,6] <- x$CI95[,2][1:tabincidrow]
  iscolnum <- sapply(tabincidout, is.numeric)
  tabincidout[iscolnum] <- lapply(tabincidout[iscolnum], round, 3)

  # Table for the latency part
  tablatencyrow <- x$q
  dimlat <- x$K+x$p+x$q
  tablatencyout <- as.data.frame(matrix(0, nrow = tablatencyrow, ncol = 6))
  rownames(tablatencyout) <- colnames(x$Z)
  colnames(tablatencyout) <- tabcolnames
  tablatencyout[,1] <- x$gammahat
  tablatencyout[,2] <- sqrt(diag(x$Covhat)[(K + tabincidrow + 1):dimlat])
  tablatencyout[,3] <- as.numeric(x$CI90[,1])[(tabincidrow + 1):(dimlat-K)]
  tablatencyout[,4] <- x$CI90[,2][(tabincidrow + 1):(dimlat-K)]
  tablatencyout[,5] <- x$CI95[,1][(tabincidrow + 1):(dimlat-K)]
  tablatencyout[,6] <- x$CI95[,2][(tabincidrow + 1):(dimlat-K)]
  iscolnum <- sapply(tablatencyout, is.numeric)
  tablatencyout[iscolnum] <- lapply(tablatencyout[iscolnum], round, 3)

  # Print output table

  cat("Fitting mixture cure model with Laplacian-P-splines \n")
  cat(paste(rep("-",50),collapse = ""),"\n")
  cat("Sample size: ", length(x$ftime), "\n")
  cat("No. of B-splines: ", K, "\n")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat("                                  (Incidence)                    \n")
  cat(paste(rep("-",90),collapse = ""),"\n")
  print.table(as.matrix(tabincidout), digits = 3, justify = "left")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat("                                   (Latency)                     \n")
  cat(paste(rep("-",90),collapse = ""),"\n")
  print.table(as.matrix(tablatencyout), digits = 3, justify = "left")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat(paste0("'Real' elapsed time: ",x$timer, " seconds.\n"))

}
