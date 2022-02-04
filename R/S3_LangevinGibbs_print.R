#' Print a LangevinGibbs object.
#'
#' @description Print method for a \code{LangevinGibbs} object.
#'
#' @param x An object of class \code{LangevinGibbs}.
#' @param ... Further arguments to be passed to print routine.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @export


print.LangevinGibbs <- function(x,...){

  tabcolnames <- c("Estimate","sd","CI90%.low","CI90.up%",
                   "CI95%.low","CI95%.up")
  K <- x$K

  # Table for the incidence part
  tabincidrow <- x$p
  tabincidout <- as.data.frame(matrix(0, nrow = tabincidrow, ncol = 6))
  rownames(tabincidout) <- colnames(x$X)
  colnames(tabincidout) <- tabcolnames
  tabincidout[,1] <- x$betahat
  tabincidout[,2] <- x$sdpost[1:x$p]
  tabincidout[,3] <- x$CI90[,1][1:tabincidrow]
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
  tablatencyout[,2] <- x$sdpost[(x$p+1):length(x$sdpost)]
  tablatencyout[,3] <- x$CI90[,1][(x$p+1):length(x$sdpost)]
  tablatencyout[,4] <- x$CI90[,2][(x$p+1):length(x$sdpost)]
  tablatencyout[,5] <- x$CI95[,1][(x$p+1):length(x$sdpost)]
  tablatencyout[,6] <- x$CI95[,2][(x$p+1):length(x$sdpost)]
  iscolnum <- sapply(tablatencyout, is.numeric)
  tablatencyout[iscolnum] <- lapply(tablatencyout[iscolnum], round, 3)

  # Print output table

  cat("Fitting mixture cure model with Langevin-Gibbs sampler \n")
  cat(paste(rep("-",55),collapse = ""),"\n")
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
  cat(paste0("'Real' elapsed time: ",x$timer, " seconds\n"))
  cat(paste0("MCMC chain length: ",x$mcmcsample, "\n"))
  cat(paste0("Burn-in length: ",x$burnin,"\n"))
  cat(paste0("MCMC acceptance rate: ",x$acceptrate,"%.\n"))

}
