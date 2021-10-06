#' Estimated cure proportion
#'
#' @description
#' Computes the estimated cure proportion based on a mixture cure model fit
#' with \code{lpsmc}. Point estimates and approximate 90% and 95% credible
#' intervals are shown.
#'
#' @param x A lpsmc object.
#' @param covarprofile The covariate profile on which to compute the
#' cure proportion.
#'
#' @return A table with the estimated cure proportion.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @examples
#' ### Application on breast cancer data
#' rm(list=ls())
#' data("breastcancer")
#' formula <- Surv(tobs, delta) ~ inci(AGE + ER) + late(AGE + ER)
#' fitcancer <- lpsmc(formula = formula, data = breastcancer, K = 20)
#' covarprofile <- matrix(c(1, 30, 1, 1, 40, 0), nrow = 2 , byrow = TRUE)
#' fitcure <- curefit(fitcancer, covarprofile)
#' fitcure$estimcure
#'
#' @export


curefit <- function(x, covarprofile){

  betahat <- x$betahat
  phat <- x$px(betahat, covarprofile)
  p <- x$p
  K <- x$K
  nprofiles <- nrow(covarprofile)

  CIcure <- function(y, alpha){
    gbhat     <- log(log(as.numeric(1 + exp(y %*% betahat))))
    Sigmabhat <- x$Covhat[(K + 1):(K + p), (K + 1):(K + p)]
    gradbhat  <- (x$px(betahat, y)/log(as.numeric(1 + exp(y %*% betahat)))) * y
    qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
    postsd    <- sqrt(as.numeric(t(gradbhat) %*% Sigmabhat %*% gradbhat))
    CIcure_alpha <- c(gbhat - qz_alpha * postsd, gbhat + qz_alpha * postsd)
    CIcure_original <- rev(exp(-exp(CIcure_alpha)))
    return(CIcure_original)
  }

  CIcuremat <- matrix(0, nrow = nprofiles , ncol = p + 5)
  colnames(CIcuremat) <- c(as.character(colnames(x$X)),"1-p(x)",
                           "CI90.low","CI90.up","CI95.low","CI95.up")
  rownames(CIcuremat) <- paste0("x.profile", seq(nprofiles))
  for(j in 1:nprofiles){
   CIcuremat[j, 1:p] <- covarprofile[j, ]
   CIcuremat[j, p+1] <- 1-phat[j]
   CIcuremat[j, (p+2):(p+3)] <- CIcure(covarprofile[j,], 0.10)
   CIcuremat[j, (p+4):(p+5)] <- CIcure(covarprofile[j,], 0.05)
  }

  outlist <- list(estimcure = CIcuremat)

}
