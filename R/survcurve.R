#' Survival curves for the mixture cure model.
#'
#' @description
#' This routine creates objects related to the baseline survival
#' curve and the survival curve of the uncured subjects in a mixture
#' cure model fitted with \code{lpsmc}.
#'
#' @param x An object of type \code{lpsmc}.
#' @param type Either baseline or uncured.
#' @param covarprofile The matrix of covariate profiles.
#' @param cred.int The level for building the credible intervals.
#' @param themetype The theme of the plot either "classic", "gray","light"
#'  or "dark".
#' @param outlisted Should a list be produced as an output?
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @examples
#' ## Plot survival curves of uncured subjects for AGE=(30,50) and ER=1
#' rm(list=ls())
#' data("breastcancer")
#' formula <- Surv(tobs, delta) ~ inci(AGE + ER) + late(AGE + ER)
#' fitcancer <- lpsmc(formula = formula, data = breastcancer, K = 20)
#' profilematrix <- matrix(c(30, 1, 50, 1), nrow = 2, byrow = TRUE)
#' survcurve(fitcancer, type = "uncured", covarprofile = profilematrix)
#'
#' @export

survcurve <- function(x, type=c("baseline","uncured"), covarprofile=NULL,
                      cred.int = 0.95,
                      themetype = c("classic","gray","light","dark"),
                      outlisted = FALSE){

  themetype <- match.arg(themetype)
  if(themetype == "classic"){
    themeval<- eval(parse(text="ggplot2::theme_classic()"))
  } else if(themetype == "gray"){
    themeval <- eval(parse(text="ggplot2::theme_gray()"))
  } else if (themetype == "light"){
    themeval <- eval(parse(text="ggplot2::theme_light()"))
  } else if (themetype == "dark"){
    themeval <- eval(parse(text="ggplot2::theme_dark()"))
  }

  # Extract variables
  K <- x$K
  thetahat <- x$thetahat
  gammahat <- x$gammahat
  tu <- x$tu

  Sdom <- seq(0, tu, length = 100) # Domain on which to plot functions

  basehaz <- function(t){
    val <- exp(as.numeric(Rcpp_cubicBspline(t, lower = 0,
                    upper = tu, K = K) %*%  thetahat))
    return(val)
  }

  basesurv <- function(t){
    val <- exp(-stats::integrate(basehaz, lower = 0, upper = t)$value)
    return(val)
  }

  S0hat <- sapply(Sdom,basesurv)
  plot_type <- match.arg(type)

  if(plot_type=="baseline"){

  # Function for CI for S0 at t
  S0CIfun <- function(t, alpha){
      cum_tt <- x$cumult(t = t, theta = thetahat, k = 0, l = 0)
       gstar  <- log(cum_tt)
      grad_g <- x$cumult(t = t, theta = thetahat, k = 1, l = 0)[-K] *
        (1 / cum_tt)
      Sig_thetahat <- x$Covhat[1:(K-1), 1:(K-1)]
      qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
      post_sd <- sqrt(as.numeric(t(grad_g) %*% Sig_thetahat %*% grad_g))
      CI_g <- c(gstar - qz_alpha * post_sd, gstar + qz_alpha * post_sd)
      CIS <- rev(exp(-exp(CI_g)))
      return(CIS)
  }

  S0CI <- sapply(Sdom, S0CIfun, alpha = 1 - cred.int)
  S0CIlow <- S0CI[1,]
  S0CIup <- S0CI[2,]

  S0dat <- data.frame(Sdom, S0hat, S0CIlow, S0CIup)
  S0dat <- S0dat[1:utils::head(which(S0dat$S0hat < 0.001), 1), ]
  S0dat[1,4] <- 1

  # baseline survival plot
  S0plot <-  ggplot2::ggplot(data = S0dat, ggplot2::aes(x=Sdom, y=S0hat)) +
                   ggplot2::geom_line(size=1.1, colour="red") +
                     themeval +
                       ggplot2::xlab("t") +
                         ggplot2::ylab(expression(S[0](t))) +
              ggplot2::theme(axis.title.x = ggplot2::element_text(size = 14),
                  axis.title.y = ggplot2::element_text(size = 14),
                  axis.text.x = ggplot2::element_text(size=12),
                  axis.text.y = ggplot2::element_text(size=12)) +
              ggplot2::geom_ribbon(ggplot2::aes(ymin=S0CIlow, ymax=S0CIup),
                          alpha = 0.15,fill="firebrick2")

  outlist <- list(tt = Sdom,
                  S0fit = S0hat,
                  S0plot = S0plot,
                  S0CIup = S0CIup,
                  S0CIlow = S0CIlow)

  S0plot # display plot


  } else if(plot_type=="uncured"){

  # survival of uncured subject given covariate profile

  SuCIfun <- function(t, alpha, covarprofile){
    cum_tt <- x$cumult(t = t, theta = thetahat, k = 0, l = 0)
    gstar  <- as.numeric(covarprofile %*% gammahat) + log(cum_tt)
    grad_g <- c(x$cumult(t = t, theta = thetahat, k = 1, l = 0)[-K] *
                  (1 / cum_tt), covarprofile)
    p <- x$p
    q <- x$q
    dimlat <- K+p+q
    SSighat <- matrix(0, nrow = (K - 1) + q, ncol = (K - 1) + q)
    SSighat[1:(K-1),1:(K-1)] <- x$Covhat[1:(K-1), 1:(K-1)]
    SSighat[1:(K-1),K:(K+1)] <- x$Covhat[1:(K - 1), (K + p + 1):dimlat]
    SSighat[K:(K+1),1:(K-1)] <- t(x$Covhat[1:(K - 1), (K + p + 1):dimlat])
    SSighat[K:(K+1),K:(K+1)] <- x$Covhat[(K + p + 1):dimlat, (K + p + 1):dimlat]
    qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
    post_sd <- sqrt(as.numeric(t(grad_g) %*% SSighat %*% grad_g))
    CI_g <- c(gstar - qz_alpha * post_sd, gstar + qz_alpha * post_sd)
    CIS <- rev(exp(-exp(CI_g)))
    return(CIS)
  }

  covarprofile1 <- as.numeric(covarprofile[1, ])

  Suhat <- S0hat ^ (exp(as.numeric(covarprofile1 %*% gammahat)))
  SuCI <- sapply(Sdom, SuCIfun, alpha = 1 - cred.int,
                 covarprofile = covarprofile1)
  SuCIlow <- SuCI[1,]
  SuCIup <- SuCI[2,]
  xlimup <- Sdom[utils::head(which(Suhat < 1e-6),1)]

  plot(Sdom, Suhat, type = "l", col = "black", ylab="Survival in latency part",
       xlab="Time (in days)", xlim = c(0, xlimup))
  graphics::polygon(c(Sdom,rev(Sdom)),c(SuCIlow,rev(SuCIup)),col="gray85",
                    border = "NA")
  graphics::lines(Sdom, Suhat, type="l", lwd = 2)

  if(nrow(covarprofile) > 1){
  covarprofile2 <- as.numeric(covarprofile[2, ])

  Suhat2 <- S0hat ^ (exp(as.numeric(covarprofile2 %*% gammahat)))
  SuCI <- sapply(Sdom, SuCIfun, alpha = 1 - cred.int,
                 covarprofile=covarprofile2)
  SuCIlow <- SuCI[1,]
  SuCIup <- SuCI[2,]

  graphics::lines(Sdom, Suhat2, type="l", lty=2, lwd=2)
  graphics::polygon(c(Sdom,rev(Sdom)),c(SuCIlow,rev(SuCIup)),col="NA",
          border = "black", lty=2, lwd = 1)
  }

  ctitle1 <- paste0(colnames(x$X)[2:x$p], "=", covarprofile1)
  ctitle1 <- paste(ctitle1[1], ctitle1[2], sep="; ")
  ctitle2 <- paste0(colnames(x$X)[2:x$p], "=", covarprofile2)
  ctitle2 <- paste(ctitle2[1], ctitle2[2], sep="; ")


  graphics::legend("topright",c(ctitle1,ctitle2), lty=c(1,2),
                   col = c("black","black"),
                    bty="n", lwd = c(2,2))

  }

  if(outlisted == TRUE){
    return(outlist)
  }

}
