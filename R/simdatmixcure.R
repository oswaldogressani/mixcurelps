#' Simulation of survival data to fit mixture cure models.
#'
#' @description
#' This routines simulates survival data with a cure fraction. The data is
#' simulated according to a mixture cure model with a logistic link for the
#' incidence part of the model. The incidence part includes two covariates
#' following a standard normal and a Bernoulli distribution with success
#' probability equal to 0.5. The latency part assumes a Weibull-Cox model with
#' two covariates (a standard normal variate and a Bernoulli variate with
#' success probability equal to 0.4).
#'
#' @param n Sample size.
#' @param wscale The positive scale parameter of the Weibull distribution used
#'  to generate the survival times of the uncured subjects.
#' @param wshape The positive shape parameter of the Weibull distribution used
#'  to generate the survival times of the uncured subjects.
#' @param setting The setting under which survival times will be generated. If
#'  \code{setting = 1}, the coefficients of the incidence part are
#'  \emph{beta0=0.70, beta1=-1.15 and beta2=0.95} and the coefficients of the
#'  latency part are \emph{gamma1=-0.10 and gamma2=0.25}. If
#'  \code{setting = 2}, the coefficients of the incidence part are
#'  \emph{beta0=1.25, beta1=-0.75 and beta2=0.45} and the coefficients of the
#'  latency part are \emph{gamma1=-0.10 and gamma2=0.20}.
#'
#' @return
#' An object of class \code{simixcure} containing different objects of the
#' simulated dataset. Details can be found by typing ?simixcure.object.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @seealso \link{simdatmixcure.object}
#'
#' @examples
#' ### Simulate a sample of size n=300 under Scenario 1.
#' set.seed(4408)
#' simdat <- simdatmixcure(n = 300, wshape = 1.45, wscale = 0.25, setting = 1)
#' plot(simdat) # Plot the baseline survival and Kaplan-Meier curve
#' simdat$info  # Print information on Cure and Censoring levels
#'
#' @export

simdatmixcure <- function(n, wscale, wshape, setting) {

  #--- Incidence part (logistic model)
  if(setting == 1){
    betas <- c(0.70,-1.15, 0.95)
    gammas <- c(-0.10, 0.25)
    censrate <- 0.16
    tau0 <- 8
  } else if (setting == 2){
    betas <- c(1.25, -0.75, 0.45)
    gammas <- c(-0.10, 0.20)
    censrate <- 0.05
    tau0 <- 8
  }
  X <- cbind(rep(1, n), stats::rnorm(n), stats::rbinom(n, 1, prob = 0.5))
  Xbetas <- as.numeric(X %*% betas)
  p <- (1 + exp(-Xbetas)) ^ (-1)

  #--- Generation of cure status
  B <- stats::rbinom(n, size = 1, prob = p) # Cure status B = 1 --> uncured

  #--- Generation of survival times for uncured (B=1) subject from Weibull CoxPH
  Z <- cbind(stats::rnorm(n), stats::rbinom(n, size = 1, prob = 0.4))

  ## Draw survival times from Weibull distribution
  weibshape <- wshape # Weibull shape parameter > 0
  weibscale <- wscale # Weibull scale parameter > 0
  S0 <- function(t) exp(-weibscale * t ^ (weibshape)) # True baseline survival
  h0 <- function(t) weibscale * weibshape * t ^ (weibshape - 1) # True h0
  U <- stats::runif(n) # Uniform draw
  Tlat <- as.numeric((-log(U) / (weibscale * exp(Z %*% gammas))) ^
                       (1 / weibshape))
  Tlat[Tlat > tau0] <- tau0    # Truncation of survival times
  Tlat[which(B == 0)] <- 20000 # Large survival time for cured subjects
  tobs <- Tlat

  #--- Censoring follows exponential distribution with rate lambda
  tau1 <- 11
  tup <- 11
  C <- stats::rexp(n, rate = censrate)
  C[C > tau1] <- tau1
  TgreatC <- which(Tlat > C)
  tobs[TgreatC] <- C[TgreatC]
  delta <- as.numeric(Tlat <= C)

  #--- Summary statistics on cure and censoring rates
  dataKM <- as.data.frame(cbind(tobs, delta))
  fitKM <- survival::survfit(survival::Surv(tobs, delta) ~ 1, data = dataKM)
  plateau <-
    fitKM$time[utils::tail(which((diff(fitKM$surv) < 0) == TRUE), 1) + 1]
  nobs_plateau <- sum(tobs > plateau)
  infomat <- matrix(0, nrow = 1, ncol = 6)
  colnames(infomat) <- c("n", "Cure level", "Cens.rate",
                         "Cens.level","% obs in plateau", "% cured in plateau")
  infomat[1, 1] <- n                                     # Sample size
  infomat[1, 2] <- round(sum(B == 0) / n * 100, 3)       # Cure rate
  infomat[1, 3] <- censrate
  infomat[1, 4] <- round(sum(1 - delta) / n * 100, 3)    # Censoring level
  infomat[1, 5] <- round(nobs_plateau / n * 100, 3)      # % of obs. in plateau
  infomat[1, 6] <- round(sum(B[which(tobs > plateau)] == 0) /
                           nobs_plateau * 100, 3)

  # Extract variables
  simdata <- data.frame(tobs, delta, X, Z)
  colnames(simdata) <- c("tobs","event","Intercept","x1","x2","z1","z2")

  outlist <- list(n = n,              # sample size
                  X = X,              # covariate matrix of incidence part
                  Z = Z,              # covariate matrix of latency part
                  betas = betas,      # coeffs. of incidence model
                  gammas = gammas,    # coeffs. of latency model
                  tobs = tobs,        # observed failure or censoring time
                  delta = delta,      # event indicator 1 --> event occurred
                  h0 = h0,            # true baseline hazard
                  S0 = S0,            # true baseline survival
                  tup = tup,          # Upper bound of follow-up
                  fitKM = fitKM,      # Kaplan-Meier fit
                  plateau = plateau,  # Plateau value
                  info = infomat,     # Summary statistics
                  setting = setting,  # The chosen setting
                  wshape = wshape,    # Shape parameter
                  wscale = wscale,    # Scale parameter
                  simdata = simdata)  # Simlated dataframe

   attr(outlist, "class") <- "simdatmixcure"
   outlist
}




















