#' Fit a mixture cure model with Laplacian-P-splines.
#'
#' This routine fits a mixture cure model with a logistic link function for
#' the incidence part and a flexible Cox proportional hazards model for the
#' latency part where the baseline survival is approximated with penalized
#' B-splines. Laplace approximations are used to approximate the conditional
#' posterior distribution of the latent vector. A robust prior specification is
#' imposed on the roughness penalty parameter following Jullion and Lambert
#' (2007). The roughness penalty parameter is optimized and a maximum a
#' posteriori estimate is returned. The routine computes point estimates and
#' credible intervals for the latent parameters.
#'
#' @param formula A model formula of the form \code{Surv(tobs,delta)~
#' inci()+late()}.
#' @param data A data frame.
#' @param K The number of B-spline coefficients.
#' @param penorder The order of the penalty.
#' @param stepsize The stepsize taken to maximize the log posterior penalty.
#' @param deltaprior The parameters of the Gamma prior for the dispersion
#'  parameter.
#' @param v0 Initial parameter value for finding MAP of penalty parameter.
#' @param checkPD Should checks for positive definiteness be made? Default is TRUE.
#'
#' @return An object of class \code{lpsmc}.
#'
#' @seealso \link{lpsmc.object}
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @examples
#' ### Real data application ECOG e1684 clinical trial
#' data("ecog1684")
#' formula <- Surv(tobs,delta) ~ inci(SEX + TRT + AGE) + late(SEX + TRT + AGE)
#' fite1684 <- lpsmc(formula = formula, data = ecog1684)
#' fite1684
#'
#' ### Application on breast cancer data
#' rm(list=ls())
#' data("breastcancer")
#' formula <- Surv(tobs, delta) ~ inci(AGE + ER) + late(AGE + ER)
#' fitcancer <- lpsmc(formula = formula, data = breastcancer, K = 20)
#' fitcancer
#'
#' @references Jullion, A. and Lambert, P. (2007). Robust specification of the
#'  roughness penalty prior distribution in spatially adaptive Bayesian
#'  P-splines models. \emph{Computational Statistical & Data Analysis}
#'  \strong{51} (5), 2542-2558.\url{https://doi.org/10.1016/j.csda.2006.09.027}
#'
#' @export

lpsmc <- function(formula, data, K = 15, penorder = 3, stepsize = 0.2,
                  deltaprior = 1e-04, v0 = 15, checkPD = TRUE){

  #--- Start clock
  tic <- proc.time()

  #--- Extracting dimensions
  if(!inherits(formula, "formula"))
    stop("Incorrect model formula")
  # Reordering if needed
  if (as.character(stats::terms(formula)[[3]][[2]])[1] == "late" &&
      as.character(stats::terms(formula)[[3]][[3]])[1] == "inci") {
    formula <- stats::as.formula(paste(
      paste0("Surv(", all.vars(formula)[1], ",", all.vars(formula)[2], ")~"),
      attr(stats::terms(formula), "term.labels")[2],
      "+",
      attr(stats::terms(formula), "term.labels")[1],
      sep = ""))
  }
  varnames <- all.vars(formula, unique = FALSE)
  collect <- parse(text = paste0("list(", paste(varnames, collapse = ","), ")"),
                   keep.source = FALSE)

  if (missing(data)) {
    mff <- data.frame(matrix(unlist(eval(collect)),
                             ncol = length(varnames), byrow = FALSE))
  } else{
    mff <- data.frame(matrix(unlist(eval(collect, envir = data)),
                             ncol = length(varnames), byrow = FALSE))
  }
  colnames(mff) <- varnames
  incidstruct   <- attr(stats::terms(formula), "term.labels")[1]
  ncovar.inci   <- sum(strsplit(incidstruct, "+")[[1]] == "+") + 1
  latestruct    <- attr(stats::terms(formula), "term.labels")[2]
  ncovar.late   <- sum(strsplit(latestruct, "+")[[1]] == "+") + 1

  ftime <- mff[, 1]   # Event times
  event <- mff[, 2]   # Event indicator
  tu <- max(ftime)    # Upper bound of follow-up
  n <- nrow(mff)
  X <- as.matrix(cbind(rep(1, n), mff[, (3:(2 + ncovar.inci))]))
  colnames(X) <- c("(Intercept)", colnames(mff[(3:(2 + ncovar.inci))]))
  if(any(is.infinite(X)))
    stop("Data contains infinite covariates")
  if(ncol(X) == 0)
    stop("The model has no covariates for the incidence part")
  Z <- as.matrix(mff[, (2 + ncovar.inci + 1) : ncol(mff)], ncol = ncovar.late)
  if(any(is.infinite(Z)))
    stop("Data contains infinite covariates")
  if(ncol(Z) == 0)
    stop("The model has no covariates for the latency part")
  colnames(Z) <- colnames(mff[(2 + ncovar.inci + 1) : ncol(mff)])
  if (!is.vector(K, mode = "numeric") || length(K) > 1 || is.na(K))
    stop("K must be a numeric scalar")
  if (K < 10)
    stop("K must be at least 10")
  p <- ncol(X)
  q <- ncol(Z)
  dimlat <- K + p + q
  penorder <- floor(penorder)
  if(penorder < 2 || penorder > 3)
    stop("Penalty order must be either 2 or 3")

  #--- Cubic B-spline basis
  Bt <- Rcpp_cubicBspline(ftime, lower = 0, upper = tu, K = K)

  #-- Bspline basis evaluated on grid to approximate S0
  J <- 300                                                   # Number of bins
  phij <- seq(0, tu, length = J + 1)                         # Grid points
  Deltaj <- phij[2] - phij[1]                                # Interval width
  sj <- phij[1:J] + Deltaj * 0.5                             # Interval midpoint
  Bsj <- Rcpp_cubicBspline(sj, lower = 0, upper = tu, K = K) # B-spline basis

  #--- Prior precision and penalty matrix
  D <- diag(K)                       # Diagonal matrix
  for (k in 1:penorder) D <- diff(D) # Difference matrix of dimension K-r by K
  P <- t(D) %*% D                    # Penalty matrix
  P <- P + diag(1e-06,K)             # Diagonal perturbation to make P full rank
  prec_betagamma <- 1e-06            # Prior precision for the regression coeffs.
  a_delta <- deltaprior              # Prior for delta
  b_delta <- a_delta                 # Prior for delta
  nu <- 3                            # Prior for lambda

  # Precision matrix
  Qv <- function(v) {
    Qmat <- matrix(0, nrow = dimlat, ncol = dimlat)
    Qmat[1:K, 1:K] <- exp(v) * P
    Qmat[(K+1):dimlat, (K+1):dimlat] <- diag(prec_betagamma, p + q)
    return(Qmat)
  }

  BKL <- list()
  for (k in 1:K) {
    BKL[[k]] <- Bsj * matrix(rep(Bsj[, k], K), ncol = K, byrow = FALSE)
  }

  cumult <- function(t, theta, k, l) {
    bin_index <- as.integer(t / Deltaj) + 1
    bin_index[which(bin_index == (J + 1))] <- J
    h0sj <- exp(colSums(theta * t(Bsj)))

    if (k == 0 && l == 0) {
      output <- (cumsum(h0sj)[bin_index]) * Deltaj
    } else if (k == 1 && l == 0) {
      h0mat <- matrix(rep(h0sj, K), ncol = K, byrow = FALSE)
      output <- (apply(h0mat * Bsj, 2, cumsum)[bin_index,]) * Deltaj
    } else if (k == 1 && l > 0) {
      h0mat <- matrix(rep(h0sj, K), ncol = K, byrow = FALSE)
      output <- (apply(h0mat * BKL[[l]], 2, cumsum)[bin_index,]) * Deltaj
    }
    return(output)
  }

  #--- Incidence function (logistic)
  px <- function(betalat, x) as.numeric((1 + exp(-(x %*% betalat))) ^ (-1))

  #--- Population survival function
  Spop <- function(latent){
    thetalat <- latent[1:K]
    betalat <- latent[(K + 1):(K + p)]
    gammalat <- latent[(K + p + 1):dimlat]
    pX <- px(betalat,X)
    Zg <- as.numeric(Z%*%gammalat)
    Su <- exp(-exp(Zg)*cumult(ftime,thetalat,k=0,l=0))
    Spop_value <- 1-pX+pX*Su
    return(Spop_value)
  }

  #--- Log-likelihood function
  loglik <- function(latent){
    thetalat <- latent[1:K]
    betalat <- latent[(K + 1):(K + p)]
    gammalat <- latent[(K + p + 1):dimlat]
    pX <- px(betalat, X)
    Zg <- as.numeric(Z %*% gammalat)
    Btheta <- as.numeric(Bt %*% thetalat)
    cumul <- cumult(ftime, thetalat, k = 0, l = 0)
    Su <- exp(-exp(Zg) * cumul)
    loglikelihood <- sum(event * (log(pX) + Zg + Btheta - exp(Zg) * cumul) +
                           (1 - event) * log(1 - pX + pX * Su))
    return(loglikelihood)
  }

  ## Gradient
  Dloglik <- function(latent) {
    gradloglik <- c()
    thetalat <- latent[1:K]
    betalat <- latent[(K + 1):(K + p)]
    gammalat <- latent[(K + p + 1):dimlat]
    pX <- px(betalat,X)
    Zg <- as.numeric(Z %*% gammalat)
    Sdelta_ratio <- (1 - event) / Spop(latent)

    # Compute preliminary quantities
    omega_oi  <- cumult(ftime, thetalat, k = 0, l = 0)
    omega_oik <- cumult(ftime, thetalat, k = 1, l = 0)

    # Partial derivative wrt spline coefficients
    gradloglik[1:K] <- colSums(event * (Bt - exp(Zg) * omega_oik) -
                                 Sdelta_ratio * pX *
                                 exp(Zg - exp(Zg) * omega_oi) * omega_oik)

    # Partial derivative wrt beta coefficients
    gradloglik[(K + 1):(K + p)] <- colSums((event * (1 - pX) + Sdelta_ratio *
                        pX * (1 - pX) * (exp(-exp(Zg) * omega_oi) - 1)) * X)

    # Partial derivative wrt gamma coefficients
    gradloglik[(K + p + 1):dimlat] <- colSums((event * (1 - exp(Zg) *
                    omega_oi) - Sdelta_ratio * pX *exp(Zg - exp(Zg) *
                             omega_oi) * omega_oi) * Z)

    return(gradloglik)

  }

  ## Hessian
  D2loglik <- function(latent){
    thetalat <- latent[1:K]
    betalat <- latent[(K + 1):(K + p)]
    gammalat <- latent[(K + p + 1):dimlat]
    pX <- px(betalat, X)
    dpX <- pX * (1 - pX) * X
    Zg <- as.numeric(Z %*% gammalat)
    Sp <- Spop(latent)

    # Compute preliminary quantities
    omega_oi <- cumult(ftime, thetalat, k = 0, l = 0)
    omega_oik <- cumult(ftime, thetalat, k = 1, l = 0)
    ff <- exp(Zg - exp(Zg) * omega_oi) * omega_oik
    dSp_beta <- pX * (1 - pX) * (exp(-exp(Zg) * omega_oi) - 1) * X
    dSp_gamma <- (-pX * exp(Zg - exp(Zg) * omega_oi) * omega_oi) * Z
    ftilde <- pX * (1 - pX) * (exp(-exp(Zg) * omega_oi) - 1)
    dftilde_beta <- (pX * (1 - pX) * (1 - 2 * pX) *
                       (exp(-exp(Zg) * omega_oi) - 1)) * X
    fbreve <- (exp(-exp(Zg) * omega_oi) - 1)
    dfbreve_gamma <- (-exp(Zg - exp(Zg) * omega_oi) * omega_oi) * Z
    fddot <- exp(Zg - exp(Zg) * omega_oi)
    dfddot_gamma <- fddot * (1 - exp(Zg) * omega_oi) * Z

    # Block 11
    Block11 <- matrix(0, nrow = K, ncol = K)
    for (l in 1:K) {
      omega_oikl <- cumult(ftime, thetalat, k = 1, l = l)
      df <- exp(Zg - exp(Zg) * omega_oi) * omega_oikl - ff *
        exp(Zg) * omega_oik[, l]
      dSp <- (-pX * exp(Zg - exp(Zg) * omega_oi) * omega_oik[, l])

      Block11[l, ] <- colSums((-event * exp(Zg) * omega_oikl) -
                                (1 - event) * pX * (Sp ^ (-2)) *
                                (df * Sp - ff * dSp))
    }

    # Block 12
    Block12 <- t(ff) %*% (-(1 - event) * (dpX * Sp - pX * dSp_beta) * Sp ^
                            (-2))

    # Block 13
    Block13 <- t(omega_oik) %*% (-event * exp(Zg) * Z) -
      t(ff) %*% ((((1 - exp(Zg) * omega_oi) * Z * Sp) -
        ((-pX * exp(Zg - exp(Zg) * omega_oi) * omega_oi) * Z)) *
        ((1 - event) * pX * Sp ^ (-2)))

    # Block 22
    Block22 <- t(-event * pX * (1 - pX) * X) %*% X +
      (t((1 - event) * Sp ^ (-2) * (dftilde_beta * Sp - ftilde * dSp_beta)) %*%
         X)

    # Block 23
    Block23 <- t(X) %*% ((1 - event) * pX * (1 - pX) * Sp ^ (-2) *
                           (dfbreve_gamma * Sp - fbreve * dSp_gamma))

    # Block 33
    Block33 <- t(-event * exp(Zg) * omega_oi * Z) %*% Z -
      t(((1 - event) * pX * omega_oi * Sp ^ (-2)) *
          (dfddot_gamma * Sp - fddot * dSp_gamma)) %*% Z

    # Construction of Hessian matrix
    Hess <- matrix(0, nrow = dimlat, ncol = dimlat)
    Hess[(1:K),(1:K)]                       <- Block11
    Hess[(1:K),(K+1):(K+p)]                 <- Block12
    Hess[(K+1):(K+p),(1:K)]                 <- t(Block12)
    Hess[(1:K),(K+p+1):dimlat]              <- Block13
    Hess[(K+p+1):dimlat,(1:K)]              <- t(Block13)
    Hess[(K+1):(K+p),(K+1):(K+p)]           <- Block22
    Hess[(K+1):(K+p),(K+p+1):(dimlat)]      <- Block23
    Hess[(K+p+1):(dimlat),(K+1):(K+p)]      <- t(Block23)
    Hess[(K+p+1):(dimlat),(K+p+1):(dimlat)] <- Block33

    return(Hess)
  }

  ## Function to correct for positive definiteness if necessary
  PDcorrect <- function(x, eigentol = 1e-07) {
    eigvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    checkpositive <- !any(eigvalues < eigentol)
    if (checkpositive == FALSE) {
      correctPD <- x + diag(abs(min(eigvalues)) + 1e-04, ncol(x))
      note = "Matrix has been adjusted to PD"
      isPD = 0
    } else{
      correctPD <- x
      note = "Matrix is already PD"
      isPD = 1
    }
    outlist <- list(correctPD = correctPD, isPD = isPD, message = note)
    return(outlist)
  }

  ## Function to return full covariance matrix of dimension dimlat x dimlat
  Cov_full <- function(Cmatrix){
    Covfull <- matrix(0, nrow = dimlat, ncol = dimlat)
    Covfull[1:(K - 1), 1:(K - 1)] <- Cmatrix[1:(K - 1), 1:(K - 1)]
    Covfull[1:(K - 1), (K + 1):dimlat] <- Cmatrix[1:(K - 1), (K:(dimlat-1))]
    Covfull[(K + 1):dimlat, 1:(K - 1)] <- Cmatrix[K:(dimlat - 1), 1:(K - 1)]
    Covfull[(K+1):dimlat,(K+1):dimlat] <- Cmatrix[K:(dimlat-1),K:(dimlat-1)]
    return(Covfull)
  }

  #--- Posterior of (log)penalty parameter v=log(lambda)

  #--- log p(v|D) function
  logpv <- function(v, lat0){
    if(checkPD == FALSE){
      LL <- Rcpp_Laplace2(lat0 = lat0, v = v, K = K, Dloglik, D2loglik, Qv)
    } else{
      LL <- Rcpp_Laplace(lat0 = lat0, v = v, K = K, PDcorrect,
                       Dloglik, D2loglik, Qv)
    }
    latstar_cc <- LL$latstar                    # Mean of Laplace
    Covstar_c <- LL$Covstar                     # Covariance matrix of Laplace
    logdetCovstar_c <- Re(LL$logdetCovstarc)
    Q <- Qv(v)

    logpv_value <- 0.5 * logdetCovstar_c + loglik(latstar_cc) -
      0.5 * sum((latstar_cc * Q) %*% latstar_cc) +
      0.5 * (K + nu) * v - (0.5 * nu + a_delta) *
      log(0.5 * nu * exp(v) + b_delta)

    outlist <- list(value = logpv_value, latstar = latstar_cc)
    return(outlist)

  }

  #--- Exploration with golden search
  find_vmap <- function(){
    v0 <- v0
    vv <- c()
    lpvv <- c()
    vv[1] <- v0
    logpvinit <- logpv(v0, lat0 = rep(0, dimlat))
    lpvv[1] <- logpvinit$value
    lat0 <- logpvinit$latstar
    signdir <- 3
    signdirvec <- c()
    m <- 2

    while(signdir > 0){
      vv[m] <- vv[m - 1] - stepsize
      logpvm <- logpv(vv[m], lat0 = lat0)
      lpvv[m] <- logpvm$value
      lat0 <- logpvm$latstar
      signdir <- lpvv[m] - lpvv[m - 1]
      m <- m + 1
    }

    vmap <- (vv[m - 1] + vv[m - 2]) * 0.5
    outlist <- list(vmap = vmap, latstar = lat0)
    return(outlist)
  }

  vfind <- find_vmap()
  vstar <- vfind$vmap
  lat0  <- vfind$latstar
  if(checkPD == FALSE){
    LL <- Rcpp_Laplace2(lat0 = lat0, v = vstar, K = K, Dloglik, D2loglik, Qv)
  } else{
    LL <- Rcpp_Laplace(lat0 = lat0, v = vstar, K = K, PDcorrect,
                     Dloglik, D2loglik, Qv)
  }
  lathat <- LL$latstar

  #--- Estimated coefficients
  thetahat <- lathat[1:K]                  # Estimated latent vector
  betahat  <- lathat[(K + 1):(K + p)]      # Estimated betas (incidence)
  gammahat <- lathat[(K + p + 1):dimlat]   # Estimated gammas (latency)
  Covhat <- Cov_full(LL$Covstar)           # Covariance matrix dimlat x dimlat

  logpv2 <- function(v){
    if(checkPD == FALSE){
      LL <- Rcpp_Laplace2(lat0 = lathat, v = v, K = K, Dloglik, D2loglik, Qv)
    } else{
      LL <- Rcpp_Laplace(lat0 = lathat, v = v, K = K, PDcorrect,
                       Dloglik, D2loglik, Qv)
    }
    latstar_cc <- LL$latstar                    # Mean of Laplace
    Covstar_c <- LL$Covstar                     # Covariance matrix of Laplace
    logdetCovstar_c <- Re(LL$logdetCovstarc)
    Q <- Qv(v)

    logpv_value <- 0.5 * logdetCovstar_c + loglik(latstar_cc) -
      0.5 * sum((latstar_cc * Q) %*% latstar_cc) +
      0.5 * (K + nu) * v - (0.5 * nu + a_delta) *
      log(0.5 * nu * exp(v) + b_delta)

    return(logpv_value)

  }

  #---- Credible intervals of regression coefficients

  CI <- function(alpha){
    CImat <- matrix(0, nrow = (p + q), ncol = 2)
    zq <- stats::qnorm(p = .5 * alpha, lower.tail = F)
    sdcoeff <- sqrt(diag(Covhat)[(K + 1):dimlat])
    CImat[, 1] <- c(betahat, gammahat) - zq * sdcoeff
    CImat[, 2] <- c(betahat, gammahat) + zq * sdcoeff
    colnames(CImat) <- c(paste("lower.", 1 - alpha), paste("upper.", 1 - alpha))
    return(CImat)
  }

  CI90 <- CI(0.10) # 90% credible intervals for reg. coeffs.
  CI95 <- CI(0.05) # 95% credible intervals for reg. coeffs.

  #---- Credible interval for incidence p(x)

  CIp <- function(x, alpha){
    gbhat     <- log(log(as.numeric(1 + exp(-(x %*% betahat)))))
    Sigmabhat <- Covhat[(K + 1):(K + p), (K + 1):(K + p)]
    gradbhat  <- (-1) * ((1 - px(betahat, x))/
                           log(as.numeric(1 + exp(-(x %*% betahat))))) * x
    qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
    postsd    <- sqrt(as.numeric(t(gradbhat) %*% Sigmabhat %*% gradbhat))
    CIp_alpha <- c(gbhat - qz_alpha * postsd, gbhat + qz_alpha * postsd)
    CIp_original <- rev(exp(-exp(CIp_alpha)))
    return(CIp_original)
  }

  #---- Credible interval for cure rate 1-p(x)

  CIcure <- function(x, alpha){
    gbhat     <- log(log(as.numeric(1 + exp(x %*% betahat))))
    Sigmabhat <- Covhat[(K + 1):(K + p), (K + 1):(K + p)]
    gradbhat  <- (px(betahat, x)/log(as.numeric(1 + exp(x %*% betahat)))) * x
    qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
    postsd    <- sqrt(as.numeric(t(gradbhat) %*% Sigmabhat %*% gradbhat))
    CIcure_alpha <- c(gbhat - qz_alpha * postsd, gbhat + qz_alpha * postsd)
    CIcure_original <- rev(exp(-exp(CIcure_alpha)))
    return(CIcure_original)
  }

  if((colnames(X)[2] == "x1")){
    xmean <- c(1, 0, 0.5)
  } else {
   xmean <- as.numeric(apply(X,2,"mean"))
  }

  # Compute CI for incidence p(x)
  CI_incidence90  <- CIp(xmean, 0.10)
  CI_incidence95  <- CIp(xmean, 0.05)

  # Compute CI for cure rate 1-p(x)
  CI_cure90 <- CIcure(xmean, 0.10)
  CI_cure95 <- CIcure(xmean, 0.05)

  toc <- round(proc.time()-tic,3)[3]

  # List to output

  outlist <- list(
    X = X,
    Z = Z,
    latenthat = lathat,
    thetahat = thetahat,
    betahat = betahat,
    gammahat = gammahat,
    vhat = vstar,
    tu = tu,
    ftime = ftime,
    penorder = penorder,
    K = K,
    p = p,
    q = q,
    xmean = xmean,
    px = px,
    logpv2 = logpv2,
    cumult = cumult,
    CI90 = CI90,
    CI95 = CI95,
    CIcure90 = CI_cure90,
    CIcure95 = CI_cure95,
    CIincid90 = CI_incidence90,
    CIincid95 = CI_incidence95,
    Covhat = Covhat,
    timer = toc
  )

  attr(outlist,"class") <- "lpsmc"
  outlist
}

