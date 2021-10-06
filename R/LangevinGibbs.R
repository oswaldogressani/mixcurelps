#' Langevin-Gibbs sampler for inference in a mixture cure model
#'
#'
#' @description
#' This routine implements a Metropolis-Langevin-within-Gibbs sampler to draw
#' samples from the posterior distribution of a mixture cure model. A
#' Metropolis-adjusted Langevin algorithm is used to sample from the
#' conditional posterior distribution of the latent vector. MCMC samples for
#' the roughness penalty parameter and the dispersion parameter are obtained
#' via a Gibbs sampler.
#'
#' @param formula A model formula of the form \code{Surv(tobs,delta)~
#' inci()+late()}.
#' @param data A data frame.
#' @param K The number of B-spline coefficients.
#' @param penorder The order of the penalty associated to the B-spline
#'  coefficients.
#' @param deltaprior The parameters of the Gamma prior for the dispersion
#'  parameter.
#' @param mcmcsample The length of the MCMC chain.
#' @param burnin The length of the burnin.
#' @param mcmcseed The seed to be used (for reproducibility).
#' @param progbar Should a progress bar be shown? Default is yes.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @export



LangevinGibbs <- function(formula, data, K = 15, penorder = 3, deltaprior = 1e-04,
                         mcmcsample = 10000, burnin = 2000,
                         mcmcseed = NULL, progbar = c("yes","no")){

  tic <- proc.time()

  showprogbar <- match.arg(progbar)
  if(showprogbar == "yes"){
    programbar <- TRUE
  } else if (showprogbar == "no"){
    programbar <- FALSE
  } else{
    print("Error progbar should be either 'yes' or 'no'.")
  }

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
  incidstruct <- attr(stats::terms(formula), "term.labels")[1]
  ncovar.inci <- sum(strsplit(incidstruct, "+")[[1]] == "+") + 1
  latestruct <- attr(stats::terms(formula), "term.labels")[2]
  ncovar.late <- sum(strsplit(latestruct, "+")[[1]] == "+") + 1

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
  P <- t(D) %*% D                    # Penalty matrix of dimension K-1 by K-1
  P <- P + diag(1e-06,K)             # Diagonal perturbation to make P full rank
  prec_betagamma <- 1e-06            # Prior precision for the regression coeffs.
  a_delta <- deltaprior              # Prior for delta
  b_delta <- a_delta                 # Prior for delta
  nu <- 3                            # Parameter in lambda prior

  # Precision matrix
  Qv <- function(v) {
    Qmat <- matrix(0, nrow = dimlat, ncol = dimlat)
    Qmat[1:K, 1:K] <- exp(v) * P
    Qmat[(K+1):dimlat, (K+1):dimlat] <- diag(prec_betagamma, p + q)
    return(Qmat)
  }

  BKL <- list()
  for (k in 1:K) BKL[[k]] <- Bsj * matrix(rep(Bsj[, k], K),
                                          ncol = K, byrow = FALSE)

  cumult <- function(t, theta,k,l){
    bin_index <- as.integer(t / Deltaj) + 1
    bin_index[which(bin_index == (J + 1))] <- J
    h0sj <- exp(colSums(theta * t(Bsj)))

    if(k==0 && l==0){
      output <- (cumsum(h0sj)[bin_index]) * Deltaj
    } else if (k==1 && l==0){
      h0mat <- matrix(rep(h0sj, K), ncol = K, byrow = FALSE)
      output <- (apply(h0mat * Bsj, 2, cumsum)[bin_index, ]) * Deltaj
    } else if (k==1 && l>0){
      h0mat <- matrix(rep(h0sj, K), ncol = K, byrow = FALSE)
      output <- (apply(h0mat * BKL[[l]], 2, cumsum)[bin_index, ]) * Deltaj
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

  #--- Gradient of log-likelihood function
  Dloglik <- function(latent) {

    gradloglik <- c()
    thetalat <- latent[1:K]
    betalat <- latent[(K + 1):(K + p)]
    gammalat <- latent[(K + p + 1):dimlat]
    pX <- px(betalat,X)
    Zg <- as.numeric(Z %*% gammalat)
    Sdelta_ratio <- (1 - event) / Spop(latent)

    # Compute preliminary quantities
    omega_oi <- cumult(ftime, thetalat, k=0, l=0)
    omega_oik <- cumult(ftime, thetalat, k=1, l=0)

    # Partial derivative wrt spline coefficients
    gradloglik[1:K] <- colSums(event * (Bt - exp(Zg) * omega_oik) -
                                 Sdelta_ratio * pX *
                                 exp(Zg - exp(Zg) * omega_oi) * omega_oik)

    # Partial derivative wrt beta coefficients
    gradloglik[(K + 1):(K + p)] <- colSums((event * (1 - pX) + Sdelta_ratio * pX *
                                              (1 - pX) * (exp(-exp(Zg) * omega_oi) - 1)) * X)

    # Partial derivative wrt gamma coefficients
    gradloglik[(K + p + 1):dimlat] <- colSums((event * (1 - exp(Zg) * omega_oi) -
                                                 Sdelta_ratio * pX *exp(Zg - exp(Zg) *
                                                                          omega_oi) * omega_oi) * Z)

    return(gradloglik)

  }

  #---------------------- Metropolis-Langevin-within-Gibbs sampler

  # log posterior of conditional latent vector (given roughness penalty)
  logtar <- function(latent, lambda) {
    v <- log(lambda)
    Q <- Qv(v)
    val <-  loglik(latent) - 0.5 * sum((latent * Q) %*% latent)
    return(val)
  }

  # Gradient of logtar
  Dlogtar <- function(latent, lambda){
    v <- log(lambda)
    Q <- Qv(v)
    val <- as.numeric(Dloglik(latent) - (Q %*% latent))
    return(val)
  }

  # Extract covariance matrix of lpsmc fit
  fitlps <- lpsmc(formula = formula, data = data, K = K,
                  penorder = penorder, deltaprior = deltaprior)
  Sigprop <- fitlps$Covhat[-K,-K]
  Pblock <- matrix(0, nrow = dimlat, ncol = dimlat)
  Pblock[1:K, 1:K] <- P

  Langevin <- function(M){

    # Initial values for latent vector and hyperparameters
    lat0 <-  c(fitlps$thetahat, fitlps$betahat, fitlps$gammahat)
    lambda0 <- exp(fitlps$vhat)
    delta0 <- (nu * 0.5 + a_delta)/ (nu * 0.5 * lambda0 + b_delta)

    # Chain hosting
    thetachain <- matrix(0, nrow = M, ncol = K)
    betachain <- matrix(0, nrow = M, ncol = p)
    gammachain <- matrix(0, nrow =M, ncol = q)
    lambdachain<- c()
    deltachain <- c()
    naccept <- 0

    progbar <- progress::progress_bar$new(
      format = crayon::white$green("Langevin-Gibbs sampling [:elapsed :spin] [:bar] :percent"),
      total = M,
      clear = FALSE
    )

    # Variance parameter for Langevin dynamics
    Ldelta <- 0.25

    for (m in 1:M) {

      meanprop <- lat0[-K] + as.numeric(0.5 * Ldelta *
                                    (Sigprop %*% Dlogtar(lat0, lambda0)[-K]))
      Covprop <- Ldelta * Sigprop
      latnew <- MASS::mvrnorm(n = 1, mu = meanprop, Sigma = Covprop)
      latnew <- c(latnew[1:(K - 1)], 1, latnew[K:(dimlat - 1)])

      Gradprop <- Dlogtar(latnew,lambda0)
      Gradcurr <- Dlogtar(lat0,lambda0)

      logprop_ratio <- as.numeric((-0.5) * (t(Gradprop + Gradcurr)[-K]) %*%
        as.numeric((latnew[-K] - lat0[-K]) + (0.25 * Ldelta) * Sigprop %*%
                     (Gradprop - Gradcurr)[-K]))

      ldiff <- (logtar(latnew, lambda0) - logtar(lat0, lambda0)) +
                    logprop_ratio

      if (ldiff < 0) {
        u <- stats::runif(1)
        if (log(u) < ldiff) {
          thetachain[m,] <- latnew[1:K]
          betachain[m, ] <- latnew[(K+1):(K+p)]
          gammachain[m, ] <- latnew[(K+p+1):dimlat]
          lat0 <- latnew
          naccept <- naccept + 1
        } else{
          thetachain[m,] <- lat0[1:K]
          betachain[m, ] <- lat0[(K+1):(K+p)]
          gammachain[m,] <- lat0[(K+p+1):dimlat]
        }
      } else{
        thetachain[m,] <- latnew[1:K]
        betachain[m, ] <- latnew[(K+1):(K+p)]
        gammachain[m, ] <- latnew[(K+p+1):dimlat]
        lat0 <- latnew
        naccept <- naccept + 1
      }

      lambda_shape <- 0.5 * (K + nu)
      lambda_rate <- 0.5 * (sum((lat0 * Pblock) %*% lat0) + nu * delta0)
      lambda0 <- stats::rgamma(n = 1, shape = lambda_shape, rate = lambda_rate)
      lambdachain[m] <- lambda0

      # Draw dispersion parameter from Gamma distribution
      delta_shape <- 0.5 * nu + a_delta
      delta_rate <- 0.5 * nu * lambda0 + b_delta
      delta0 <- stats::rgamma(n = 1, shape = delta_shape, rate = delta_rate)
      deltachain[m] <- delta0

      # Automatic tuning of Langevin algorithm
      accept_prob <- min(c(1, exp(ldiff)))
      heval <- sqrt(Ldelta) + (1 / m) * (accept_prob - 0.57)

      hfun <- function(x) {
        epsil <- 1e-04
        Apar <- 10 ^ 4
        if (x < epsil) {
          val <- epsil
        } else if (x >= epsil && x <= Apar) {
          val <- x
        } else{
          val <- Apar
        }
        return(val)
      }

      Ldelta <- (hfun(heval)) ^ 2

      if(programbar == TRUE){
        progbar$tick()
      }
    }

    accept_rate <- round((naccept / M) * 100, 3)
    outlist <- list(accept_rate = accept_rate,
                    thetachain = thetachain,
                    betachain = betachain,
                    gammachain = gammachain,
                    lambdachain = lambdachain,
                    deltachain = deltachain,
                    M = M)
    return(outlist)

  }

  if(!is.null(mcmcseed)){
    set.seed(mcmcseed)
  }

  # Run Metropolis-Langevin-within-Gibbs sampler
  MCMCG <- Langevin(M = mcmcsample)

  # Extract MCMC output
  thetachain <- coda::as.mcmc(MCMCG$thetachain)
  betachain <- coda::as.mcmc(MCMCG$betachain)
  gammachain <- coda::as.mcmc(MCMCG$gammachain)
  lambdachain <- coda::as.mcmc(MCMCG$lambdachain)
  deltachain <- coda::as.mcmc(MCMCG$deltachain)
  acceptrate <- MCMCG$accept_rate

  # Geweke diagnostics
  Geweke <- matrix(
    c((coda::geweke.diag(thetachain[(burnin + 1):mcmcsample,])$z)[-K],
      coda::geweke.diag(betachain[(burnin + 1):mcmcsample,])$z,
      coda::geweke.diag(gammachain[(burnin + 1):mcmcsample,])$z,
      coda::geweke.diag(lambdachain[(burnin + 1):mcmcsample])$z,
      coda::geweke.diag(deltachain[(burnin + 1):mcmcsample])$z),
      ncol = 1)
  colnames(Geweke) <- "Geweke z-scores"
  rownames(Geweke) <- c(paste0("theta",seq_len(K-1)),
                        c("Intercept",paste0("beta",seq_len(p-1))),
                        paste0("gamma",seq_len(q)),
                        "lambda","delta")

  thetahat <- colMeans(thetachain[(burnin + 1):mcmcsample, ])
  betahat <- colMeans(betachain[(burnin + 1):mcmcsample, ])
  gammahat <- colMeans(matrix(gammachain[(burnin + 1):mcmcsample,], ncol = q))
  vhat <- mean(log(lambdachain[(burnin + 1):mcmcsample]))
  sdpost <- apply(cbind(betachain[(burnin + 1):mcmcsample, ],
                        gammachain[(burnin + 1):mcmcsample, ]), 2, "sd")

  # Quantile-based CI
  CI90 <- matrix(0, nrow = (p+q), ncol = 2)
  CI95 <- matrix(0, nrow = (p+q), ncol = 2)
  bgmat <- cbind(betachain[(burnin + 1):mcmcsample, ],
                 gammachain[(burnin + 1):mcmcsample, ])

  for (j in 1:(p + q)) {
    CI90[j,] <- stats::quantile(bgmat[,j], probs = c(0.05, 0.95))
    CI95[j,] <- stats::quantile(bgmat[,j], probs = c(0.025, 0.975))
  }

  #--- Recording time
  toc <- round(proc.time()-tic,3)[3]
  toc <- toc - fitlps$timer

  outlist <- list(
    X = X,
    Z = Z,
    thetahat = thetahat,
    thetachain = thetachain,
    betahat = betahat,
    betachain = betachain,
    gammahat = gammahat,
    gammachain = gammachain,
    vhat = vhat,
    lambdachain = lambdachain,
    deltachain = deltachain,
    tu = tu,
    ftime = ftime,
    penorder = penorder,
    K = K,
    p = p,
    q = q,
    px = px,
    cumult = cumult,
    CI90 = CI90,
    CI95 = CI95,
    sdpost = sdpost,
    acceptrate = acceptrate,
    Geweke = Geweke,
    timer = toc,
    burnin = burnin,
    mcmcsample = mcmcsample
  )

  attr(outlist,"class") <- "LangevinGibbs"
  outlist

}
