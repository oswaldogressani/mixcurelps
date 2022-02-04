#' Assessing the statistical performance of LangevinGibbs.
#'
#' @description
#' Simulates the LangevinGibbs routine.
#'
#' @param n Sample size.
#' @param K Number of B-spline basis functions.
#' @param scenario Either 1 or 2.
#' @param S Total number of replications.
#' @param chainlength The length of the MCMC chain.
#' @param burn Burnin length.
#' @param exactrep Exactly replicate results.
#' @param simnum A seed for reproducibility.
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#'
#' @export

simLangevin <- function(n = 300, K = 15, scenario = 1, S = 500,
                       chainlength = 1500, burn = 500, exactrep = FALSE,
                       simnum = NULL){


  if(exactrep == TRUE) {
    # Setting seeds (for reproducibility)
    if (scenario == 1 && n == 300) {
      set.seed(1986) # FIFA World cup in Mexico
    } else if (scenario == 2 && n == 300) {
      set.seed(1789) # France initiates revolution
    } else if (scenario == 1 && n == 600) {
      set.seed(1080) # Battle on the Elster
    } else if (scenario == 2 && n == 600) {
      set.seed(1055) # Death of Constantine IX
    }
  }

  # Prepare matrix to host esimates
  thetamat <- matrix(0, nrow = S, ncol = K) # Matrix to host estim. spline coef.
  betamat  <- matrix(0, nrow = S, ncol = 3) # Matrix to host estimated betas
  gammamat <- matrix(0, nrow = S, ncol = 2) # Matrix to host estimated gammas
  coverage90 <- matrix(0, nrow = S, ncol = 5) # Matrix for 90% coverage proba.
  coverage95 <- matrix(0, nrow = S, ncol = 5) # Matrix for 95% coverage proba.
  pvec <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
  S0cover90 <- matrix(0, nrow = S, ncol = length(pvec)) # For 90% coverage of S0
  S0cover95 <- matrix(0, nrow = S, ncol = length(pvec)) # For 95% coverage of S0
  Sucover90 <- matrix(0, nrow = S, ncol = length(pvec)) # For 90% coverage of Su
  Sucover95 <- matrix(0, nrow = S, ncol = length(pvec)) # For 95% coverage of Su
  S090CIwidth <- matrix(0, nrow = S, ncol = length(pvec))
  S095CIwidth <- matrix(0, nrow = S, ncol = length(pvec))
  Su90CIwidth <- matrix(0, nrow = S, ncol = length(pvec))
  Su95CIwidth <- matrix(0, nrow = S, ncol = length(pvec))
  Diff_S0 <- matrix(0, nrow = S, ncol = length(pvec))
  Gewekemat <- matrix(0, nrow = S, ncol = (K+3+2+2-1))
  acceptvec <- c()

  if(scenario == 1){
    betas_true  <- c(0.70,-1.15, 0.95)
    gammas_true <- c(-0.10, 0.25)
    regcoeffs_true <- c(betas_true, gammas_true)
  } else if (scenario == 2){
    betas_true  <- c(1.25, -0.75, 0.45)
    gammas_true <- c(-0.10, 0.20)
    regcoeffs_true <- c(betas_true, gammas_true)
  } else{
    stop("Scenario must be either 1 or 2")
  }

  #-- Progress bar
  progbar <- progress::progress_bar$new(
    format = crayon::white$yellow("Simulation in progress [:elapsed :spin] [:bar] :percent"),
    total = S,
    clear = FALSE
  )

  tic <- proc.time() #start clock

  if(!is.null(simnum)){
    set.seed(simnum)
  }

  for (s in 1:S) {

    # Generate survival data from mixture cure model
    simdat <- simdatmixcure(n = n, wshape = 1.45, wscale = 0.25,
                            setting = scenario)

    # Extract variables and construct formula
    simdatframe <- simdat$simdata
    formula <- Surv(tobs, event) ~ inci(x1 + x2) + late(z1 + z2)

    # Fit model with Griddy-Gibbs
    fitmcmc <- LangevinGibbs(formula = formula, data = simdatframe, K = K,
                           mcmcseed = NULL, penorder = 3,
                           mcmcsample = chainlength, burnin = burn,
                           progbar = "no")

    Gewekemat[s,] <- fitmcmc$Geweke
    acceptvec[s] <- fitmcmc$acceptrate

    # Fill host matrices
    thetamat[s,] <- fitmcmc$thetahat
    betamat[s,]  <- fitmcmc$betahat
    gammamat[s,] <- fitmcmc$gammahat

    CI90 <- fitmcmc$CI90
    CI95 <- fitmcmc$CI95

    for (j in 1:5) {
      coverage90[s, j] <- regcoeffs_true[j] >= CI90[j, 1] &&
        regcoeffs_true[j] <= CI90[j, 2]
      coverage95[s, j] <- regcoeffs_true[j] >= CI95[j, 1] &&
        regcoeffs_true[j] <= CI95[j, 2]
    }

    # Compute coverage probabilities for baseline survival at selected quantiles

    ## Inverse of true baseline survival to find quantiles of interest tq
    quantileS0 <- function(p) {
      val <- ((-1 / simdat$wscale) * log(p)) ^ (1 / simdat$wshape)
      return(val)
    }

    tq <- sapply(pvec, quantileS0)
    S0tq <- sapply(tq, simdat$S0)
    tu <- fitmcmc$tu

    basehaz <- function(thetahat, t){
      val <- exp(as.numeric(Rcpp_cubicBspline(t, lower = 0,
                          upper = tu, K = K) %*%  thetahat))
      return(val)
    }

    basesurv <- function(thetahat, t){
      val <- exp(-stats::integrate(basehaz, lower = 0, upper = t,
                                   thetahat = thetahat)$value)
      return(val)
    }

    chainl <- fitmcmc$mcmcsample-fitmcmc$burnin
    S0chains <- matrix(0, nrow = chainl, ncol = length(pvec))

    for(j in 1:ncol(S0chains)){
      S0chains[,j] <- apply(fitmcmc$thetachain[(fitmcmc$burnin+1):
                          fitmcmc$mcmcsample,],MARGIN = 1, basesurv, t=tq[j])
    }

    Diff_S0[s, ] <- colMeans(S0chains)-S0tq

    CI_S090 <- apply(S0chains, 2, "quantile", probs = c(0.05, 0.95))
    CI_S095 <- apply(S0chains, 2, "quantile", probs = c(0.025, 0.975))

    for(j in 1:length(tq)){
      S0cover90[s,j] <- pvec[j] >= CI_S090[1,j] && pvec[j] <= CI_S090[2,j]
      S0cover95[s,j] <- pvec[j] >= CI_S095[1,j] && pvec[j] <= CI_S095[2,j]
    }

    # Measuring CI width for baseline survival
    S090CIwidth[s, ] <- apply(CI_S090, 2, "diff")
    S095CIwidth[s,] <- apply(CI_S095, 2, "diff")

    zprofile <- c(0, 0.4)

    ## Inverse of survival for uncured to find quantiles of interest tqq
    quantileSu <- function(p) {
      val <- ((-1 / (simdat$wscale*exp(as.numeric(zprofile%*%simdat$gammas))))
              * log(p)) ^ (1 / simdat$wshape)
      return(val)
    }

    tqq <- sapply(pvec, quantileSu)
    Sutq <- sapply(tqq, simdat$S0)^(exp(as.numeric(zprofile%*%simdat$gammas)))


    Susurv <- function(thetagamm, t){
      thetahat <- thetagamm[1:K]
      gammahat <- thetagamm[((K+1):(K+2))]
      exponent <- as.numeric(exp(zprofile%*%gammahat))
      val <- (basesurv(thetahat, t))^(exponent)
      return(val)
    }
    Suchains <- matrix(0, nrow = chainl, ncol = length(pvec))

    for(j in 1:ncol(Suchains)){

      matthetagamm <- matrix(0, nrow = chainl, ncol = (K+2))
      matthetagamm[,(1:K)] <- fitmcmc$thetachain[(fitmcmc$burnin+1):
                                                   fitmcmc$mcmcsample,]
      matthetagamm[,((K+1):(K+2))] <- fitmcmc$gammachain[(fitmcmc$burnin+1):
                                                   fitmcmc$mcmcsample,]
      Suchains[,j] <- apply(matthetagamm, MARGIN = 1, Susurv, t=tqq[j])
    }

    CI_Su90 <- apply(Suchains, 2, "quantile", probs = c(0.05, 0.95))
    CI_Su95 <- apply(Suchains, 2, "quantile", probs = c(0.025, 0.975))

    for(j in 1:length(tqq)){
      Sucover90[s,j] <- pvec[j] >= CI_Su90[1,j] && pvec[j] <= CI_Su90[2,j]
      Sucover95[s,j] <- pvec[j] >= CI_Su95[1,j] && pvec[j] <= CI_Su95[2,j]
    }

    # Measuring CI width for survival of uncured
    Su90CIwidth[s, ] <- apply(CI_Su90, 2, "diff")
    Su95CIwidth[s,] <- apply(CI_Su95, 2, "diff")

     progbar$tick()

     realtimecoverage90 <- round(apply(coverage90,2,"sum")/s,3) * 100
     realtimecoverage95 <- round(apply(coverage95,2,"sum")/s,3) * 100
     realtimeS090 <- round(apply(S0cover90,2, "sum")/s,3) * 100
     realtimeS095 <- round(apply(S0cover95,2, "sum")/s,3) * 100
     realtimeSu90 <- round(apply(Sucover90,2, "sum")/s,3) * 100
     realtimeSu95 <- round(apply(Sucover95,2, "sum")/s,3) * 100

  }

  toc <- proc.time() - tic
  toc <- round(toc[3],3)

  ## Compute output statistics Mean, Bias, ESE, RMSE, CP90%, CP95%

  # Define metrics

  # Bias
  Bias <- function(bmat, btrue, S) {
    btruemat <- matrix(rep(btrue, S), ncol = length(btrue), byrow = T)
    diff <- bmat - btruemat
    biasvec <- apply(diff, 2, "mean")
    return(biasvec)
  }

  # Empirical standard error (ESE)
  ESE <- function(bmat, btrue, S) {
    bmean <- apply(bmat, 2, "mean")
    bmean_mat <- matrix(rep(bmean, S), ncol = length(btrue), byrow = T)
    ESEvec <- sqrt((1 / (S - 1)) * colSums((bmat - bmean_mat) ^ 2))
    return(ESEvec)
  }

  # Root mean square error (RMSE)
  RMSE <- function(bmat, btrue,S){
    btruemat <- matrix(rep(btrue, S), ncol = length(btrue), byrow = T)
    RMSEvec <- sqrt(apply((bmat - btruemat) ^ 2, 2, "mean"))
    return(RMSEvec)
  }

  # Mean CI width
  meanS0CI90 <- colMeans(S090CIwidth)
  meanS0CI95 <- colMeans(S095CIwidth)
  meanSuCI90 <- colMeans(Su90CIwidth)
  meanSuCI95 <- colMeans(Su95CIwidth)

  # Create output matrix
  simulres <- matrix(0, nrow = 5 , ncol = 8)
  colnames(simulres) <- c("Scenario","Parameters", "Mean", "Bias", "ESE",
                          "RMSE", "CP90","CP95")
  rownames(simulres) <- c("beta0","beta1","beta2","gamma1","gamma2")
  simulres[, 1] <- rep(scenario, 5)
  simulres[, 2] <- regcoeffs_true
  simulres[, 3] <- round(colMeans(cbind(betamat,gammamat)),3)
  simulres[, 4] <- round(c(Bias(betamat,betas_true,S),
                           Bias(gammamat,gammas_true,S)),3)
  simulres[, 5] <- round(c(ESE(betamat,betas_true,S),
                           ESE(gammamat,gammas_true,S)),3)
  simulres[, 6] <- round(c(RMSE(betamat, betas_true, S),
                           RMSE(gammamat,gammas_true, S)),3)
  simulres[, 7] <- round(colMeans(coverage90) * 100, 2)
  simulres[, 8] <- round(colMeans(coverage95) * 100, 2)

  # Bias of S0 at selected quantiles
  Bias_S0 <- colMeans(Diff_S0)

  # Compute coverage probabilities for baseline survival

  S0CP90 <- colMeans(S0cover90) * 100
  S0CP95 <- colMeans(S0cover95) * 100
  baselinecoverage <- matrix(0, nrow = 2, ncol = length(pvec))
  colnames(baselinecoverage) <- paste0("t", pvec)
  rownames(baselinecoverage) <- c("CP90%","CP95%")
  baselinecoverage[1, ] <- S0CP90
  baselinecoverage[2, ] <- S0CP95


  # Compute coverage probabilities for survival of uncured

  SuCP90 <- colMeans(Sucover90) * 100
  SuCP95 <- colMeans(Sucover95) * 100

  Sucoverage <- matrix(0, nrow = 2, ncol = length(pvec))
  colnames(Sucoverage) <- paste0("t", pvec)
  rownames(Sucoverage) <- c("CP90%","CP95%")
  Sucoverage[1, ] <- SuCP90
  Sucoverage[2, ] <- SuCP95

  # Print on screen
  cat(paste(rep("-",50),collapse = ""),"\n")
  cat("Simulation results for the Langevin-Gibbs sampler \n")
  cat(paste(rep("-",50),collapse = ""),"\n")
  cat("Scenario: ", scenario,          "\n")
  cat("Sample size: ", n,              "\n")
  cat("No. of B-splines: ", K, "\n")
  cat("No. of replications: ", S, "\n")
  cat(paste(rep("-",90),collapse = ""),"\n")
  print.table(as.matrix(simulres), digits = 3, justify = "left")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat("Estimated coverage probabilities of S0 at selected quantiles \n")
  cat(paste(rep("-",65),collapse = ""),"\n")
  print.table(baselinecoverage, digits = 3, justify = "left")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat("Estimated coverage probabilities of S(uncured) for z=(0,0.4) \n")
  cat(paste(rep("-",65),collapse = ""),"\n")
  print.table(Sucoverage, digits = 3, justify = "left")
  cat(paste(rep("-",90),collapse = ""),"\n")
  cat(paste0("Total elapsed time: ",toc, " seconds.\n"))
  cat(paste0("MCMC chain length: ",chainlength, ".\n"))
  cat(paste0("Burnin length: ",burn, ".\n"))

  outlist <- list(simulres = simulres,
                  S0coverage = baselinecoverage,
                  Sucoverage = Sucoverage,
                  elapsed = toc,
                  Geweke = Gewekemat,
                  acceptrate = acceptvec,
                  meanS0CI90 =  meanS0CI90,
                  meanS0CI95 = meanS0CI95,
                  meanSuCI90 = meanSuCI90,
                  meanSuCI95 = meanSuCI95,
                  Bias_S0 = Bias_S0)

}


















