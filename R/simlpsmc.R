#' Assessing the statistical performance of lpsmc.
#'
#'
#' @description
#' This routine can be used to assess the statistical performance of the
#' \code{lpsmc} routine in a simulation setting. The scenario (1 or 2)
#' determines the data generating process and the total number of replications
#' can be fixed by the user.
#'
#' @param n Sample size.
#' @param K Number of B-spline basis functions.
#' @param scenario Either 1 or 2.
#' @param S Total number of replications.
#' @param exactrep Exactly replicate results.
#' @param simnum A seed for reproducibility.
#' @param themetype The theme of the plot either "classic", "gray","light"
#'  or "dark".
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @examples
#' ### Scenario 1,  n=300
#' sim1 <- simlpsmc(n = 300, scenario = 1, S = 10, exactrep = TRUE)
#' suppressWarnings(print(sim1$S0plot))
#' sim1$ASEplot
#'
#' @export

simlpsmc <- function(n = 300, K = 15, scenario = 1, S = 500, exactrep = FALSE,
                     simnum = NULL,
                     themetype = c("classic","gray","light","dark")){

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

  if(exactrep == TRUE) {
    # Setting seeds (for reproducibility)
    if (scenario == 1 && n == 300) {
      set.seed(1989) # Fall of the Berlin wall
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

  ASEpvec <- c()
  coverageincid90 <- c()
  coverageincid95 <- c()
  covercure90 <- c()
  covercure95 <- c()

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

    # Fit model
    fitlpsmc <- lpsmc(formula = formula, data = simdatframe, K = K)

    # Fill host matrices
    thetamat[s,]   <- fitlpsmc$thetahat
    betamat[s,]    <- fitlpsmc$betahat
    gammamat[s,]   <- fitlpsmc$gammahat

    CI90 <- fitlpsmc$CI90
    CI95 <- fitlpsmc$CI95

    for (j in 1:5) {
      coverage90[s, j] <- regcoeffs_true[j] >= CI90[j, 1] &&
        regcoeffs_true[j] <= CI90[j, 2]
      coverage95[s, j] <- regcoeffs_true[j] >= CI95[j, 1] &&
        regcoeffs_true[j] <= CI95[j, 2]
    }

    # Compute coverage probability of incidence p(x)
    xmean_true <- c(1, 0, 0.5)

    # True incidence at mean covariate vector x
    p_meancovar <- fitlpsmc$px(simdat$betas,xmean_true)
    coverageincid90[s] <- as.numeric(p_meancovar >= fitlpsmc$CIincid90[1] &&
                                    p_meancovar <= fitlpsmc$CIincid90[2])
    coverageincid95[s] <- as.numeric(p_meancovar >= fitlpsmc$CIincid95[1] &&
                                    p_meancovar <= fitlpsmc$CIincid95[2])


    # Compute coverage probability of cure rate 1-p(x)
    pcure <- 1-p_meancovar
    covercure90[s] <- as.numeric(pcure >= fitlpsmc$CIcure90[1] &&
                                   pcure <= fitlpsmc$CIcure90[2])
    covercure95[s] <- as.numeric(pcure >= fitlpsmc$CIcure95[1] &&
                                   pcure <= fitlpsmc$CIcure95[2])

    # Compute coverage probabilities for baseline survival at selected quantiles

    ## Inverse of true baseline survival to find quantiles of interest tq
    quantileS0 <- function(p) {
      val <- ((-1 / simdat$wscale) * log(p)) ^ (1 / simdat$wshape)
      return(val)
    }

    tq <- sapply(pvec, quantileS0)
    S0tq <- sapply(tq, simdat$S0)

    # Function for CI for S0 at t
    S0CI <- function(t, alpha, thetahat){
      cum_tt <- fitlpsmc$cumult(t = t, theta = thetahat, k = 0, l = 0)
      gstar  <- log(cum_tt)
      grad_g <- fitlpsmc$cumult(t = t, theta = thetahat, k = 1, l = 0)[-K] *
        (1 / cum_tt)
      Sig_thetahat <- fitlpsmc$Covhat[1:(K-1), 1:(K-1)]
      qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
      post_sd <- sqrt(as.numeric(t(grad_g) %*% Sig_thetahat %*% grad_g))
      CI_g <- c(gstar - qz_alpha * post_sd, gstar + qz_alpha * post_sd)
      CIS <- rev(exp(-exp(CI_g)))
      return(CIS)
    }

    # Compute 90% and 95% coverage probability

    CI_S090 <- sapply(tq, S0CI, alpha = 0.10, thetahat = fitlpsmc$thetahat)
    colnames(CI_S090) <- pvec
    CI_S095 <- sapply(tq, S0CI, alpha = 0.05, thetahat = fitlpsmc$thetahat)
    colnames(CI_S095) <- pvec

    for(j in 1:length(tq)){
      S0cover90[s,j] <- pvec[j] >= CI_S090[1,j] && pvec[j] <= CI_S090[2,j]
      S0cover95[s,j] <- pvec[j] >= CI_S095[1,j] && pvec[j] <= CI_S095[2,j]
    }

    # Measuring CI width for baseline survival
    S090CIwidth[s, ] <- apply(CI_S090, 2, "diff")
    S095CIwidth[s,] <- apply(CI_S095, 2, "diff")

    # Compute coverage probabilities for survival of uncured at selected
    # quantiles and for covariate profile z=(0,0.4)

    zprofile <- c(0, 0.4)

    ## Inverse of survival for uncured to find quantiles of interest tqq
    quantileSu <- function(p) {
      val <- ((-1 / (simdat$wscale*exp(as.numeric(zprofile%*%simdat$gammas))))
              * log(p)) ^ (1 / simdat$wshape)
      return(val)
    }

    tqq <- sapply(pvec, quantileSu)
    Sutq <- sapply(tqq, simdat$S0)^(exp(as.numeric(zprofile%*%simdat$gammas)))

    SuCI <- function(t, alpha, thetahat){
      cum_tt <- fitlpsmc$cumult(t = t, theta = thetahat, k = 0, l = 0)
      gstar  <- as.numeric(zprofile %*% fitlpsmc$gammahat) + log(cum_tt)
      grad_g <- c(fitlpsmc$cumult(t = t, theta = thetahat, k = 1, l = 0)[-K] *
                    (1 / cum_tt), zprofile)
      p <- fitlpsmc$p
      q<- fitlpsmc$q
      dimlat <- K+p+q
      SSighat <- matrix(0, nrow = (K - 1) + q, ncol = (K - 1) + q)
      SSighat[1:(K-1),1:(K-1)] <- fitlpsmc$Covhat[1:(K-1), 1:(K-1)]
      SSighat[1:(K-1),K:(K+1)] <- fitlpsmc$Covhat[1:(K - 1), (K + p + 1):dimlat]
      SSighat[K:(K + 1), 1:(K - 1)] <-
        t(fitlpsmc$Covhat[1:(K - 1), (K + p + 1):dimlat])
      SSighat[K:(K + 1), K:(K + 1)] <-
        fitlpsmc$Covhat[(K + p + 1):dimlat, (K + p + 1):dimlat]
      qz_alpha  <- stats::qnorm(alpha * 0.5, lower.tail = FALSE)
      post_sd <- sqrt(as.numeric(t(grad_g) %*% SSighat %*% grad_g))
      CI_g <- c(gstar - qz_alpha * post_sd, gstar + qz_alpha * post_sd)
      CIS <- rev(exp(-exp(CI_g)))
      return(CIS)
    }

    CI_Su90 <- sapply(tqq, SuCI, alpha = 0.10, thetahat = fitlpsmc$thetahat)
    colnames(CI_Su90) <- pvec
    CI_Su95 <- sapply(tqq, SuCI, alpha = 0.05, thetahat = fitlpsmc$thetahat)
    colnames(CI_Su95) <- pvec

    for(j in 1:length(tqq)){
      Sucover90[s,j] <- pvec[j] >= CI_Su90[1,j] && pvec[j] <= CI_Su90[2,j]
      Sucover95[s,j] <- pvec[j] >= CI_Su95[1,j] && pvec[j] <= CI_Su95[2,j]
    }

    # Measuring CI width for survival of uncured
    Su90CIwidth[s, ] <- apply(CI_Su90, 2, "diff")
    Su95CIwidth[s,] <- apply(CI_Su95, 2, "diff")

    progbar$tick()
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

  # Create plot for baseline survival

  S0tdom <- seq(0, fitlpsmc$tu, length = 200)
  S0tar <- sapply(S0tdom, simdat$S0)
  S0hat <- matrix(0, nrow = S, ncol = length(S0tdom))

  basehaz <- function(t, thetahat){
    val <- exp(as.numeric(Rcpp_cubicBspline(t, lower = 0,
             upper = fitlpsmc$tu, K = K) %*%  thetahat))
    return(val)
  }

  basesurv <- function(t, thetahat){
    val <- exp(-stats::integrate(basehaz, lower = 0, upper = t,
                                 thetahat = thetahat)$value)
    return(val)
  }

  for(s in 1:S){
    Shats <-  sapply(S0tdom, basesurv, thetahat = thetamat[s, ])
    Shats[1] <- 1
    S0hat[s,] <- Shats
  }

  S0hat <- t(S0hat)
  S0dat <- data.frame(cbind(S0tdom,S0tar,S0hat,apply(S0hat,1,"median")))
  S0baseplot <- ggplot2::ggplot(data = S0dat, ggplot2::aes(x=S0tdom))
  S0colnames <- names(S0dat)

  S0baseplot <- S0baseplot +
    ggplot2::geom_line(ggplot2::aes(y=S0dat[,2]), colour="black", size = 1.2) +
     ggplot2::xlim(0,8.5)


  for(s in 3:(S+2)){
    S0baseplot <-
      S0baseplot + eval(parse(
        text = paste(
          "ggplot2::geom_line(ggplot2::aes(y=",
          S0colnames[s],
          "),colour='#c1cad7')",
          sep = ""
        )
      ))
  }

  S0baseplot <- S0baseplot + ggplot2::geom_line(ggplot2::aes(y=S0dat[,2]),
                                                colour="black", size = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y=S0dat[,ncol(S0dat)]),
                       colour="firebrick2", size = 1.1, linetype = "dashed") +
    themeval +
    ggplot2::xlab("t") +
    ggplot2::ylab(expression(S[0](t))) +
    ggplot2::ggtitle(paste0("Scenario ",scenario,", n=",n))+
    ggplot2:: theme(axis.title.x = ggplot2::element_text(size = 14),
                    axis.title.y = ggplot2::element_text(size = 14),
                    plot.title = ggplot2::element_text(hjust = 0.5, size = 15),
                    axis.text.x = ggplot2::element_text(size=12),
                    axis.text.y = ggplot2::element_text(size=12))

  # Compute bias of baseline survival at quantiles tq

  Diff_S0 <- matrix(0, nrow = S, ncol = length(pvec))
  for(s in 1:S){
    Diff_S0[s, ]<-  sapply(tq, basesurv, thetahat = thetamat[s, ]) - S0tq
  }
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

  # Average squared error for proportion of uncured (incidence)
  x1m <- seq(-1.5, 1.5, by = 0.001)
  x2m <- c(rep(0, length(x1m)), rep(1, length(x1m)))
  Xx <- cbind(rep(1, length(x2m)), rep(x1m, 2), x2m)
  ptrue <- fitlpsmc$px(simdat$betas, Xx)

  for (s in 1:S) {
    phat <- fitlpsmc$px(betamat[s,], Xx)
    ASEpvec[s] <- mean((phat - ptrue) ^ 2)
  }

  ASEdatframe <- data.frame(as.factor(rep(1, S)), ASEpvec)
  ASEplot <- ggplot2::ggplot(data = ASEdatframe,
                             ggplot2::aes(x = ASEdatframe[,1], y = ASEpvec)) +
    ggplot2::geom_boxplot(color = "black", fill = "#DE4E4E", width = 0.6) +
    ggplot2::xlab("") +
    ggplot2::ylab("ASE(p)") +
    ggplot2::ggtitle(paste0("Scenario ",scenario,", n=",n)) +
    themeval +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title.y = ggplot2::element_text(size = 14),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 15))


  # Compute coverage of incidence and cure
  CI90_incidence <- round(mean(coverageincid90) * 100, 2)
  CI95_incidence <- round(mean(coverageincid95) * 100, 2)
  CI90_cure <- round(mean(covercure90) * 100, 2)
  CI95_cure <- round(mean(covercure95) * 100, 2)


  # Print on screen
  cat(paste(rep("-",50),collapse = ""),"\n")
  cat("Simulation results for LPSMC \n")
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

  outlist <- list(simulres = simulres,
                  S0plot = S0baseplot,
                  ASEplot = ASEplot,
                  ASEpvec = ASEpvec,
                  S0coverage = baselinecoverage,
                  Sucoverage = Sucoverage,
                  cover_inci90 = CI90_incidence,
                  cover_inci95 = CI95_incidence,
                  cover_cure90 = CI90_cure,
                  cover_cure95 = CI95_cure,
                  meanS0CI90 =  meanS0CI90,
                  meanS0CI95 = meanS0CI95,
                  meanSuCI90 = meanSuCI90,
                  meanSuCI95 = meanSuCI95,
                  Bias_S0 = Bias_S0,
                  elapsed = toc)

}



