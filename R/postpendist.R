#' Plot the normalized approximate posterior of the penalty parameter.
#'
#' @description
#' Plots the normalized approximate postrior of the penalty parameter based
#' on a an object of class \code{lpsmc}.
#'
#' @param x An object of class \code{lpsmc}.
#' @param low The lower bound on the x-axis (in log scale).
#' @param up The upper bound on the x-axis (in log scale).
#' @param themetype The theme, either "classic", "gray", "light" or "dark".
#'
#' @examples
#' ### Posterior penalty distribution for breast cancer dataset
#' rm(list=ls())
#' data("breastcancer")
#' formula <- Surv(tobs, delta) ~ inci(AGE + ER) + late(AGE + ER)
#' fitcancer <- lpsmc(formula = formula, data = breastcancer, K = 20,
#'                     stepsize = 0.1)
#' postpendist(fitcancer, 5, 12, themetype = "gray")
#'
#' @author Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
#'
#' @export

postpendist <- function(x, low, up, themetype=c("classic","gray","light","dark")){

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

  vgrid <- seq(low, up, length = 200)
  dv <- vgrid[2] - vgrid[1]
  logpvgrid <- unlist(lapply(lapply(vgrid, x$logpv2), "[[", 1))
  pvv <- exp(logpvgrid - max(logpvgrid))
  cnorm <- 1 / sum(pvv * dv)
  pv <- cnorm * pvv
  pvimg <- data.frame(vgrid, pv)
  skplot <- ggplot2::ggplot(data = pvimg, ggplot2::aes(x = vgrid, y = pv))
  skplot + ggplot2::geom_line(colour="darkblue", size=1.2) +
    ggplot2::xlab("v") +
    ggplot2::ylab("Normalized approximate posterior to p(v|D)") +
    ggplot2:: theme(axis.title.x = ggplot2::element_text(size = 14),
          axis.title.y = ggplot2::element_text(size = 14)) +
    ggplot2::geom_vline(xintercept = x$vhat, linetype = "dashed", size = 1) +
    themeval

}



