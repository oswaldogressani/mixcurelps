#' @method plot LangevinGibbs
#' @export


# Plot method for an object of class LangevinGibbs

plot.LangevinGibbs <- function(x, param="beta1",
                             themetype = c("classic","gray","light","dark"),
                             tracecol = "darkblue",...){


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

  betadim <- ncol(x$X)
  gammadim <- ncol(x$Z)
  burn <- x$burnin
  M <- x$mcmcsample

  paramtype <- substr(param, 1, 1)
  if (paramtype == "b") {
    paramtype <- "beta"
    paramidx <- as.numeric(substr(param, 5, 5))
  } else if (paramtype == "g") {
    paramtype <- "gamma"
    paramidx <- as.numeric(substr(param, 6, 6))
  } else if (paramtype == "G"){
    paramtype <- "Geweke"
  } else if (paramtype == "l"){
    paramtype <- "lambda"
  }

  idxseq <- seq_len(M - burn)
  if (paramtype == "beta") {
    chain <- x$betachain[(burn + 1):M, (paramidx + 1)]
  } else if(paramtype == "gamma"){
    chain <- x$gammachain[(burn+1):M, paramidx]
  } else if (paramtype == "Geweke"){
    chain <- x$Geweke
    idxseq <- seq_len(length(x$Geweke))
  } else if (paramtype == "lambda"){
    chain <- x$lambdachain[(burn+1):M]
    paramidx <- ""
  }


  chaindat <- data.frame(idxseq, chain)

  if(paramtype == "Geweke"){
    ggplot2::ggplot(data = chaindat, ggplot2::aes(x=idxseq, y = chain)) +
      ggplot2::xlab("Latent index") +
      ggplot2::ylab("Geweke z-score") +
      themeval +
      ggplot2::ggtitle("Geweke diagnostics") +
      ggplot2::geom_point(size = 2, colour=tracecol) +
      ggplot2::geom_hline(yintercept = -1.96, linetype = "dashed",
                          size = 1.1, colour="red") +
      ggplot2::geom_hline(yintercept = 1.96, linetype = "dashed",
                          size = 1.1, colour="red") +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12),
        axis.title.y = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  } else {
    ggplot2::ggplot(data = chaindat, ggplot2::aes(x=idxseq, y = chain)) +
      ggplot2::geom_line(colour=tracecol) +
      ggplot2::xlab("Iterations") +
      ggplot2::ylab("") +
      themeval +
      ggplot2::ggtitle(paste0("Trace of ",paramtype," ",paramidx)) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  }


}
