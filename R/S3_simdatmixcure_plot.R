#' @method plot simdatmixcure
#' @export

# Plot method for an object of class simdatmixcure
plot.simdatmixcure <- function(x, ...) {


  # tdom <- seq(0, x$tup, length = 1000)

  # graphics::plot(x$fitKM, mark.time = TRUE, mark = "x", xlab = "t",
  #                ylab = expression(S[0](t)), main = "Baseline survival",
  #                cex.main = 0.9)
  # graphics::abline(v = x$plateau, lty = 2, lwd = 2, col = "orange")
  # graphics::lines(tdom, sapply(tdom, x$S0), type = "l", col = "blue")
  # graphics::legend("topright", lty = c(1,1,2),
  #                  col = c("black", "blue", "orange"),
  #                  c("Kaplan-Meier", "Weibull baseline survival",
  #                    "Start of plateau"), bty = "n", cex = 0.8)


  # With survminer
  tobs <- x$tobs
  status <- x$delta
  dataKapM <- data.frame(tobs, status)
  fitKapM <- survival::survfit(survival::Surv(tobs, status) ~ 1,
                               data = dataKapM)
  plotsurv <- survminer::ggsurvplot(fitKapM,
                        data = dataKapM,
                        censor.shape="x",
                        censor.size = 5.5,
                        size = 1,
                        palette = "#0089FF",
                        conf.int = TRUE,
                        font.tickslab = c(14,"darkblue"),
                        font.x =c(14,"black"),
                        font.y = c(14,"black"),
                        ggtheme = ggplot2::theme_light(),
                        risk.table = "percentage",
                        risk.table.col ="darkblue",
                        legend="none",
                        legend.title="",
                        tables.theme = survminer::theme_cleantable()
                        )
  plotsurv$table <- plotsurv$table + ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size=14),
    axis.title.x = ggplot2::element_text(size=14)
  )

  plotsurv $plot<- plotsurv$plot + ggplot2::geom_vline(xintercept = x$plateau,
                                 linetype = "dashed", size = 1,
                                 colour = "#15BA57")

  plotsurv





}
