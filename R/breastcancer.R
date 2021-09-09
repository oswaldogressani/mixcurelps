#' Breast cancer data.
#'
#' @docType data
#'
#' @description Breast cancer data from the \code{breastCancerVDX} package.
#'
#' @usage data(breastcancer)
#'
#' @format A data frame with 286 rows and 4 columns.
#' \describe{
#'  \item{\code{tobs}}{Distant-metastasis-free survival (in days).}
#'  \item{\code{delta}}{Event indicator \code{1}=death or relapse, \code{0}=censored.}
#'  \item{\code{AGE}}{Age of patients.}
#'  \item{\code{ER}}{Estrogen receptor \code{0}="<=10fmol", \code{1}=">10fmol".}
#' }
#'
#'
#' @source \url{https://doi.org/doi:10.18129/B9.bioc.breastCancerVDX}
#'
#' @references  Schroeder M, Haibe-Kains B, Culhane A, Sotiriou C, Bontempi G,
#'  Quackenbush J (2021). breastCancerVDX: Gene expression datasets published
#'  by Wang et al. [2005] and Minn et al. [2007] (VDX). R package version 1.30.0.
"breastcancer"
