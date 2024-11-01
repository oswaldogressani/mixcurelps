% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simdatmixcure.R
\name{simdatmixcure}
\alias{simdatmixcure}
\title{Simulation of survival data to fit mixture cure models.}
\usage{
simdatmixcure(n, wscale, wshape, setting)
}
\arguments{
\item{n}{Sample size.}

\item{wscale}{The positive scale parameter of the Weibull distribution used
to generate the survival times of the uncured subjects.}

\item{wshape}{The positive shape parameter of the Weibull distribution used
to generate the survival times of the uncured subjects.}

\item{setting}{The setting under which survival times will be generated. If
\code{setting = 1}, the coefficients of the incidence part are
\emph{beta0=0.70, beta1=-1.15 and beta2=0.95} and the coefficients of the
latency part are \emph{gamma1=-0.10 and gamma2=0.25}. If
\code{setting = 2}, the coefficients of the incidence part are
\emph{beta0=1.25, beta1=-0.75 and beta2=0.45} and the coefficients of the
latency part are \emph{gamma1=-0.10 and gamma2=0.20}.}
}
\value{
An object of class \code{simixcure} containing different objects of the
simulated dataset. Details can be found by typing ?simixcure.object.
}
\description{
This routines simulates survival data with a cure fraction. The data is
simulated according to a mixture cure model with a logistic link for the
incidence part of the model. The incidence part includes two covariates
following a standard normal and a Bernoulli distribution with success
probability equal to 0.5. The latency part assumes a Weibull-Cox model with
two covariates (a standard normal variate and a Bernoulli variate with
success probability equal to 0.4).
}
\examples{
### Simulate a sample of size n=300 under Scenario 1.
set.seed(4408)
simdat <- simdatmixcure(n = 300, wshape = 1.45, wscale = 0.25, setting = 1)
plot(simdat) # Plot the baseline survival and Kaplan-Meier curve
simdat$info  # Print information on Cure and Censoring levels

}
\seealso{
\link{simdatmixcure.object}
}
\author{
Oswaldo Gressani \email{oswaldo_gressani@hotmail.fr} .
}
