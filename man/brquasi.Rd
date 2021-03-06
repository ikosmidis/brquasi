% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brquasi-package.R
\docType{package}
\name{brquasi}
\alias{brquasi}
\title{Bias reduction in quasi likelihood estimation}
\description{
Estimation and inference from quasi likelihoods based on various
methods for bias reduction. The 'brquasi' fitting method can
achieve reduction of bias either by using the empirical adjustments
to the quasi score functions in Kosmidis & Lunardon (2020) or by
direct subtraction of an empirical estimate of the bias of the
quasi-likelihood estimator. Estimation in all cases takes place via
a quasi Newton-Raphson algorithm.
}
\references{
Kosmidis I, Lunardon N (2020). Empirical bias-reducing adjustments
to estimating functions. *ArXiV eprints*. \url{https://arxiv.org/abs/2001.03786}
}
\seealso{
[`brquasi_fit()`], [`brglm2::brglm_fit()`]
}
\author{
Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
}
