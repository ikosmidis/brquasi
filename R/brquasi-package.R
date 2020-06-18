#' Bias reduction in quasi likelihood estimation
#'
#' Estimation and inference from quasi likelihoods based on various
#' methods for bias reduction. The 'brquasi' fitting method can
#' achieve reduction of bias either by using the empirical adjustments
#' to the quasi score functions in Kosmidis & Lunardon (2020) or by
#' direct subtraction of an empirical estimate of the bias of the
#' quasi-likelihood estimator. Estimation in all cases takes place via
#' a quasi Newton-Raphson algorithm.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso
#'
#' [`brquasi_fit()`], [`brglm2::brglm_fit()`]
#'
#' @references
#'
#' Kosmidis I, Lunardon N (2020). Empirical bias-reducing adjustments
#' to estimating functions. *ArXiV eprints*. \url{https://arxiv.org/abs/2001.03786}
#'
#' @docType package
#' @name brquasi
#' @import stats
#' @importFrom numDeriv grad
#'
NULL

## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if (getRversion() >= "2.15.1") globalVariables(c("n"))
