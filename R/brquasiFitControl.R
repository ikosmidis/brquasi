# Copyright (C) 2020 Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


#' Auxiliary function for [`glm()`] fitting using the
#' [`brquasiFit()`] method.
#'
#' Typically only used internally by `brquasiFit()`, but may be used
#' to construct a `control` argument.
#'
#' @aliases quasi_control
#' @param epsilon positive convergence tolerance epsilon. Default is
#'     `1e-04`.
#' @param maxit integer giving the maximal number of iterations
#'     allowed. Default is `500`.
#' @param trace logical indicating if output should be produced for
#'     each iteration. Default is `FALSE`.
#' @param type the type of fitting method to be used. The options are
#'     `"M"` (standard M-estimation based on quasi-likelihoods),
#'     `"iRBM"` (implicit reduced-bias M-estiamtion: empirically
#'     adjusted quasi likelihood equations for mean-bias reduction;
#'     default), `eRBM` (explicit reduced-bias M-estimation:
#'     correction of asymptotic mean bias using empirical bias
#'     estimates). `"MPQL_trace"` (maximum penalized quasi likelihood
#'     estimation; see Details).
#' @param only_beta: Should RBM estimation be used for improving
#'     estimation of the regression coefficients only? Default is
#'     `TRUE`.
#' @param slowit a positive real used as a multiplier for the
#'     stepsize. The smaller it is the smaller the steps are. Default
#'     is `1`.
#' @param max_step_factor the maximum number of step halving steps to
#'     consider. Default is `12`.
#' @param response_adjustment a (small) positive constant or a vector
#'     of such. Default is `NULL`. See Details.
#' @param a mutliple to the additive trace penalty to the quasi
#'     likelihood when `type = "MPQL_trace"`. See Details.
#' @param lambda a ridge adjustment to be added to the diagonal of the
#'     jacobian before inverting it for the computation of the step
#'     size in the quasi-Fisher iteration. Default is `1e-10`. See
#'     Details.
#' @param disp_factor factor by which to divide the sum of pearson
#'     residuals when estimating the dispersion. Defauls is `"n-p"`
#'     which corresponds to "number of observations" minus
#'     "number of parameters". See Details.
#'
#' Describe quasi Fisher iteration
#'
#' Describe what RBM
#'
#' Describe disp_factor
#'
#' @details
#'
#' @export
brquasiFitControl <- function(epsilon = 1e-04, maxit = 500,
                              trace = FALSE,
                              type = c("iRBM", "M", "eRBM", "MPQL_trace"),
                              slowit = 1,
                              response_adjustment = 0,
                              max_step_factor = 12,
                              only_beta = TRUE,
                              lambda = 1e-10,
                              disp_factor = c("n-p", "n")) {
    type <- match.arg(type)
    disp_factor <- match.arg(disp_factor)
    if (isTRUE(epsilon <= 0)) {
        stop("value of 'epsilon' must be > 0")
    }
    if (isTRUE(lambda < 0)) {
        stop("value of 'lambda' must be >= 0")
    }
    list(epsilon = epsilon,
         maxit = maxit,
         trace = trace,
         response_adjustment = response_adjustment,
         type = type,
         slowit = slowit,
         only_beta = only_beta,
         lambda = lambda,
         max_step_factor = max_step_factor,
         disp_factor = disp_factor)
}

