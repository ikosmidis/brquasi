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

#' Fitting function for [`glm()`] for reduced-bias
#' estimation and inference using quasi likelihoods
#'
#' `brquasiFit` is a fitting method for [`glm()`] for
#' mean-bias reduction in quasi-likelihood estimation, using the
#' implicit and explicit reduced-bias M-estimators in Kosmidis &
#' Lunardon (2020).
#'
#' Estimation is performed using a quasi Newton-Raphson iteration
#' (similar to the quasi Fisher-scoring iteration described, for example, in
#' `vignette("iteration", "brglm2")`, which, in the case of
#' mean-bias reduction, resembles an iterative correction of the
#' asymptotic bias of the Fisher scoring iterates.
#'
#' @inheritParams stats::glm.fit
#' @aliases brquasi_fit
#' @param x `x` is a design matrix of dimension `n * p`.
#' @param y `y` is a vector of observations of length `n`.
#' @param family any one of the [`quasi()`] families.
#' @param control a list of parameters controlling the fitting
#'     process. See [`brquasiControl()`] for details.
#' @param start starting values for the parameters in the linear
#'     predictor. If `NULL` (default) then the maximum likelihood
#'     estimates are calculated and used as starting values.
#' @param mustart applied only when start is not `NULL`. Starting
#'     values for the vector of means to be passed to
#'     `glm.fit` when computing starting values using
#'     maximum likelihood.
#' @param etastart applied only when start is not
#'     `NULL`. Starting values for the linear predictor to be
#'     passed to `glm.fit` when computing starting values
#'     using maximum likelihood.
#' @param singular.ok logical. If `FALSE`, a singular model is an
#'     error.
#'
#' @details
#'
#' @author Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [`glm.fit()`] [`glm()`] [`quasi()`]
#'
#' @export
brquasiFit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                        mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                        control = list(), intercept = TRUE,
                        singular.ok = TRUE)
{

    ## key_quantities, grad, info and bias are ALWAYS in beta, dispersion parameterization
    key_quantities <- function(pars, y, level = 0, qr = TRUE) {
        betas <- pars[seq.int(nvars)]
        dispersion <- pars[nvars + 1]
        prec <- 1/dispersion
        etas <- drop(x %*% betas + offset)
        mus <- linkinv(etas)
        out <- list(precision = prec,
                    betas = betas,
                    dispersion = dispersion,
                    etas = etas,
                    mus = mus)
        mean_quantities <- function(out) {
            d1mus <- mu.eta(etas)
            d2mus <- d2mu.deta(etas)
            d3mus <- d3mu.deta(etas)
            varmus <- variance(mus)
            working_weights <- weights * d1mus^2 / varmus
            wx <- sqrt(working_weights) * x
            out$d1mus <- d1mus
            out$d2mus <- d2mus
            out$d3mus <- d3mus
            out$varmus <- varmus
            out$d1varmus <- d1variance(mus)
            out$d2varmus <- d2variance(mus)
            out$working_weights <- working_weights
            if (qr) out$qr_decomposition <- qr(wx)
            out
        }
        dispersion_quantities <- function(out) {
            zetas <- -weights * prec
            out$zetas <- zetas
            ## Evaluate the derivatives of the a function only for
            ## observations with non-zero weight
            out$deviance_residuals <- dev.resids(y, mus, weights)
            out
        }
        if (level == 0) {
            out <- mean_quantities(out)
        }
        if (level == 1) {
            out <- dispersion_quantities(out)
        }
        if (level > 1) {
            out <- mean_quantities(out)
            out <- dispersion_quantities(out)
        }
        out
    }

    gradient <- function(pars, fit = NULL, contributions = FALSE, only_beta = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
        }
        with(fit, {
            score_components_beta <- precision * weights * d1mus  * (y - mus) / varmus * x
            if (only_beta) {
                score_components <- score_components_beta
            }
            else {
                score_components_dispersion <- weights / varmus * (y - mus)^2 -
                    dispersion * (disp_factor)/nobs
                score_components <- cbind(score_components_beta, score_components_dispersion)
            }
            if (contributions) {
                score_components
            }
            else {
                colSums(score_components)
            }
        })
    }

    ## The gradient of the set of estimating functions that gives rise
    ## to the pearson esitmator of the dispersion parameter
    jmat <- function(pars, fit = NULL, only_beta = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
        }
        with(fit, {
            bmus <- weights * d1mus / varmus
            cmus <- weights / varmus
            d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
            d1cmus <- - weights * d1mus * d1varmus / varmus^2
            emus <- (y - mus)
            qmus <- (bmus * d1mus - d1bmus * emus) / dispersion
            pmus <- bmus * emus / dispersion^2
            rmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2
            jbetabeta <- crossprod(x * qmus, x)
            if (only_beta) {
                return(jbetabeta)
            }
            else {
                jbetaphi_components <- pmus * x
                jphibeta_components <- rmus * x
                jphiphi_components <- disp_factor
                jbetaphi <- .colSums(jbetaphi_components, nobs, nvars, TRUE)
                jphibeta <- .colSums(jphibeta_components, nobs, nvars, TRUE)
                jphiphi <- sum(jphiphi_components)
                out <- rbind(cbind(jbetabeta, jbetaphi),
                             c(jphibeta, jphiphi))
                return(out)
            }
        })
    }

    emat <- function(pars, fit = NULL, only_beta = FALSE) {
        grad <- gradient(pars, fit = fit, contributions = TRUE, only_beta = only_beta)
        crossprod(grad)
    }

    umat <- function(pars, dispersion_block = FALSE, fit = NULL, only_beta = FALSE){
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
        }
        with(fit, {
            emus <- (y - mus)
            if (dispersion_block) {
                cmus <- weights / varmus
                d1cmus <- - weights * d1mus * d1varmus / varmus^2
                d2cmus <- - weights * (d2varmus * d1mus^2 + d1varmus * d2mus - 2 * d1mus^2 * d1varmus^2 / varmus) / varmus^2
                gmus <- d2cmus * emus^2 - 4 * d1cmus * d1mus * emus - 2 * cmus * d2mus * emus + 2 * cmus * d1mus^2
                ubetabeta <- crossprod(x * gmus, x)
                if (only_beta) {
                    out <- ubetabeta
                }
                else {
                    out <- cbind(rbind(ubetabeta, 0), 0)
                }
            }
            else {
                bmus <- weights * d1mus / varmus
                d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
                d2bmus <- weights * (d3mus / varmus - (3 * d1mus * d2mus * d1varmus + d1mus^3 * d2varmus) / varmus^2 +  d1mus^3 * d1varmus^2 / varmus^3)
                qmus <- (bmus * d1mus - d1bmus * emus) / dispersion
                pmus <- bmus * emus / dispersion^2
                d1qmus <- (2 * d1bmus * d1mus + bmus * d2mus - d2bmus * emus) / dispersion

                if (only_beta) {
                    out <- lapply(seq.int(nvars), function(j) {
                        - crossprod(x * d1qmus * x[, j], x)
                    })
                }
                else {
                    out <- lapply(seq.int(nvars), function(j) {
                        ubetabeta <-  - crossprod(x * d1qmus * x[, j], x)
                        ubetaphi <- crossprod(x * qmus, x)[j, ] / dispersion
                        uphiphi <- 2 * .colSums(pmus * x, nobs, nvars, TRUE)[j] / dispersion
                        rbind(cbind(ubetabeta, ubetaphi),
                              c(ubetaphi, uphiphi))
                    })
                }
            }
            return(out)
        })
    }

    dmat <- function(pars, dispersion_block = FALSE, fit = NULL, only_beta = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
        }
        with(fit, {
            emus <- (y - mus)
            if (only_beta) {
                score_components <- precision * weights * d1mus  * emus / varmus * x
            }
            else {
                score_components <- cbind(precision * weights * d1mus  * emus / varmus * x,
                                          weights / varmus * emus^2 - dispersion * (disp_factor)/nobs)
            }
            if (dispersion_block) {
                cmus <- weights / varmus
                d1cmus <- - weights * d1mus * d1varmus / varmus^2
                rmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2
                if (only_beta) {
                    out <- - crossprod(x * rmus, score_components)
                }
                else {
                    out <- - crossprod(cbind(x * rmus, (disp_factor)/nobs), score_components)
                }
            }
            else {
                bmus <- weights * d1mus / varmus
                d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
                qmus <- (bmus * d1mus - d1bmus * emus) / dispersion
                pmus <- bmus * emus / dispersion^2
                if (only_beta) {
                    out <- lapply(seq.int(nvars), function(r) {
                        - crossprod(x * qmus  * x[, r], score_components)
                    })
                }
                else {
                    out <- lapply(seq.int(nvars), function(r) {
                        - crossprod(cbind(x * qmus, pmus)  * x[, r], score_components)
                    })
                }
            }
            return(out)
        })
    }

    hessian <- function(pars, fit = NULL, only_beta = FALSE, lambda = 1e-06) {
        j <- jmat(pars, fit = fit, only_beta = only_beta)
        -(j + lambda * diag(nvars + !only_beta))
    }

    RBM_adjustment <- function(pars, fit = NULL,
                                             jm = NULL, em = NULL, sandwich = NULL, only_beta = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = 2, qr = TRUE)
        }
        inv_jm <- chol2inv(chol(jmat(pars, fit = fit, only_beta = only_beta)))
        em <- emat(pars, fit = fit, only_beta = only_beta)
        sandm <- inv_jm %*% em %*% inv_jm
        dm <- dmat(pars, dispersion_block = FALSE, fit = fit, only_beta = only_beta)
        um <- umat(pars, dispersion_block = FALSE, fit = fit, only_beta = only_beta)
        abeta <- sapply(seq.int(nvars), function(r) {
            -sum(inv_jm * t(dm[[r]])) - sum(sandm * t(um[[r]])) / 2
        })
        if (only_beta) {
            c(abeta)
        }
        else {
            dm <- dmat(pars, dispersion_block = TRUE, fit = fit, only_beta = only_beta)
            um <- umat(pars, dispersion_block = TRUE, fit = fit, only_beta = only_beta)
            apsi <- - sum(inv_jm * t(dm)) - sum(sandm * t(um)) / 2
            c(abeta, apsi)
        }
    }

    control <- do.call("brquasiControl", control)

    is_correction <- control$type == "eRBM"
    adjustment_function <- switch(control$type,
                                  "iRBM" = RBM_adjustment,
                                  "eRBM" = RBM_adjustment,
                                  "M" = function(pars, ...) 0,
                                  "MPQL_trace" = NA)

    ## Some useful quantities
    is_QL <- control$type == "M"

    ## Ensure x is a matrix, extract variable names, observation
    ## names, nobs, nvars, and initialize weights and offsets if
    ## needed
    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    converged <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)

    disp_factor <- switch(control$disp_factor,
                          "n-p"= nobs - nvars,
                          "n" = nobs)

    EMPTY <- nvars == 0
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missing_offset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }

    if (!isTRUE(family$family %in% c("quasi", "quasibinomial", "quasipoisson", "wedderburn"))) {
        stop("`brquasiFit` does not support families other than `quasi`, `quasipoisson` and `quasibinomial`.")
    }

    ## Enrich family
    family <- enrichwith::enrich(family, with = c("d1afun", "d2afun", "d3afun", "d1variance", "d2variance"))
    ## Enrich the link object with d2mu.deta and update family object
    linkglm <- make.link(family$link)
    linkglm <- enrichwith::enrich(linkglm, with = c("d2mu.deta", "d3mu.deta"))
    ## Put everything into the family object
    family[names(linkglm)] <- linkglm


    ## Extract functions from the enriched family object
    variance <- family$variance
    d1variance <- family$d1variance
    d2variance <- family$d2variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    if (!is.function(variance) || !is.function(linkinv)) {
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    }
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    ## If the family is custom then d2mu.deta cannot survive when
    ## passing throguh current family functions. But mu.eta does; so
    ## we compute d2mu.deta numerically; this allows also generality,
    ## as the users can then keep their custom link implementations
    ## unaltered. Issue is scalability, due to the need of evaluating
    ## n numerical derivatives

    if (is.null(family$d2mu.deta)) {
        d2mu.deta <- function(eta) {
            numDeriv::grad(mu.eta, eta)
        }
        d3mu.deta <- function(eta) {
            numDeriv::grad(d2mu.deta, eta)
        }
    }
    else {
        d2mu.deta <- family$d2mu.deta
        d3mu.deta <- family$d3mu.deta
    }
    ## d1afun <- family$d1afun
    ## d2afun <- family$d2afun
    ## d3afun <- family$d3afun
    simulate <- family$simulate

    ## Check for invalid etas and mus
    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu, function(mu) TRUE)

    ## FIXME: mustart and etastart set to NULL by default
    mustart <- NULL
    etastart <- NULL

    ## Initialize as prescribed in family
    eval(family$initialize)

    ## If there are no covariates in the model then evaluate only the offset
    if (EMPTY) {
        etas <- rep.int(0, nobs) + offset
        if (!valid_eta(etas))
            stop("invalid linear predictor values in empty model", call. = FALSE)
        mus <- linkinv(etas)
        if (!valid_mu(mus))
            stop("invalid fitted means in empty model", call. = FALSE)
        ## deviance <- sum(dev.resids(y, mus, weights))
        working_weights <- ((weights * mu.eta(etas)^2)/variance(mus))^0.5
        residuals <- (y - mus)/mu.eta(etas)
        keep <- rep(TRUE, length(residuals))
        boundary <- converged <- TRUE
        betas_all <- numeric()
        rank <- 0
        iter <- 0L
    }
    else {
        boundary <- converged <- FALSE
        ## Detect aliasing
        qrx <- qr(x)
        rank <- qrx$rank
        is_full_rank <- rank == nvars

        if (!singular.ok && !is_full_rank) {
            stop("singular fit encountered")
        }
        if (!isTRUE(is_full_rank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            X_all <- x
            x <- x[, -aliased]
            nvars_all <- nvars
            nvars <- ncol(x)
            betas_names_all <- betas_names
            betas_names <- betas_names[-aliased]
        }
        else {
            nvars_all <- nvars
            betas_names_all <- betas_names
        }
        betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
        keep <- weights > 0
        nkeep <- sum(keep)
        df_residual <- nkeep - rank
        ## Handle starting values
        ## If start is NULL then start at the ML estimator else use start
        if (is.null(start)) {
            ## Adjust counts if binomial or Poisson in order to avoid infinite estimates
            adj <- control$response_adjustment
            if (is.null(adj)) {
                adj <- nvars/nobs
            }
            if (family$family == "quasibinomial") {
                weights.adj <- weights + (!(is_correction)) * adj
                y.adj <- (weights * y + (!(is_correction)) * 0.5 * adj)/weights.adj
            }
            else {
                weights.adj <- weights
                y.adj <- y + if (family$family == "quasipoisson") (!(is_correction)) * 0.5 * adj else 0
            }
            ## ML fit to get starting values
            warn <- getOption("warn")
            ## Get startng values and kill warnings whilst doing that
            options(warn = -1)
            tempFit <- glm.fit(x = x, y = y.adj, weights = weights.adj,
                               etastart = etastart, mustart = mustart,
                               offset = offset, family = family,
                               control = list(epsilon = control$epsilon,
                                              maxit = 10000, trace = FALSE),
                               intercept = intercept)

            ## Set warn to its original value
            options(warn = warn)
            betas <- coef(tempFit)
            names(betas) <- betas_names

            dispersion <- sum((tempFit$weights * tempFit$residuals^2)[tempFit$weights >  0])/tempFit$df.residual
        }
        else {
            if (length(start) != nvars_all + 1) {
                stop(paste(paste(gettextf("length of 'start' should be equal to %d and correspond to initial betas for %s", nvars_all, paste(deparse(betas_names_all), collapse = ", "), "or", gettextf("to %d and also include a starting value for the transformed dispersion", nvars_all + 1)))), domain = NA_real_)
            }
            betas_all <- start[seq.int(nvars_all)]
            names(betas_all) <- betas_names_all
            if (!isTRUE(is_full_rank)) {
                betas_all[aliased] <- NA_real_
                betas <- betas_all[-aliased]
            }
            else {
                betas <- betas_all
            }
            dispersion <- start[nvars_all + 1]
        }

        adjusted_grad_all <- rep(NA_real_, nvars_all + 1)
        names(adjusted_grad_all) <- c(betas_names_all, "dispersion")


        compute_step <- function(pars, fit = NULL, only_beta = FALSE, lambda = 1e-06) {
            if (is.null(fit)) {
                fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
            }
            hess <- hessian(pars, fit = fit, only_beta = only_beta, lambda = lambda)
            grad <- gradient(pars, fit = fit, contributions = FALSE, only_beta = only_beta)
            if (!is_QL) {
                adj <- RBM_adjustment(pars, fit = fit, only_beta = only_beta)
                grad <- grad + adj
            }
            list(step = drop(solve(hess) %*% grad), grad = grad, fit = fit)
        }


        ee <- function(pars, only_beta = FALSE) {
            quantities <- try(key_quantities(pars, y = y, level = 2, qr = FALSE),
                              silent = TRUE)
            grad <- try(gradient(pars, fit = quantities, contributions = FALSE, only_beta = only_beta),
                        silent = TRUE)
            if (inherits(quantities, "try-error") || inherits(grad, "try-error")) {
                return(rep(NA_real_, nvars + !only_beta))
            }
            if (!is_QL) {
                adj <- try(RBM_adjustment(pars, fit = quantities, only_beta = only_beta),
                           silent = TRUE)
                if (inherits(adj, "try-error")) {
                    return(rep(NA_real_, nvars + !only_beta))
                }
                grad <- grad + adj
            }
            if (any(is.na(grad))) {
            }
            else {
                grad
            }
        }

        penalty <- function(pars, fit = NULL, only_beta = TRUE) {
            if (!only_beta) {
                stop("A quasi-likelihood for both beta and dispersion does not exist")
            }
            if (is.null(fit)) {
                fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
            }
            inv_jm <- chol2inv(chol(jmat(pars, fit = fit, only_beta = only_beta)))
            emat <- emat(pars, fit = fit, only_beta = only_beta)
            - sum(inv_jm * t(emat))/2
        }

        ## Allow for the use of generic non-linear equation solvers through a control argument

        ############

        slowit <- if (is_correction) 1 else control$slowit
        only_beta <- control$only_beta
        lam <- if (is_correction) 0 else control$lambda
        maxit <- if (is_correction) 1 else control$maxit
        max_step_factor <- if (is_correction) 1 else control$max_step_factor

        epsilon <- control$epsilon

        ## Evaluate at the starting values
        theta <- c(betas, dispersion)

        names(theta) <- c(betas_names, "dispersion")

        step0_l1 <- Inf
        failed <- FALSE
        ## Get step at the starting values
        quantities <- key_quantities(theta, y = y, level = 2, qr = FALSE)
        step_components <- try(compute_step(theta, fit = quantities, only_beta = only_beta, lambda = lam),
                               silent = TRUE)
        if (inherits(step_components, "try-error")) {
            failed <- TRUE
            step_l1 <- NA
        }
        else {
            step <- step_components$step
            step_l1 <- sum(abs(step))
        }
        for (iter in seq.int(maxit)) {
            if (failed) {
                break
            }
            test_step <- TRUE
            step_factor <- 0
            step0_l1 <- step_l1
            while (test_step & step_factor < max_step_factor) {
                if (only_beta) {
                    betas <- betas - slowit * step / 2^step_factor
                    theta <- c(betas, dispersion)
                }
                else {
                    theta <- theta - slowit * step / 2^step_factor
                }
                quantities <- key_quantities(theta, y = y, level = 2, qr = FALSE)
                step_components <- try(compute_step(theta, fit = quantities, only_beta = only_beta, lambda = lam),
                                       silent = TRUE)
                if (inherits(step_components, "try-error")) {
                    failed <- TRUE
                    break
                }
                step <- step_components$step
                grad <- step_components$grad
                step_l1 <- sum(abs(step))
                test_step <- step_l1 > step0_l1
                if (control$trace) {
                    cat("iter:", iter, "step", step_factor, "step length:", step_l1, sep = " | ",  "\n")
                }
                step_factor <- step_factor + 1
            }
            if (step_l1 < epsilon) break
        }

        if (failed) {
            grad <- rep(NA_real_, nvars + !only_beta)
        }
        else {
            grad <- step_components$grad
        }

        if (only_beta) {
            ## re-estimate dispersion
            dispersion <- theta[nvars + 1] <- sum(weights * (y - quantities$mus)^2 / quantities$varmus) / df_residual
        }


        names(theta) <- c(betas_names, "dispersion")

        adjusted_grad_all[betas_names] <- grad[betas_names]
        adjusted_grad_all["dispersion"] <- grad[nvars + 1]
        betas_all[betas_names] <- theta[betas_names]
        dispersion <- theta["dispersion"]

        ## Convergence analysis
        if ((failed || iter >= control$maxit) & !(is_correction)) {
            warning("brquasiFit: algorithm did not converge", call. = FALSE)
            converged <- FALSE
        }
        else {
            converged <- TRUE
        }

        ## QR decomposition and fitted values are at the final value
        ## for the coefficients
        ## QR decomposition for cov.unscaled
        if (!isTRUE(is_full_rank)) {
            x <- X_all
            betas <- betas_all
            betas[is.na(betas)] <- 0
            nvars <- nvars_all
        }

        qr.Wx <- try(qr(sqrt(quantities$working_weights) * x), silent = TRUE)
        if (inherits(qr.Wx, "try-error")) {
            qr.Wx <- NULL
        }

        mus <- quantities$mus
        etas <- quantities$etas
        ## Residuals
        residuals <- with(quantities, (y - mus)/d1mus)
        working_weights <- quantities$working_weights

        eps <- 10 * .Machine$double.eps
        if (family$family == "quasibinomial") {
            if (any(mus > 1 - eps) || any(mus < eps)) {
                warning("brquasiFit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
                boundary <- TRUE
            }
        }
        if (family$family == "quasipoisson") {
            if (any(mus < eps)) {
                warning("brquasiFit: fitted rates numerically 0 occurred", call. = FALSE)
                boundary <- TRUE
            }
        }
        if (df_residual == 0) {
            dispersion <- NA_real_
        }
    }

    ## Working weights
    wt <- rep.int(0, nobs)
    wt[keep] <- working_weights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <- names(weights) <- names(y) <- ynames
    ## For the null deviance:
    ##
    ## If there is an intercept but not an offset then the ML fitted
    ## value is the weighted average and is calculated easily below if
    ## ML is used
    ##
    control0 <- control
    control0$maxit <- 1000

    if (intercept & missing_offset) {
        nullFit <- brquasiFit(x = x[, "(Intercept)", drop = FALSE], y = y, weights = weights,
                              offset = rep(0, nobs), family = family, intercept = TRUE,
                              control = control0[c("epsilon", "maxit", "type", "slowit", "lambda", "only_beta")],
                              start = c(linkfun(mean(y)), 1))
        ## FIX: Starting values above are hard-coded. Change in future versions
        nullmus <- nullFit$fitted
    }
    ## If there is an offset but not an intercept then the fitted
    ## value is the inverse link evaluated at the offset
    ##
    ## If there is neither an offset nor an intercept then the fitted
    ## values is the inverse link at zero (and hence covered by
    ## linkinv(offset) because offset is zero
    if (!intercept) {
        nullmus <- linkinv(offset)
    }
    ## If there is an intercept and an offset then, for calculating
    ## the null deviance glm will make a call to the fitter to fit the
    ## glm with intercept and the offset
    if (intercept & !missing_offset) {
        nullmus <- mus
        ## doen't really matter what nullmus is set to. glm will make
        ## a new call to brquasiFit and use the deviance from that call
        ## as null
    }
    nulldev <- sum(dev.resids(y, nullmus, weights))
    nulldf <- nkeep - as.integer(intercept)
    deviance <- sum(dev.resids(y, mus, weights))
    aic.model <- aic(y, n, mus, weights, deviance) + 2 * rank

    list(coefficients = betas_all,
         residuals = residuals,
         fitted.values = mus,
         ## TODO: see effects?
         ## effects = if (!EMPTY) effects,
         R = if (!EMPTY && !is.null(qr.Wx)) qr.R(qr.Wx) else NULL,
         rank = rank,
         qr = if (!EMPTY) structure(qr.Wx[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
         family = family,
         linear.predictors = etas,
         deviance = deviance,
         aic = aic.model,
         null.deviance = nulldev,
         iter = iter,
         weights = wt,
         prior.weights = weights,
         df.residual = df_residual,
         df.null = nulldf,
         y = y,
         converged = converged,
         boundary = boundary,
         dispersion = dispersion,
         grad = adjusted_grad_all,
         transformation = control$transformation,
         type = control$type,
         class = "brquasiFit")
}

#' @export
coef.brquasiFit <- function(object, model = c("mean", "full", "dispersion"), ...) {
    model <- match.arg(model)
    switch(model,
           "mean" = {
               object$coefficients
           },
           "dispersion" = {
               object$dispersion
           },
           "full" = {
               c(object$coefficients, object$dispersion)
           })
}

#' `summary()` method for `brquasiFit()` objects
#'
#' @inheritParams stats::summary.glm
#'
#' @details The interface of the summary method for
#'     `brquasiFit()` objects is identical to that of
#'     `glm()` objects. The summary method for
#'     `brquasiFit()` objects computes the p-values of the
#'     individual Wald statistics based on the standard normal
#'     distribution, unless the family is Gaussian, in which case a t
#'     distribution with appropriate degrees of freedom is used.
#'
#' @seealso `summary.glm()` and `glm()`
#'
#' @examples
#' ## For examples see `examples(brquasiFit)`
#'
#' @export
summary.brquasiFit <- function(object, dispersion = NULL,
                            correlation = FALSE, symbolic.cor = FALSE,
                             ...) {
    if (is.null(dispersion)) {
        if (object$family$family == "Gaussian") {
            dispersion <- NULL
        }
        else {
            dispersion <- object$dispersion
        }
    }
    summary.glm(object, dispersion = dispersion,
                correlation = correlation,
                symbolic.cor = symbolic.cor, ...)
}

#' Method for computing confidence intervals for one or more
#' regression parameters in a \code{\link{brquasi}} object
#'
#' @inheritParams stats::confint
#'
#' @export
confint.brquasiFit <- function(object, parm, level = 0.95, ...) {
    confint.default(object, parm, level, ...)
}


DD <- function(expr,name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr,name)
    else DD(D(expr, name), name, order - 1)
}
