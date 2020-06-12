library("numDeriv")
data("lizards", package = "brglm")
source("~/Dropbox/brquasi/brquasi.R")
source("~/Dropbox/brquasi/brquasiControl.R")
source("~/Dropbox/brquasi/utils.R")

lizardsBR_mean <- glm(cbind(grahami, opalinus) ~ height + diameter +
                      light + time, family = quasibinomial("logit"), data = lizards,
                      method = "brquasi")


gradient1 <- function(pars, fit = NULL, contributions = FALSE, only_beta = FALSE) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    }
    with(fit, {
        bm <- weights * d1mus / varmus
        gm <- bm * (y - mus)
        score_components_beta <- precision * gm * x
        if (only_beta) {
                score_components <- score_components_beta
        }
        else {
            cm <- weights / varmus
            km <- cm * (y - mus)^2
            score_components_dispersion <- km - dispersion * (disp_factor)/nobs
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
jmat1 <- function(pars, fit = NULL, only_beta = FALSE) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    }
    with(fit, {
        bmus <- weights * d1mus / varmus
        cmus <- weights / varmus
        d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
        d1cmus <- - weights * d1mus * d1varmus / varmus^2
        emus <- (y - mus)
        qmus <- (bmus * d1mus - d1bmus * emus)
        gmus <- bmus * emus
        fmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2
        jbetabeta <- precision * crossprod(x * qmus, x)
        if (only_beta) {
            return(jbetabeta)
        }
        else {
            jbetaphi_components <- precision^2 * gmus * x
            jphibeta_components <- fmus * x
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


## The gradient of the set of estimating functions that gives rise
## to the pearson esitmator of the dispersion parameter
emat1 <- function(pars, fit = NULL, only_beta = FALSE) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    }
    with(fit, {
        bmus <- weights * d1mus / varmus
        cmus <- weights / varmus
        d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
        d1cmus <- - weights * d1mus * d1varmus / varmus^2
        emus <- (y - mus)
        qmus <- (bmus * d1mus - d1bmus * emus)
        gmus <- bmus * emus
        fmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2
        cmus <- weights / varmus
        kmus <- cmus * (y - mus)^2
        ebetabeta <- precision^2 * crossprod(x * gmus^2, x)
        if (only_beta) {
            return(jbetabeta)
        }
        else {
            ebetaphi <- precision * t(x * gmus * (kmus - dispersion)) %*% rep(1, nobs)
            ephiphi <- sum((kmus  - dispersion)^2)
            out <- rbind(cbind(ebetabeta, ebetaphi),
                         c(ebetaphi, ephiphi))
            return(out)
        }
    })
}


umat1 <- function(pars, dispersion_block = FALSE, fit = NULL, only_beta = FALSE){
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    }
    with(fit, {
        emus <- (y - mus)
        if (dispersion_block) {
            cmus <- weights / varmus
            d1cmus <- - weights * d1mus * d1varmus / varmus^2
            d2cmus <- - weights * (d2varmus * d1mus^2 + d1varmus * d2mus - 2 * d1mus^2 * d1varmus^2 / varmus) / varmus^2
            umus <- d2cmus * emus^2 - 4 * d1cmus * d1mus * emus - 2 * cmus * d2mus * emus + 2 * cmus * d1mus^2
            ubetabeta <- crossprod(x * umus, x)
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
            d2bmus <- weights * (d3mus / varmus - (3 * d1mus * d2mus * d1varmus + d1mus^3 * d2varmus) / varmus^2 + 2 * d1mus^3 * d1varmus^2 / varmus^3)
            qmus <- (bmus * d1mus - d1bmus * emus)
            pmus <- bmus * emus
            d1qmus <- (2 * d1bmus * d1mus + bmus * d2mus - d2bmus * emus)
            if (only_beta) {
                out <- lapply(seq.int(nvars), function(j) {
                    - precision * t(x) %*% diag(d1qmus * x[, j]) %*% x## crossprod(x * d1qmus * x[, j], x)
                })
            }
            else {
                out <- lapply(seq.int(nvars), function(j) {
                    ubetabeta <-  - precision * crossprod(x * d1qmus * x[, j], x)
                    ubetaphi <- precision^2 * crossprod(x * qmus, x)[j, ] 
                    uphiphi <- 2 * precision^3 * .colSums(pmus * x, nobs, nvars, TRUE)[j]
                    rbind(cbind(ubetabeta, ubetaphi),
                          c(ubetaphi, uphiphi))
                })
            }
        }
        return(out)
    })
}


dmat1 <- function(pars, dispersion_block = FALSE, fit = NULL, only_beta = FALSE) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = 2, qr = FALSE)
    }
    with(fit, {
        emus <- (y - mus)
                    cmus <- weights / varmus
            bmus <- weights * d1mus / varmus
            cmus <- weights / varmus
            d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
            d1cmus <- - weights * d1mus * d1varmus / varmus^2
            emus <- (y - mus)
            qmus <- (bmus * d1mus - d1bmus * emus)
            gmus <- bmus * emus
            fmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2
            cmus <- weights / varmus
            kmus <- cmus * (y - mus)^2
            d1cmus <- - weights * d1mus * d1varmus / varmus^2
            fmus <- 2 * cmus * d1mus * emus - d1cmus * emus^2

        if (dispersion_block) {                
            dbb <- - precision * crossprod(x * fmus * gmus, x)
            dbp <- - colSums(x * fmus * (kmus - dispersion))
            dpb <- - precision * colSums(gmus * x)
            dpp <- - sum(kmus - dispersion)            
            out <- rbind(cbind(dbb, dbp), c(dpb, dpp))
        }
        else {
            bmus <- weights * d1mus / varmus
            d1bmus <- weights * (d2mus / varmus - d1mus^2 * d1varmus / varmus^2)
            qmus <- (bmus * d1mus - d1bmus * emus) 
            gmus <- bmus * emus
            if (only_beta) {
                out <- lapply(seq.int(nvars), function(r) {
                    - crossprod(x * qmus  * x[, r], score_components)
                })
            }
            else {
                out <- lapply(seq.int(nvars), function(r) {
                    dbb <- - precision^2 * t(x) %*% diag(qmus * x[, r] * gmus) %*% x
                    dbp <- - precision * t(x) %*% diag(qmus * x[, r] * (kmus - dispersion)) %*% rep(1, nobs)
                    dpb <- - precision^3 * colSums(x * x[, r] * gmus^2)
                    dpp <- - precision^2 * sum(gmus * x[, r] * (kmus - dispersion))
                    rbind(cbind(dbb, dbp), c(dpb, dpp))
                })
            }
        }
        return(out)
    })
}



set.seed(111)
theta0 <- theta
theta0[1:6] <- rexp(6)
theta0[7] <- 2.3

disp_factor <- nobs
max(abs(gradient(theta0) - gradient1(theta0)))
max(abs(jmat(theta0) - jmat1(theta0)))
max(abs(emat(theta0) - emat1(theta0)))

max(abs(umat(theta0, dispersion_block = TRUE) - umat1(theta0, dispersion_block = TRUE)))
lapply(1:6, function(j) {
    max(abs(umat(theta0)[[j]] - umat1(theta0)[[j]]))
})

max(abs(dmat(theta0, dispersion_block = TRUE) - dmat1(theta0, dispersion_block = TRUE)))
lapply(1:6, function(j) {
    max(abs(dmat(theta0)[[j]] - dmat1(theta0)[[j]]))
})


