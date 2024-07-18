#### mc_evol (documentation) ####
#' Monte Carlo sampling for (average) evolvability measures
#'
#' These functions generate a Monte Carlo sample for
#' evolvability measures, whose average is an estimate of
#' the average evolvability measure
#'
#' All these functions are wrappers of either
#' \code{qfratio::\link[qfratio]{rqfr}()} or
#' \code{qfratio::\link[qfratio]{rqfmr}()} for random number generation of
#' simple and multiple ratios of quadratic forms in normal variables.  See
#' documentations of those functions for technical details.
#'
#' When \code{G} is singular and the distribution of the selection gradients
#' is not within its range, the conditional evolvability and autonomy are 0
#' with probability 1 (see Watanabe, 2023).  Accordingly, \code{mc_cevo()} and
#' \code{mc_auto()} simply yield a vector of 1's in this case.
#'
#' @inheritParams avr_evol
#'
#' @param nit
#'   Number of iteration
#' @param mu,Sigma
#'   Mean vector and covariance matrix, respectively, of selection gradients;
#'   default zero vector and identity matrix.  Passed to
#'   \code{qfratio::rqfr()} or \code{qfratio::rqfmr()}.
#' @param ...
#'   Additional arguments are passed to \code{qfratio::rqfr()} or
#'   \code{qfratio::rqfmr()}.  The only useful arguments will be \code{mu} and
#'   \code{Sigma} for the mean vector and covariance matrix, respectively,
#'   of selection gradients.
#'
#' @return
#' Numeric vector of length \code{nit}
#'
#' @references
#' Cheverud, J. M. (1996) Quantitative genetic analysis of cranial morphology
#'   in the cotton-top (*Saguinus oedipus*) and saddle-back(*S. fuscicollis*)
#'   tamarines. *Journal of Evolutionary Biology*, **9**, 5--42.
#'   \doi{10.1046/j.1420-9101.1996.9010005.x}.
#'
#' Cheverud, J. M. and Marroig, G. (2007) Comparing covariance matrices: random
#'   skewers method compared to the common principal components model.
#'   *Genetics and Molecular Biology*, **30**, 461--469.
#'   \doi{10.1590/S1415-47572007000300027}.
#'
#' Marroig, G., Shirai, L. T., Porto A., de Oliveira, F. B. and De Conto, V.
#'   (2009) The evolution of modularity in the mammalian skull II:
#'   evolutionary consequences. *Evolutionary Biology*, **36**, 136--148.
#'   \doi{10.1007/s11692-009-9051-1}.
#'
#' Watanabe, J. (2023) Exact expressions and numerical evaluation of average
#'   evolvability measures for characterizing and comparing **G** matrices.
#'   *Journal of Mathematical Biology*, **86**, 95.
#'   \doi{10.1007/s00285-023-01930-8}.
#'
#' @seealso
#' \code{\link[qfratio]{rqfr}}
#'
#' @name mc_evol
#'
#' @examples
#' ## Simple covariance matrices
#' nv <- 4
#' G1 <- diag(nv:1)
#' G2 <- diag(sqrt(1:nv))
#' nit <- 1000
#'
#' ## Monte Carlo sample of conditional evolvability
#' ## under spherical distribution of beta,
#' ## and its mean as an estimate of average conditional evolvability
#' ## plus its 95% CI
#' mcsample_cevo <- mc_cevo(nit, G1)
#' mean(mcsample_cevo)
#' mean(mcsample_cevo) + sd(mcsample_cevo) / sqrt(nit) *
#'     qt(c(0.025, 0.975), nit - 1)
#'
#' ## Monte Carlo estimate of average response difference
#' mean(mc_rdif(nit, G1, G2))
#'
#' ## Monte Carlo estimate of average response correlation,
#' ## aka random skewers correlation
#' mean(mc_rcor(nit, G1, G2))
#'
#' ## Advanced: Monte Carlo estimation of average evolvability under
#' ## non-spherical distribution of selection gradient
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#' mean(mc_evol(nit, G1, mu = mu, Sigma = Sigma))
#'
#' ## For singular G, these yield vectors of 1's/0's
#' ## (unless appropriate Sigma is provided)
#' Gsing <- diag(nv:1 - 1)
#' mc_cevo(5, Gsing)
#' mc_auto(5, Gsing)
#' mc_inte(5, Gsing)
#'
NULL

#### mc_evol ####
#' Monte Carlo sampling for (average) evolvability
#'
#' \code{mc_evol()}: (Average) evolvability
#'
#' @rdname mc_evol
#' @export
mc_evol <- function(nit, G, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    return(qfratio::rqfr(nit = nit, A = G, B = diag(dim(G)[1]),
                         p = 1, q = 1, ...))
}

#### mc_cevo ####
#' Monte Carlo sampling for (average) conditional evolvability
#'
#' \code{mc_cevo()}: (Average) conditional evolvability
#'
#' @rdname mc_evol
#' @export
mc_cevo <- function(nit, G, mu = rep.int(0, nvar),
                    Sigma = diag(nvar), tol_qr = 1e-7, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    qrG <- qr(G, tol = tol_qr)
    if (qrG$rank < nvar) {
        proj_RG <- tcrossprod(qr.Q(qrG)[, seq_len(qrG$rank)])
        mu_ <- proj_RG %*% c(mu)
        Sigma_ <- tcrossprod(proj_RG, tcrossprod(proj_RG, Sigma))
        beta_in_RG <- qfratio:::iseq(mu, mu_) && qfratio:::iseq(Sigma, Sigma_)
        if (beta_in_RG) {
            Gi <- MASS::ginv(G)
            return(qfratio::rqfr(nit = nit, A = diag(nvar), B = Gi,
                                 p = 1, q = 1, mu = mu, Sigma = Sigma, ...))
        } else {
            message("Singular covariance matrix; ",
                    "conditional evolvability is 0 with probability 1")
            return(rep.int(0, nit))
        }
    } else {
        Gi <- chol2inv(chol(G))
        return(qfratio::rqfr(nit = nit, A = diag(nvar), B = Gi,
                             p = 1, q = 1, mu = mu, Sigma = Sigma, ...))
    }
}

#### mc_resp ####
#' Monte Carlo sampling for (average) respondability
#'
#' \code{mc_resp()}: (Average) respondability
#'
#' @rdname mc_evol
#' @export
mc_resp <- function(nit, G, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    return(qfratio::rqfr(nit = nit, A = crossprod(G), B = diag(nvar),
                         p = 1/2, q = 1/2, ...))
}

#### mc_flex ####
#' Monte Carlo sampling for (average) flexibility
#'
#' \code{mc_flex()}: (Average) flexibility
#'
#' @rdname mc_evol
#' @export
mc_flex <- function(nit, G, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    return(qfratio::rqfmr(nit = nit, A = G, B = crossprod(G), D = diag(nvar),
                          p = 1, q = 1/2, r = 1/2, ...))
}

#### mc_auto ####
#' Monte Carlo sampling for (average) autonomy
#'
#' \code{mc_auto()}: (Average) autonomy
#'
#' @rdname mc_evol
#' @export
mc_auto <- function(nit, G, mu = rep.int(0, nvar),
                    Sigma = diag(nvar), tol_qr = 1e-7, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    qrG <- qr(G, tol = tol_qr)
    if (qrG$rank < nvar) {
        proj_RG <- tcrossprod(qr.Q(qrG)[, seq_len(qrG$rank)])
        mu_ <- proj_RG %*% c(mu)
        Sigma_ <- tcrossprod(proj_RG, tcrossprod(proj_RG, Sigma))
        beta_in_RG <- qfratio:::iseq(mu, mu_) && qfratio:::iseq(Sigma, Sigma_)
        if (beta_in_RG) {
            Gi <- MASS::ginv(G)
            return(qfratio::rqfmr(nit = nit, A = diag(nvar), B = G, D = Gi,
                                  p = 2, q = 1, r = 1,
                                  mu = mu, Sigma = Sigma, ...))
        } else {
            message("Singular covariance matrix; ",
                    "autonomy is 0 with probability 1")
            return(rep.int(0, nit))
        }
    } else {
        Gi <- chol2inv(chol(G))
        return(qfratio::rqfmr(nit = nit, A = diag(nvar), B = G, D = Gi,
                              p = 2, q = 1, r = 1,
                              mu = mu, Sigma = Sigma, ...))
    }
}

#### mc_inte ####
#' Monte Carlo sampling for (average) integration
#'
#' \code{mc_inte()}: (Average) integration
#'
#' @rdname mc_evol
#' @export
mc_inte <- function(nit, G, ...) {
    1 - mc_auto(nit = nit, G = G, ...)
}

#### mc_cons ####
#' Monte Carlo sampling for (average) constraints
#'
#' \code{mc_cons()}: (Average) constraints
#'
#' @rdname mc_evol
#' @export
mc_cons <- function(nit, G, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    Lsq <- eigen(G, symmetric = TRUE, only.values = TRUE)$values ^ 2
    nvar <- length(Lsq)
    Gsq1 <- diag(c(Lsq[1], rep.int(0, nvar - 1)))
    Gsq <- diag(Lsq)
    qfratio::rqfr(nit = nit, A = Gsq1, B = Gsq, p = 1/2, q = 1/2, ...)
}

#### mc_rdif ####
#' Monte Carlo sampling for (average) response difference
#'
#' \code{mc_rdif()}: (Average) response difference
#'
#' @rdname mc_evol
#' @export
mc_rdif <- function(nit, G1, G2, ...) {
    stopifnot(
        "G1 and G2 should be symmetric" = isSymmetric(G1) && isSymmetric(G2),
        "G1 and G2 should have the same dimension" = all(dim(G1) == dim(G2))
    )
    nvar <- dim(G1)[1]
    return(qfratio::rqfr(nit = nit, A = crossprod(G1 - G2), B = diag(nvar),
                         p = 1/2, q = 1/2, ...))
}

#### mc_rcor ####
#' Monte Carlo sampling for (average) response correlation
#'
#' \code{mc_rcor()}: (Average) response correlation
#' (aka random skewers correlation)
#'
#' @rdname mc_evol
#' @export
mc_rcor <- function(nit, G1, G2, ...) {
    stopifnot(
        "G1 and G2 should be symmetric" = isSymmetric(G1) && isSymmetric(G2),
        "G1 and G2 should have the same dimension" = all(dim(G1) == dim(G2))
    )
    return(qfratio::rqfmr(nit = nit, A = crossprod(G1, G2),
                          B = crossprod(G1), D = crossprod(G2),
                          p = 1, q = 1/2, r = 1/2, ...))
}
