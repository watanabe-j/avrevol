#### avr_evol (documentation) ####
#' Series evaluation of average evolvability measures
#'
#' These functions evaluate series expressions of average evolvability measures
#' as their partial sums.
#'
#' All these functions are wrappers of either
#' \code{qfratio::\link[qfratio]{qfrm}()} or
#' \code{qfratio::\link[qfratio]{qfmrm}()} for moments of simple and multiple
#' ratios of quadratic forms in normal variables.  See documentations of
#' those functions for technical details.
#'
#' When \eqn{\mathbf{G}}{G} is singular, average conditional evolvability
#' and average autonomy can be nonzero only when the distribution of
#' \eqn{\bm{\beta}}{\beta} is within the range of \eqn{\mathbf{G}}{G}.  In
#' \code{avr_cevo()} and \code{avr_auto()}, this condition is assessed by
#' the QR decomposition of \eqn{\mathbf{G}}{G} (with pivoting).  When
#' \eqn{\mathbf{G}}{G} is not of full-rank, \eqn{\bm{\mu}}{\mu} and
#' \eqn{\bm{\Sigma}}{\Sigma} are projected onto the range of \eqn{\mathbf{G}}{G}
#' using the Q part of the QR decomposition.  If these are equal to the
#' original ones, then the arguments are passed to evaluation of the moment
#' using quadratic  forms; otherwise, these functions return 0 (formatted as
#' an \code{qfrm} object).
#'
#' @param G,G1,G2
#'   Target covariance matrix/matrices.  Assumed validly structured
#'   (symmetric, nonnegative definite), although symmetry is checked.
#' @param m
#'   Order of evaluation \eqn{m} at which partial sum is truncated
#'   (\eqn{M} in, e.g., Hillier et al., 2014).  Adjust this value
#'   depending on desired accuracy.
#' @param mu,Sigma
#'   Mean vector and covariance matrix, respectively, of selection gradients;
#'   default zero vector and identity matrix with dimension \code{nvar}
#'   deduced from that of \code{G} (or \code{G1}).  Passed to
#'   \code{qfratio::qfrm()} or \code{qfratio::qfmrm()}.
#' @param cpp_method
#'   Option to specify \code{C++} algorithm to avoid numerical
#'   overflow/underflow in \code{qfratio::qfrm()} and
#'   \code{qfratio::qfmrm()}.  Default \code{"coef_wise"}
#'   is typically the most robust and modestly fast option.  See documentation
#'   of \code{qfratio::\link[qfratio]{qfrm}()} for details.
#' @param check_convergence
#'   Option to specify how numerical convergence is checked.  Default for
#'   \code{avr_auto()}, \code{avr_cons()}, and \code{avr_rcor()} is a strict one
#'   (\code{"strict_relative"}) because no error bound is available for these
#'   measures.  Other options are: \code{"relative"}, which is far less strict,
#'   and \code{"none"} (or \code{FALSE}) for no checking.
#' @param tol,tol_qr,tol_ginv
#'   Tolerance to determine singularity of \code{G}.  By default, \code{tol}
#'   is copied into \code{tol_qr} and \code{tol_ginv}, which in turn are
#'   passed to \code{\link[base]{qr}()} and \code{\link[MASS]{ginv}()},
#'   respectively.  These can be specified separately if desired
#'   (mainly for debugging purposes).
#' @param ...
#'   Additional arguments are passed to \code{qfratio::qfrm()} or
#'   \code{qfratio::qfmrm()}.  Notable arguments involve \code{mu} and
#'   \code{Sigma} for the mean vector and covariance matrix, respectively,
#'   of selection gradients.
#'
#' @return
#' \code{qfrm} object defined by the package \code{qfratio}, which is a list
#' including \code{$statistic} (the partial sum) and \code{$error_bound}.  That
#' package also defines \code{print} and \code{plot} methods.
#'
#' @references
#' Cheverud, J. M. (1996) Quantitative genetic analysis of cranial morphology
#'   in the cotton-top (*Saguinus oedipus*) and saddle-back (*S. fuscicollis*)
#'   tamarines. *Journal of Evolutionary Biology*, **9**, 5--42.
#'   \doi{10.1046/j.1420-9101.1996.9010005.x}.
#'
#' Hansen, T. F. and Houle, D. (2008) Measuring and comparing evolvability and
#'   constraint in multivariate characters. *Journal of Evolutionary Biology*,
#'   **21**, 1201--1219.
#'   \doi{10.1111/j.1420-9101.2008.01573.x}.
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
#' \code{\link[qfratio]{qfrm}} and \code{\link[qfratio]{qfmrm}}
#'
#' @name avr_evol
#'
#' @examples
#' ## Simple covariance matrices
#' nv <- 4
#' G1 <- diag(nv:1)
#' G2 <- diag(sqrt(1:nv))
#' nit <- 1000
#'
#' ## Average conditional evolvability using series expression
#' (res_cevo <- avr_cevo(G1))
#'
#' ## Inspect its plot to check for convergence
#' plot(res_cevo)
#'
#' ## Average response difference using series expression
#' (res_rdif <- avr_rdif(G1, G2))
#' plot(res_rdif)
#'
#' ## Average response correlation using series expression
#' ## whose Monte Carlo estimate is known as "random skewers" correlation
#' (res_rcor <- avr_rcor(G1, G2, m = 100)) # throws warning for non-convergence
#'
#' (res_rcor <- avr_rcor(G1, G2, m = 500)) # good with larger m
#' plot(res_rcor)
#'
#' ## Larger problem
#' nv_l <- 20
#' G_l <- diag((nv_l:1)^2)
#' (res_cevo_l <- avr_cevo(G_l, m = 100)) # throws warning for non-convergence
#' plot(res_cevo_l) # Note the growing partial sum
#'
#' ## Use larger m to evaluate further
#' (res_cevo_l <- avr_cevo(G_l, m = 1000))
#' plot(res_cevo_l)
#'
#' ## Advanced: Average evolvability under
#' ## non-spherical distribution of selection gradient
#' ## (the same works for other evolvability measures as well)
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#' (res_evol_nsph <- avr_evol(G1, mu = mu, Sigma = Sigma))
#' plot(res_evol_nsph)
#'
#' ## For singular G, these yield 0/1 ...
#' Gsing <- diag(nv:1 - 1)
#' avr_cevo(Gsing)
#' avr_auto(Gsing)
#' avr_inte(Gsing)
#'
#' ## ... unless distribution of selection gradients is restricted to G's range
#' Sigmasing <- diag(c(rep.int(1, nv - 1), 0))
#' avr_cevo(Gsing, Sigma = Sigmasing)
#' avr_auto(Gsing, Sigma = Sigmasing)
#' avr_inte(Gsing, Sigma = Sigmasing)
#'
NULL

#### avr_evol ####
#' Series evaluation of average evolvability
#'
#' \code{avr_evol()}: Average evolvability
#'
#' @rdname avr_evol
#' @export
avr_evol <- function(G, m = 100L, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    qfratio::qfrm(A = G, B = diag(dim(G)[1]), p = 1, q = 1, m = m, ...)
}

#### avr_cevo ####
#' Series evaluation of average conditional evolvability
#'
#' \code{avr_cevo()}: Average conditional evolvability
#'
#' @rdname avr_evol
#' @export
avr_cevo <- function(G, m = 100L, mu = rep.int(0, nvar),
                     Sigma = diag(nvar), tol = sqrt(.Machine$double.eps),
                     tol_qr = tol, tol_ginv = tol, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    qrG <- qr(G, tol = tol_qr)
    if (qrG$rank < nvar) {
        proj_RG <- tcrossprod(qr.Q(qrG)[, seq_len(qrG$rank)])
        mu_ <- proj_RG %*% c(mu)
        Sigma_ <- tcrossprod(proj_RG, tcrossprod(proj_RG, Sigma))
        beta_in_RG <- qfratio:::iseq(mu, mu_) && qfratio:::iseq(Sigma, Sigma_)
        if (beta_in_RG) {
            Gi <- MASS::ginv(G, tol = tol_ginv)
            return(qfratio::qfrm(A = diag(nvar), B = Gi, p = 1, q = 1,
                                 m = m, mu = mu, Sigma = Sigma, ...))
        } else {
            message("Singular covariance matrix; ",
                    "average conditional evolvability is 0")
            return(qfratio:::new_qfrm(0, 0, 0, 0, exact = TRUE))
        }
    } else {
        Gi <- chol2inv(chol(G))
        return(qfratio::qfrm(A = diag(nvar), B = Gi, p = 1, q = 1, m = m,
                             mu = mu, Sigma = Sigma, ...))
    }
}

#### avr_resp ####
#' Series evaluation of average respondability
#'
#' \code{avr_resp()}: Average respondability
#'
#' @rdname avr_evol
#' @export
avr_resp <- function(G, m = 100L, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    Gsq <- crossprod(G)
    qfratio::qfrm(A = Gsq, B = diag(nvar), p = 1/2, q = 1/2, m = m, ...)
}

#### avr_flex ####
#' Series evaluation of average flexibility
#'
#' \code{avr_flex()}: Average flexibility
#'
#' @rdname avr_evol
#' @export
avr_flex <- function(G, m = 100L, ...) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    Gsq <- crossprod(G)
    qfratio::qfmrm(A = G, B = Gsq, D = diag(nvar),
                   p = 1, q = 1/2, r = 1/2, m = m, ...)
}

#### avr_auto ####
#' Series evaluation of average autonomy
#'
#' \code{avr_auto()}: Average autonomy
#'
#' @rdname avr_evol
#' @export
avr_auto <- function(G, m = 100L, mu = rep.int(0, nvar),
                     Sigma = diag(nvar),
                     cpp_method = c("coef_wise", "double", "long_double"),
                     check_convergence = "strict_relative",
                     tol = sqrt(.Machine$double.eps),
                     tol_qr = tol, tol_ginv = tol, ...) {
    cpp_method <- match.arg(cpp_method)
    stopifnot("G should be symmetric" = isSymmetric(G))
    nvar <- dim(G)[1]
    qrG <- qr(G, tol = tol_qr)
    if (qrG$rank < nvar) {
        proj_RG <- tcrossprod(qr.Q(qrG)[, seq_len(qrG$rank)])
        mu_ <- proj_RG %*% c(mu)
        Sigma_ <- tcrossprod(proj_RG, tcrossprod(proj_RG, Sigma))
        beta_in_RG <- qfratio:::iseq(mu, mu_) && qfratio:::iseq(Sigma, Sigma_)
        if (beta_in_RG) {
            Gi <- MASS::ginv(G, tol = tol_ginv)
            return(qfratio::qfmrm(A = diag(nvar), B = G, D = Gi,
                                  p = 2, q = 1, r = 1, m = m,
                                  mu = mu, Sigma = Sigma,
                                  cpp_method = cpp_method, ...))
        } else {
            message("Singular covariance matrix; average autonomy is 0")
            return(qfratio:::new_qfrm(0, 0, 0, 0, exact = TRUE))
        }
    } else {
        Gi <- chol2inv(chol(G))
        return(qfratio::qfmrm(A = diag(nvar), B = G, D = Gi,
                              p = 2, q = 1, r = 1, m = m,
                              mu = mu, Sigma = Sigma,
                              cpp_method = cpp_method,
                              check_convergence = check_convergence, ...))
    }
}

#### avr_inte ####
#' Series evaluation of average integration
#'
#' \code{avr_inte()}: Average integration
#'
#' @rdname avr_evol
#' @export
avr_inte <- function(G, ...) {
    auto <- avr_auto(G = G, ...)
    qfratio:::new_qfrm(statistic = 1 - auto$statistic,
                       terms = c(1 - auto$terms[1],
                                 diff(1 - cumsum(auto$terms))),
                       error_bound = -1 * auto$error_bound,
                       seq_error = -1 * auto$seq_error)
}

#### avr_cons ####
#' Series evaluation of average constraints
#'
#' \code{avr_cons()}: Average constraints
#'
#' @rdname avr_evol
#' @export
avr_cons <- function(G, m = 100L,
                     cpp_method = c("coef_wise", "double", "long_double"),
                     check_convergence = "strict_relative", ...) {
    cpp_method <- match.arg(cpp_method)
    stopifnot("G should be symmetric" = isSymmetric(G))
    Lsq <- eigen(G, symmetric = TRUE, only.values = TRUE)$values ^ 2
    nvar <- length(Lsq)
    Gsq1 <- diag(c(Lsq[1], rep.int(0, nvar - 1)), nrow = nvar, ncol = nvar)
    Gsq <- diag(Lsq, nrow = nvar, ncol = nvar)
    qfratio::qfrm(A = Gsq1, B = Gsq, p = 1/2, q = 1/2, m = m,
                  cpp_method = cpp_method,
                  check_convergence = check_convergence, ...)
}

#### avr_rdif ####
#' Series evaluation of average response difference
#'
#' \code{avr_rdif()}: Average response difference
#'
#' @rdname avr_evol
#' @export
avr_rdif <- function(G1, G2, m = 100L, ...) {
    stopifnot(
        "G1 and G2 should be symmetric" = isSymmetric(G1) && isSymmetric(G2),
        "G1 and G2 should have the same dimension" = all(dim(G1) == dim(G2))
    )
    nvar <- dim(G1)[1]
    Gdsq <- crossprod(G1 - G2)
    qfratio::qfrm(A = Gdsq, B = diag(nvar), p = 1/2, q = 1/2, m = m, ...)
}

#### avr_rcor ####
#' Series evaluation of average response correlation
#'
#' \code{avr_rcor()}: Average response correlation
#'
#' @rdname avr_evol
#' @export
avr_rcor <- function(G1, G2, m = 100L,
                     cpp_method = c("coef_wise", "double", "long_double"),
                     check_convergence = "strict_relative",
                     ...) {
    cpp_method <- match.arg(cpp_method)
    stopifnot(
        "G1 and G2 should be symmetric" = isSymmetric(G1) && isSymmetric(G2),
        "G1 and G2 should have the same dimension" = all(dim(G1) == dim(G2))
    )
    G12 <- crossprod(G1, G2)
    G1sq <- crossprod(G1)
    G2sq <- crossprod(G2)
    qfratio::qfmrm(A = G12, B = G1sq, D = G2sq, p = 1, q = 1/2, r = 1/2, m = m,
                   cpp_method = cpp_method,
                   check_convergence = check_convergence, ...)
}
