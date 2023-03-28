#### hh_evol (documentation) ####
#' Delta method approximation of average evolvability measures
#'
#' These functions calculate Hansen and Houle's (2008, 2009)
#' delta method approximations for average evolvability measures.
#'
#' The notation \eqn{\mathrm{Var}}{Var} in Hansen and Houle's (2008, 2009)
#' seems to represent average sum of squares from the average (divided by
#' \eqn{k}, the number of variables) and not the ordinary unbiased estimator
#' of variance for infinite populations (divided by \eqn{k - 1}).  The
#' package \code{evolvability} v 2.0.0 incorrectly uses the latter divisor.
#'
#' \code{svd(..., nu = 0L, nv = 0L)$d} is used to extract eigenvalues because
#' it is slightly faster than \code{eigen(..., symmetric = TRUE)$values}.
#' Note that the squared singular values of \code{A} is equal to the eigenvalues
#' of \code{crossprod(A)}.
#'
#' Evolvability measures not treated in Hansen and Houle (2008) (flexibility and
#' response correlation) have not been implemented (although this should be
#' doable).  All the expressions are for the spherical normal (or uniform on
#' hypersphyere) distribution of selection gradients.
#'
#' @inheritParams avr_evol
#'
#' @references
#' Hansen, T. F. and Houle, D. (2008) Measuring and comparing evolvability and
#'   constraint in multivariate characters. *Journal of Evolutionary Biology*,
#'   **21**, 1201--1219.
#'   \doi{10.1111/j.1420-9101.2008.01573.x}.
#'
#' Hansen, T. F. and Houle, D. (2009) Corrigendum \[of ``Measuring and comparing
#'   evolvability and constraint in multivariate characters''\].
#'   *Journal of Evolutionary Biology*,
#'   **22**, 913--915.
#'   \doi{10.1111/j.1420-9101.2009.01715.x}.
#'
#' @name hh_evol
#'
#' @examples
#' ## Simple covariance matrices
#' nv <- 4
#' G1 <- diag(nv:1)
#' G2 <- diag(sqrt(1:nv))
#'
#' ######
#' ## Hansen & Houle's (2008) delta method approximations:
#' 
#' ## - average conditional evolvability
#' hh_cevo(G1)
#'
#' ## - average respondability
#' hh_resp(G1)
#'
#' ## - average autonomy
#' hh_auto(G1)
#'
#' ## - average response difference
#' hh_rdif(G1, G2)
#'
NULL

#### hh_cevo ####
#' Average conditional evolvability using Hansen--Houle delta method
#'
#' \code{hh_cevo()}: Average conditional evolvability
#'
#' @rdname hh_evol
#' @export
hh_cevo <- function(G) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    I <- function(x) sum((x - mean(x)) ^ 2) / mean(x)^2 / length(x)
    L <- svd(G, nu = 0L, nv = 0L)$d
    p <- length(L)
    (1 + 2 * I(1 / L) / (p + 2)) / mean(1 / L)
}

#### hh_resp ####
#' Average respondability using Hansen--Houle delta method
#'
#' \code{hh_resp()}: Average respondability
#'
#' @rdname hh_evol
#' @export
hh_resp <- function(G) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    I <- function(x) sum((x - mean(x)) ^ 2) / mean(x)^2 / length(x)
    L <- svd(G, nu = 0L, nv = 0L)$d
    p <- length(L)
    sqrt(mean(L ^ 2)) * (1 - I(L ^ 2) / 4 / (p + 2))
}

#### hh_auto ####
#' Average autonomy using Hansen--Houle delta method
#'
#' \code{hh_auto()}: Average autonomy
#'
#' @rdname hh_evol
#' @export
hh_auto <- function(G) {
    stopifnot("G should be symmetric" = isSymmetric(G))
    I <- function(x) sum((x - mean(x)) ^ 2) / mean(x)^2 / length(x)
    H <- function(x) 1 / mean(1 / x)
    L <- svd(G, nu = 0L, nv = 0L)$d
    p <- length(L)
    H(L) / mean(L) *
        ( 1 + 2 * (I(L) + I(1 / L) - 1 + H(L) / mean(L) +
         2 * I(L) * I(1 / L) / (p + 2)) / (p + 2) )
}

#### hh_inte ####
#' Average integration using Hansen--Houle delta method
#'
#' \code{hh_inte()}: Average integration
#'
#' @rdname hh_evol
#' @export
hh_inte <- function(G) {
    1 - hh_auto(G)
}

#### hh_rdif ####
#' Average response difference using Hansen--Houle delta method
#'
#' \code{hh_rdif()}: Average response difference
#'
#' @rdname hh_evol
#' @export
hh_rdif <- function(G1, G2) {
    stopifnot(
        "G1 and G2 should be symmetric" = isSymmetric(G1) && isSymmetric(G2),
        "G1 and G2 should have the same dimension" = all(dim(G1) == dim(G2))
    )
    I <- function(x) sum((x - mean(x)) ^ 2) / mean(x)^2 / length(x)
    L12sq <- svd(G1 - G2, nu = 0L, nv = 0L)$d ^ 2
    p <- length(L12sq)
    sqrt(mean(L12sq)) * (1 - I(L12sq) / 4 / (p + 2))
}
