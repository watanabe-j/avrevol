## Documentation of the package
#' avrevol: Average Evolvability Measures
#'
#' This small package provides wrapper functions to evaluate average
#' evolvability measures in evolutionary quantitative genetics, based on the
#' implementation of recursion-based evaluation of moments of ratios of
#' quadratic forms in the package \code{qfratio}.
#'
#' The package \code{qfratio} provides functions to evaluate moments of
#' ratios of quadratic forms in normal variables using recursive algorithms.
#' That package was originally developed for evaluating average evolvability
#' measures (Watanabe, 2022), but is capable of evaluating moments in rather
#' general conditions beyond those.
#' The idea of this package is to provide a simple, convenient interface
#' specifically aiming at average evolvability measures, by passing
#' appropriately specified arguments to functions from the \code{qfratio}
#' package.
#' All average evolvability measures treated by Watanabe (2022) are implemented,
#' accommodating arbitrary mean and covariance for the selection gradients.
#'
#' The DESCRIPTION file:
#' \packageDESCRIPTION{avrevol}
#' \packageIndices{avrevol}
#'
#' @section Author/Maintainer:
#' Junya Watanabe <jw2098@cam.ac.uk>
#'
#' @references
#' Cheverud, J. M. (1996) Quantitative genetic analysis of cranial morphology
#'   in the cotton-top (*Saguinus oedipus*) and saddle-back (*S. fuscicollis*)
#'   tamarines. *Journal of Evolutionary Biology*, **9**, 5--42.
#'   doi:[10.1046/j.1420-9101.1996.9010005.x](https://doi.org/10.1046/j.1420-9101.1996.9010005.x).
#'
#' Cheverud, J. M. & Marroig, G. (2007) Comparing covariance matrices: random
#'   skewers method compared to the common principal components model.
#'   *Genetics and Molecular Biology*, **30**, 461--469.
#'   doi:[10.1590/S1415-47572007000300027](https://doi.org/10.1590/S1415-47572007000300027).
#'
#' Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and
#'   constraint in multivariate characters. *Journal of Evolutionary Biology*,
#'   **21**, 1201--1219.
#'   doi:[10.1111/j.1420-9101.2008.01573.x](https://doi.org/10.1111/j.1420-9101.2008.01573.x).
#'
#' Marroig, G., Shirai, L. T., Porto A., de Oliveira, F. B., & De Conto, V.
#'   (2009) The evolution of modularity in the mammalian skull II:
#'   evolutionary consequences. *Evolutionary Biology*, **36**, 136--148.
#'   doi:[10.1007/s11692-009-9051-1](https://doi.org/10.1007/s11692-009-9051-1).
#'
#' Watanabe, J. (2022). Exact expressions and numerical evaluation of average
#'   evolvability measures for characterizing and comparing **G** matrices.
#'   *bioRxiv* preprint, 2022.11.02.514929.
#'   doi:[10.1101/2022.11.02.514929](https://doi.org/10.1101/2022.11.02.514929).
#'
#' @seealso
#'   \code{\link{avr_evol}}: Average evolvability measures using
#'                           series expression
#'
#'   \code{\link{hh_evol}}: Approximate average evolvability measures using
#'                          the delta method of Hansen & Houle
#'
#'   \code{\link{mc_evol}}: Monte Carlo sampling for (average) evolvability
#'                          measures
#'
#' @examples
#' ## Simple covariance matrices
#' nv <- 4
#' G1 <- diag(nv:1)
#' G2 <- diag(sqrt(1:nv))
#' nit <- 1000
#'
#' ## Average conditional evolvability using series expression
#' ## Inspect plot to check for convergence
#' (res_cevo <- avr_cevo(G1))
#' plot(res_cevo)
#'
#' ## Hansen & Houle's (2008) delta method approximation of the same
#' hh_cevo(G1)
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
#' ## Average response difference using series expression,
#' ## Hansen-Houle delta method approximation, and
#' ## Monte Carlo estimate
#' (res_rdif <- avr_rdif(G1, G2))
#' plot(res_rdif)
#' hh_rdif(G1, G2)
#' mean(mc_rdif(nit, G1, G2))
#'
#' ## Average response correlation using series expression and
#' ## its Monte Carlo estimate, aka "random skewers" correlation
#' (res_rcor <- avr_rcor(G1, G2))
#' plot(res_rcor)
#' mean(mc_rcor(nit, G1, G2))
#'
#' ## Advanced: Average evolvability under
#' ## non-spherical distribution of selection gradient
#' ## (the same works for other evolvability measures as well)
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#' (res_evol_nsph <- avr_evol(G1, mu = mu, Sigma = Sigma))
#' plot(res_evol_nsph)
#' mean(mc_evol(nit, G1, mu = mu, Sigma = Sigma))
#'
#' @docType package
#' @name avrevol-package
# #' @aliases avrevol
#'
NULL
