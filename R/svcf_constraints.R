#' Constrain an SVCF candidate to its biological parameter space
#'
#' Internal helpers for the overlapping-CNV estimator. The unconstrained
#' candidate is preserved by `calc_svcf()`; these helpers supply its constrained
#' value and an auditable boundary status.
#'
#' @param x Numeric vector of unconstrained SVCF candidates.
#'
#' @return `constrain_svcf()` returns a numeric vector in `[0, 1]`, with
#'   nonfinite values changed to `NA`. `svcf_constraint_status()` returns one of
#'   `in_range`, `boundary_low`, `boundary_high`, or `nonfinite`.
#' @keywords internal
constrain_svcf <- function(x) {
  ifelse(is.finite(x), pmin(pmax(x, 0), 1), NA_real_)
}

#' @rdname constrain_svcf
#' @keywords internal
svcf_constraint_status <- function(x) {
  dplyr::case_when(
    !is.finite(x) ~ "nonfinite",
    x < 0 ~ "boundary_low",
    x > 1 ~ "boundary_high",
    TRUE ~ "in_range"
  )
}
