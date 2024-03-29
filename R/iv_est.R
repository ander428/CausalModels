#' @title Standard Instrumental Variable Estimator
#' @description `iv_est` calculates the standard IV estimand using the conditional means on a given instrumental variable.
#'
#' @param IV the instrumental variable to be used in the conditional means. Must be a factor with no more than 2 levels.
#' It is assumed the second level is the positive level, i.e., the binary equivalent of the second factor level should be 1
#' and the first should be 0.
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param n.boot an integer value that indicates number of bootstrap iterations to calculate standard error.
#'
#' @returns \code{iv_est} returns a data frame containing the standard IV estimate, standard error, and Wald 95% CI.
#'
#' @export
#'
#' @examples
#' library(causaldata)
#' data(nhefs)
#' nhefs.nmv <- nhefs[which(!is.na(nhefs$wt82)), ]
#' nhefs.nmv$qsmk <- as.factor(nhefs.nmv$qsmk)
#'
#' confounders <- c(
#'   "sex", "race", "age", "education", "smokeintensity",
#'   "smokeyrs", "exercise", "active", "wt71"
#' )
#' nhefs.iv <- nhefs[which(!is.na(nhefs$wt82) & !is.na(nhefs$price82)), ]
#' nhefs.iv$highprice <- as.factor(ifelse(nhefs.iv$price82 >= 1.5, 1, 0))
#' nhefs.iv$qsmk <- as.factor(nhefs.iv$qsmk)

#' init_params(wt82_71, qsmk,
#'             covariates = confounders,
#'             data = nhefs.iv)
#'
#' iv_est("highprice", nhefs.iv)

iv_est <- function(IV, data, n.boot = 50) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])
  IV_levels <- levels(data[[IV]])

  IV <- as.character(params$IV)
  if (!is.factor(data[[IV]])) {
    stop("Instrumental variable must be a factor")
  } else if (length(IV_levels) > 2) {
    stop("Instrumental variable must be binary")
  }

  est_func <- function(data, indices, ...) {
    data <- data[indices, ]
    # manual calculation of standard IV Estimator
    numer_1 <- mean(data[data[[IV]] == IV_levels[[2]], ][[pkg.env$outcome]], na.rm = TRUE)
    numer_0 <- mean(data[data[[IV]] == IV_levels[[1]], ][[pkg.env$outcome]], na.rm = TRUE)
    denom_1 <- mean(as.numeric(data[data[[IV]] == IV_levels[[2]], ][[pkg.env$treatment]]), na.rm = TRUE)
    denom_0 <- mean(as.numeric(data[data[[IV]] == IV_levels[[1]], ][[pkg.env$treatment]]), na.rm = TRUE)

    return((numer_1 - numer_0) / (denom_1 - denom_0))
  }

  boot_result <- boot(data = data, statistic = est_func, R = n.boot)

  # calculate 95% CI
  beta <- boot_result$t0
  SE <- sd(boot_result$t)
  ATE <- data.frame(
    "ATE" = beta,
    "SE" = SE,
    conf_int(beta, SE),
    check.names = FALSE
  )

  return(ATE)
}
