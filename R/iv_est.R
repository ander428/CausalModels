#' @exportClass iv_est
setClass("iv_est")

#' @title Propensity Scores
#' @description `propensity_scores` builds a logistic regression with the target as the treatment variable
#' and the covariates as the independent variables.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:binomial]{binomial}}
#' NOTE: if this is changed, the outcome of the model may not be the probabilities and the results will not be valid.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{propensity_scores} returns an object of \code{\link[base::class]{class} "propensity_scores"}
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm model.
#'
#' An object of class \code{"propensity_scores"} is a list containing the following:
#'
#' \tabular{ll}{
#'  \code{call} \tab the matched call. \cr
#'  \tab \cr
#'  \code{formula} \tab the formula used in the model. \cr
#'  \tab \cr
#'  \code{model} \tab the underlying glm model. \cr
#'  \tab \cr
#'  \code{p.scores} \tab the estimated propensity scores.\cr
#' }
#'
#' @export
#'
#' @examples
#' library(causaldata)
#' data(nhefs)
#' nhefs.nmv <- nhefs[which(!is.na(nhefs$wt82)),]
#' nhefs.nmv$qsmk <- as.factor(nhefs.nmv$qsmk)
#'
#' confounders <- c("sex", "race", "age", "education", "smokeintensity",
#'                      "smokeyrs", "exercise", "active", "wt71")
#'
#' init_params(wt82_71, qsmk,
#'             covariates = confounders,
#'             data = nhefs.nmv)
#'
#' p.score <- propensity_scores(nhefs.nmv)
#' p.score

iv_est <- function(IV, data) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])
  IV_levels <- levels(data[[IV]])

  IV <- as.character(params$IV)
  if(!is.factor(data[[IV]])) {
    stop("Instrumental variable must be a factor")
  }
  else if(length(IV_levels) > 2) {
    stop("Instrumental variable must be binary")
  }

  # manual calculation of standard IV Estimator
  numer_1 <- mean(data[data[[IV]] == IV_levels[[2]],][[pkg.env$outcome]], na.rm = TRUE)
  numer_0 <- mean(data[data[[IV]] == IV_levels[[1]],][[pkg.env$outcome]], na.rm = TRUE)
  denom_1 <- mean(as.numeric(data[data[[IV]] == IV_levels[[2]],][[pkg.env$treatment]]), na.rm = TRUE)
  denom_0 <- mean(as.numeric(data[data[[IV]] == IV_levels[[1]],][[pkg.env$treatment]]), na.rm = TRUE)

  return((numer_1 - numer_0) / (denom_1 - denom_0))
}
