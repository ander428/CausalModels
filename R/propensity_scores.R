#' @title Propensity Scores
#' @description `propensity_scores` builds a logistic regression with the target as the treatment variable
#' and the covariates as the independent variables.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:binomial]{binomial}}
#' NOTE: if this is changed, the outcome of the model may not be the probabilities and the results will not be valid.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{propensity_scores} returns an object of \code{\link[base:class]{class} "propensity_scores"}
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} model.
#'
#' An object of class \code{"propensity_scores"} is a list containing the following:
#'
#'  \item{call}{the matched call.}
#'  \item{formula}{the formula used in the model.}
#'  \item{model}{the underlying glm model.}
#'  \item{p.scores}{the estimated propensity scores.}
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
#'
#' init_params(wt82_71, qsmk,
#'   covariates = confounders,
#'   data = nhefs.nmv
#' )
#'
#' p.score <- propensity_scores(nhefs.nmv)
#' p.score
#'
propensity_scores <- function(data, f = NA, simple = pkg.env$simple, family = binomial(), ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if no formula provided
  if (is.na(as.character(f))[1]) {
    # override simple
    if (simple != pkg.env$simple) {
      f <- build_formula(
        out = pkg.env$treatment, cov = pkg.env$covariates,
        data = data, simple = simple
      )
    }
    # use default
    else {
      f <- formula(pkg.env$f_tr)
    }
  }

  # build model and generate scores
  model <- glm(f, data = data, family = family, ...)

  scores <- predict(model, type = "response")

  # probability of treatment for the untreated
  # scores[which(data[[pkg.env$treatment]] == 0)] = 1 - scores[which(data[[pkg.env$treatment]] == 0)]

  model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"
  output <- list("call" = model$call, "formula" = model$call$formula, "model" = model, "p.scores" = scores)

  class(output) <- "propensity_scores"
  return(output)
}

#' @export
print.propensity_scores <- function(x, ...) {
  print(x$model, ...)
}

#' @export
summary.propensity_scores <- function(object, ...) {
  summary(object$model, ...)
}

#' @export
predict.propensity_scores <- function(object, ...) {
  return(predict(object$model, type = "response", ...))
}
