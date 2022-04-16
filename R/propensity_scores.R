#' @exportClass propensity_scores
setClass("propensity_scores")

#' @title Propensity Scores
#' @description `propensity_scores` builds a logistic regression with the target as the treatment variable
#' and the covariates as the independent variables.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple (optional) a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param family (optional) the family to be used in the general linear model.
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
#' data(nhefs)
#' nhefs$cens <- ifelse(is.na(nhefs$wt82), 1, 0)
#'
#' nhefs.nmv <- nhefs[which(!is.na(nhefs$wt82)),] # provisionally ignore subjects with missing values for weight in 1982
#' nhefs.nmv$qsmk <- as.factor(nhefs.nmv$qsmk)
#'
#' confounders <- c("sex", "race", "age", "education", "smokeintensity", "smokeyrs", "exercise", "active", "wt71")
#' init_params(wt82_71, qsmk,
#'             covariates = confounders,
#'             data = nhefs.nmv, simple = T)
#'
#' p.score <- propensity_scores(nhefs.nmv)
#' p.score

propensity_scores <- function(data, f = NA,  simple = pkg.env$simple, family = binomial(), ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if formula provided, then override default
  if(is.na(as.character(f))[1]) {
    f <- formula(pkg.env$f_tr)
  }

  # build model and generate scores
  model <- glm(f, data = data, family = family, ...)

  scores <- predict(model, type="response")

  # probability of treatment for the untreated
  scores[which(data[[pkg.env$treatment]] == 0)] = 1 - scores[which(data[[pkg.env$treatment]] == 0)]

  model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"
  output <- list("call" = model$call, "formula" = model$call$formula, "model" = model, "p.scores" = scores)

  class(output) <- "propensity_scores"
  return(output)
}

#' @export
print.propensity_scores <- function(x) {
  print(x$model)
}

#' @export
summary.propensity_scores <- function(x) {
  summary(x$model)
}

#' @export
predict.propensity_scores <- function(x, newdata=NULL) {
  if(is.null(newdata)) {
    return(predict(x$model, type = 'response'))
  }
  else{
    return(predict(x$model, newdata=newdata, type = 'response'))
  }
}