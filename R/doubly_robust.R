#' @title Parametric IP Weighting
#' @description `doubly_robust` uses the \code{\link[=propensity_scores]{propensity_scores}} function to generate inverse probability
#' weights. The weights can either be standardized weights or non-standardized weights. The weights are used to train a
#' general linear model whose coefficient for treatment represents the average treatment effect on the additive scale.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param p.f (optional) an object of class "formula" that overrides the default formula for the denominator of the IP
#' weighting function.
#' @param p.simple a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.family the family to be used in the underlying propensity model.
#' By default, this is set to \code{\link[stats:gaussian]{binomial}}.
#' @param p.scores (optional) use calculated propensity scores for the weights. If using standardized weights,
#' the numerator will still be modeled.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{doubly_robust} returns an object of \code{\link[base:class]{class} "doubly_robust"}.
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} model.
#'
#' An object of class \code{"doubly_robust"} is a list containing the following:
#'
#' \tabular{ll}{
#'  \code{call} \tab the matched call. \cr
#'  \tab \cr
#'  \code{formula} \tab the formula used in the model. \cr
#'  \tab \cr
#'  \code{model} \tab the underlying glm model. \cr
#'  \tab \cr
#'  \code{weights} \tab the estimated IP weights. \cr
#'  \tab \cr
#'  \code{ATE} \tab the estimated average treatment effect (risk difference). \cr
#'  \tab \cr
#'  \code{ATE.summary} \tab a data frame containing the ATE, SE, and 95\% CI of the ATE. \cr
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
#' # model using all defaults
#' model <- doubly_robust(data = nhefs.nmv)
#' summary(model)
#'
#' # Model using calculated propensity scores and manual outcome formula
#' p.scores <- propensity_scores(nhefs.nmv)$p.scores
#' model <- doubly_robust(wt82_71 ~ qsmk, p.scores = p.scores, data = nhefs.nmv)
#' summary(model)

doubly_robust <- function(data, f = NA, family = gaussian(), simple = pkg.env$simple,
                          p.f = NA, p.simple = pkg.env$simple, p.family = binomial(), p.scores = NA, ...) {

  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # generate both models
  prop.model <- propensity_scores(f=p.f, data = data, family=p.family, simple=p.simple)
  std.model <- standardization(f=f, data = data, simple = simple, family = family)

  est <- doubly_robust_est(std.model, prop.model, data)

  # output <- list("call" = model$call, "formula" = model$call$formula, "model" = model,
  #                "weights" = data$sw, "ATE" = beta, "ATE.summary" = ATE)
  #
  # class(output) <- "doubly_robust"
  return(est)
}

#' @export
print.doubly_robust <- function(x, ...) {
  print(x$model, ...)
  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat("Estimate - ", x$ATE, "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.doubly_robust <- function(object, ...) {
  s <- summary(object$model, ...)
  s$ATE <- object$ATE.summary
  class(s) <- "summary.doubly_robust"
  return(s)
}

#' @export
print.summary.doubly_robust <- function(x, ...) {
  class(x) <- "summary.glm"
  print(x, ...)
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  print(x$ATE)
  cat("\r\n")
}

#' @export
predict.doubly_robust <- function(object, ...) {
  return(predict(object$model, ...))
}

doubly_robust_est <- function(std.mod, prop.mod, data) {
  S <- predict(std.mod)
  W <- prop.mod$p.scores
  n <- length(S)
  tr <- as.numeric(levels(data[[pkg.env$treatment]]))[data[[pkg.env$treatment]]]

  return(sum(S - ((tr * (data[[pkg.env$outcome]]-S)) / W)) / n)
}
