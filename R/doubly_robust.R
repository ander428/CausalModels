#' @title Doubly Robust Model
#' @description \code{`doubly_robust`} trains both an outcome model and a propensity model to generate predictions
#' for the outcome and probability of treatment respectively. By default, the model uses
#' \code{\link[=standardization]{standardization}} and \code{\link[=propensity_scores]{propensity_scores}} to form a
#' doubly-robust model between standardization and IP weighting. Alternatively, any outcome and treatment
#' models can be provided instead, but must be compatible with the \code{\link[stats:predict]{predict}} generic function in R.
#' Since many propensity models may not predict probabilities without additional arguments into the
#' predict function, the predictions themselves can be given for both the outcome and propensity scores.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param out.mod (optional) a regression model that predicts the outcome. NOTE: the model given
#' must be compatible with the \code{\link[stats:predict]{predict}} generic function.
#' @param p.mod (optional) a propensity model that predicts the probability of treatment. NOTE: the model given
#' must be compatible with the \code{\link[stats:predict]{predict}} generic function.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param scores (optional) use calculated outcome estimates.
#' @param p.f (optional) an object of class "formula" that overrides the default formula for the denominator of the IP
#' weighting function.
#' @param p.simple a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.family the family to be used in the underlying propensity model.
#' By default, this is set to \code{\link[stats:gaussian]{binomial}}.
#' @param p.scores (optional) use calculated propensity scores.
#' @param n.boot an integer value that indicates number of bootstrap iterations to calculate standard error.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{doubly_robust} returns an object of \code{\link[base:class]{class}} "doubly_robust".
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} model.
#'
#' An object of class \code{"doubly_robust"} is a list containing the following:
#'
#'  \item{out.call}{the matched call of the outcome model.}
#'  \item{p.call}{the matched call of the propensity model.}
#'  \item{out.model}{the underlying outcome model.}
#'  \item{p.model}{the underlying propensity model.}
#'  \item{y_hat}{the estimated outcome values.}
#'  \item{p.scores}{the estimated propensity scores.}
#'  \item{ATE}{the estimated average treatment effect (risk difference).}
#'  \item{ATE.summary}{a data frame containing the ATE, SE, and 95\% CI of the ATE. }
#'  \item{data}{the data frame used to train the model.}
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
#' # use alternative outcome model
#' out.mod <- propensity_matching(data = nhefs.nmv)
#' db.model <- doubly_robust(out.mod = out.mod,
#'                           data = nhefs.nmv)
#' db.model
#'
#' # give calculated outcome predictions and give formula for propensity scores
#' db.model <- doubly_robust(scores = predict(out.mod),
#'                           p.f = qsmk ~ sex + race + age,
#'                           data = nhefs.nmv)
#' db.model


doubly_robust <- function(data, out.mod = NULL, p.mod = NULL, f = NA,
                          family = gaussian(), simple = pkg.env$simple, scores = NA,
                          p.f = NA, p.simple = pkg.env$simple, p.family = binomial(),
                          p.scores = NA, n.boot = 50, ...) {

  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  if(anyNA(scores)) {
    # generate standardization if outcome model not specified
    if(is.null(out.mod)) {
      out.mod <- standardization(f=f, data = data, simple = simple, family = family)
    }

    scores <- predict(out.mod)
  }


  if(anyNA(p.scores)) {
    # generate propensity model if treatment model not specified
    if(is.null(p.mod)) {
      p.mod <- propensity_scores(f=p.f, data = data, family=p.family, simple=p.simple)
    }

    p.scores <- predict(p.mod)
  }

  boot_result <- boot(data=data, statistic = function(data, indices) {
    data <- data[indices,]
    scores <- scores[indices]
    p.scores <- p.scores[indices]
    return(doubly_robust_est(scores, p.scores, data))
  }, R = n.boot)

  # calculate 95% CI
  beta <- boot_result$t0
  SE <- sd(boot_result$t)
  ATE <- data.frame(
    "ATE" = beta,
    "SE" = SE,
    conf_int(beta, SE),
    check.names=FALSE
  )

  output <- list("out.call" = out.mod$call, "p.call" = p.mod$call, "out.model" = out.mod,
                 "p.model" = p.mod,  "y_hat" = scores, "p.scores" = p.scores,
                 "ATE" = beta, "ATE.summary" = ATE, "data" = data)

  class(output) <- "doubly_robust"
  return(output)
}

#' @export
print.doubly_robust <- function(x, ...) {
  cat("Outcome Model\r\n\r\n")
  cat("Call:\r\n")
  print(x$out.call, ...)
  cat("\r\nPredictions:\r\n")
  print(summary(x$y_hat))
  cat("\r\nPropensity Model\r\n\r\n")
  cat("Call:\r\n")
  print(x$p.call, ...)
  cat("\r\nPredictions:\r\n")
  print(summary(x$p.scores))
  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat("Estimate - ", x$ATE, "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.doubly_robust <- function(object, ...) {
  s1 <- summary(object$out.model, ...)
  s1$ATE <- object$ATE
  s2 <- summary(object$p.model, ...)
  s2$ATE <- object$ATE
  s <- list("out.summary" = s1, "p.summary" = s2)
  class(s) <- "summary.doubly_robust"
  return(s)
}

#' @export
print.summary.doubly_robust <- function(x, ...) {
  s1 <- x$out.summary
  s2 <- x$p.summary
  cat("Outcome Model Summary\r\n")
  print(s1, ...)
  cat("Propensity Model Summary\r\n")
  print(s2, ...)
}

#' @export
predict.doubly_robust <- function(object, ...) {
  return(doubly_robust_est(object$y_hat, object$p.scores, object$data))
}

doubly_robust_est <- function(S, W, data) {
  n <- length(S)
  tr <- as.numeric(levels(data[[pkg.env$treatment]]))[data[[pkg.env$treatment]]]

  return(sum(S - ((tr * (data[[pkg.env$outcome]]-S)) / W)) / n)
}
