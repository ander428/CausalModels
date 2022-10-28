#' @title Parametric IP Weighting
#' @description `ipweighting` uses the \code{\link[=propensity_scores]{propensity_scores}} function to generate inverse probability
#' weights. The weights can either be standardized weights or non-standardized weights. The weights are used to train a
#' general linear model whose coefficient for treatment represents the average treatment effect on the additive scale.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param p.f (optional) an object of class "formula" that overrides the default formula for the denominator of the IP
#' weighting function.
#' @param p.simple a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.family the family to be used in the underlying propensity model.
#' By default, this is set to \code{\link[stats:binomial]{binomial}}.
#' @param p.scores (optional) use calculated propensity scores for the weights. If using standardized weights,
#' the numerator will still be modeled.
#' @param SW a boolean indicator to indicate the use of standardized weights. By default, this is set to true.
#' @param n.boot (optional) an integer value that indicates number of bootstrap iterations to calculate standard error.
#' If no value is given, the standard error from the underlying linear model will be used.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{ipweighting} returns an object of \code{\link[base:class]{class} "ipweighting"}.
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} model.
#'
#' An object of class \code{"ipweighting"} is a list containing the following:
#'
#'  \item{call}{the matched call.}
#'  \item{formula}{the formula used in the model.}
#'  \item{model}{the underlying glm model.}
#'  \item{weights}{the estimated IP weights.}
#'  \item{ATE}{the estimated average treatment effect (risk difference).}
#'  \item{ATE.summary}{a data frame containing the ATE, SE, and 95\% CI of the ATE. }
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
#' model <- ipweighting(data = nhefs.nmv)
#' summary(model)
#'
#' # Model using calculated propensity scores and manual outcome formula
#' p.scores <- propensity_scores(nhefs.nmv)$p.scores
#' model <- ipweighting(wt82_71 ~ qsmk, p.scores = p.scores, data = nhefs.nmv)
#' summary(model)

ipweighting <- function(data, f = NA, family = gaussian(), p.f = NA, p.simple = pkg.env$simple,
                        p.family = binomial(), p.scores = NA, SW = TRUE, n.boot = 0, ...) {

  check_init()
  data$weights <- rep(1, nrow(data))

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if user gives an outcome formula
  if(is.na(as.character(f))[1]) {
    f <- as.formula(paste(pkg.env$outcome, "~", pkg.env$treatment))
  }

  # if user does not give propensity scores
  if(anyNA(p.scores)) {
    if(!is.na(as.character(p.f))[1]) {
      p.scores <- propensity_scores(p.f, data = data, family = p.family)$p.scores
    }
    # if no given propensity formula
    else {
      if(p.simple != pkg.env$simple) {
        p.f <- build_formula(out = pkg.env$treatment, cov = pkg.env$covariates,
                             data = data, simple = p.simple)
      }
      # use default
      else {
        p.f <- formula(pkg.env$f_tr)
      }

      p.scores <- propensity_scores(p.f, data = data, family = p.family)$p.scores
    }

  }
  # if user does give propensity scores
  else {
    if(!is.na(as.character(p.f))[1]) {
      message("Ignoring given propensity formula since propensity scores have been given.")
    }
    message("Using given propensity scores.")
  }

  if(SW) {
    numer_scores <- propensity_scores(as.formula(paste(pkg.env$treatment,"~1")), family = binomial(), data = data)$p.scores

    data$weights <- numer_scores / p.scores
  }

  else {
    data$weights <- 1 / p.scores
  }

  model_func <- function(data, indices, f, family, weights, ...) {
    if(!anyNA(indices)) {
      data <- data[indices,]
    }

    model <- glm(f, weights=weights, data=data, family = family, ...)
    model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"

    return(list("model" = model, "ATE" = coef(model)[[2]]))
  }

  # build model
  result <- model_func(data=data, indices=NA, f=f, family=family, weights = data$weights, ...)
  model <- result$model
  beta <- 0
  SE <- 0
  ATE <- list()
  if(n.boot > 1) {
    # build bootstrapped estimates
    boot_result <- boot(data=data, R=n.boot, f=f, family=family, weights = data$weights,
                        statistic = function(data, indices, f, family, weights, ...) {
                          model_func(data, indices, f, family, weights, ...)$ATE
                        }, ...)

    # calculate 95% CI
    beta <- boot_result$t0
    SE <- sd(boot_result$t)
    ATE <- data.frame(
      "Beta" = beta,
      "SE" = SE,
      conf_int(beta, SE),
      check.names=FALSE
    )
  }

  else {
    # calculate causal stats
    beta <- coef(model)[[2]]
    SE <- coef(summary(model))[2,2]
    ATE <- data.frame(
      "Beta" = beta,
      "SE" = SE,
      conf_int(beta, SE),
      check.names=FALSE
    )
  }


  output <- list("call" = model$call, "formula" = model$call$formula, "model" = model,
                 "weights" = data$weights, "ATE" = beta, "ATE.summary" = ATE)

  class(output) <- "ipweighting"
  return(output)
}

#' @export
print.ipweighting <- function(x, ...) {
  print(x$model, ...)
  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat("Estimate - ", x$ATE, "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.ipweighting <- function(object, ...) {
  s <- summary(object$model, ...)
  s$ATE <- object$ATE.summary
  class(s) <- "summary.ipweighting"
  return(s)
}

#' @export
print.summary.ipweighting <- function(x, ...) {
  class(x) <- "summary.glm"
  print(x, ...)
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  print(x$ATE, row.names = FALSE)
  cat("\r\n")
}

#' @export
predict.ipweighting <- function(object, ...) {
  return(predict(object$model, ...))
}

