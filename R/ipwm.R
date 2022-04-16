#' @exportClass ipwm
setClass("ipwm")

#' @exportClass summary.ipwm
setClass("summary.ipwm")

#' @title Parametric IP Weighting
#' @description `ipwm` uses the \code{\link[propensity_scores]{propensity_scores}} function to generate inverse probability
#' weights. The weights can either be standardized weights or non-standardized weights. The weights are used to train a
#' general linear model whose coefficient for treatment represents the average treatment effect on the additive scale.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family (optional) the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple (optional) a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.f (optional) an object of class "formula" that overrides the default formula for the denominator of the IP
#' weighting function.
#' @param p.simple (optional) a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.family (optional) the family to be used in the underlying propensity model.
#' By default, this is set to \code{\link[stats:gaussian]{binomial}}.
#' @param p.scores (optional) use calculated propensity scores for the weights. If using standardized weights,
#' the numerator will still be modeled.
#' @param SW (optional) a boolean indicator to indicate the use of standardized weights. By default, this is set to true.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{ipwm} returns an object of \code{\link[base::class]{class} "ipwm"}.
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm outcome model.
#'
#' An object of class \code{"ipwm"} is a list containing the following:
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
#'  \code{ATE} \tab the estimated average treatment effect. \cr
#'  \tab \cr
#'  \code{ATE.summary} \tab a data frame containing the ATE, SE, and 95\% CI of the ATE. \cr
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
#' # model using all defaults
#' model <- ipwm(data = nhefs.nmv)
#' summary(model)
#'
#' # Model using calculated propensity scores and manual outcome formula
#' p.scores <- propensity_scores(nhefs.nmv)
#' model <- ipwm(wt82_71 ~ qsmk, p.scores = p.scores, data = nhefs.nmv)
#' summary(model)

ipwm <- function(data, f = NA, family = gaussian(), simple = pkg.env$simple,
                 p.f = NA, p.simple = pkg.env$simple, p.family = binomial(),
                 p.scores = NA, SW = T, ...) {

  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if user gives an outcome formula
  if(is.na(as.character(f))[1]) {
    f <- as.formula(paste(pkg.env$outcome, "~", pkg.env$treatment))
  }

  # if user given a propensity formula
  if(!is.na(as.character(p.f))[1]) {
    # if no propensity scores
    if(anyNA(p.scores)) {
      p.scores <- propensity_scores(p.f, data = data, family = p.family)$p.scores
    }
    else {
      message("Ignoring given propensity formula since propensity scores have been given.")
      message("Using given propensity scores.")
    }
  }
  # if no given propensity formula
  else {
    # if no propensity scores
    if(anyNA(p.scores)) {
      p.scores <- propensity_scores(pkg.env$f_tr, data = data, family = p.family)$p.scores
    }
  }

  if(SW) {
    numer_scores <- propensity_scores(as.formula(paste(pkg.env$treatment,"~1")), family = binomial(), data = data)$p.scores

    data$sw <- numer_scores / p.scores
  }

  else {
    data$sw <- 1 / p.scores
  }

  model <- glm(f, data=data, weights=sw, family = family, ...)
  model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"

  # calculate causal stats
  beta <- coef(model)[[2]]
  SE <- coef(summary(model))[2,2]
  ATE <- data.frame(
    "Beta" = beta,
    "SE" = SE,
    "2.5 %" = beta-qnorm(0.975)*SE,
    "97.5 %" = beta+qnorm(0.975)*SE,
    check.names=FALSE
  )
  row.names(ATE) <- pkg.env$treatment

  output <- list("call" = model$call, "formula" = model$call$formula, "model" = model,
                 "weights" = data$sw, "ATE" = beta, "ATE.summary" = ATE)

  class(output) <- "ipwm"
  return(output)
}

#' @export
print.ipwm <- function(x) {
  print(x$model)
  cat("\r\n")
  cat(pkg.env$treatment, "ATE:", "\r\n")
  cat("Estimate - ", x$ATE, "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.ipwm <- function(x) {
  s <- summary(x$model)
  s$ATE <- x$ATE.summary
  class(s) <- "summary.ipwm"
  return(s)
}

#' @export
print.summary.ipwm <- function(s) {
  class(s) <- "summary.glm"
  print(s)
  cat("ATE:", "\r\n")
  print(s$ATE)
  cat("\r\n")
}

#' @export
predict.ipwm <- function(x, newdata=NULL) {
  if(is.null(newdata)) {
    return(predict(x$model))
  }
  else{
    return(predict(x$model, newdata=newdata))
  }
}

