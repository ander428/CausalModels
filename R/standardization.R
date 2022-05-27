#' @title Parametric Standardization
#' @description `standardization` uses a standard \code{\link[stats:glm]{glm}} linear model to perform parametric standardization
#' by adjusting bias through including all confounders as covariates. The model will calculate during training both the risk difference
#' and the risk ratio. Both can be accessed from the model as well as estimates of the counterfactuals of treatment.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param n.boot an integer value that indicates number of bootstrap iterations to calculate standard error.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{standardization} returns an object of \code{\link[base:class]{class} "standardization"}.
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} model.
#'
#' An object of class \code{"standardization"} is a list containing the following:
#'
#' \item{call}{the matched call.}
#' \item{formula}{the formula used in the model.}
#' \item{model}{the underlying glm model.}
#' \item{ATE}{a data frame containing estimates of the treatment effect
#'  of the observed, counterfactuals, and risk metrics.}
#' \item{ATE.summary}{a data frame containing the ATE, SE, and 95\% CI of the ATE. }
#'
#' @export
#'
#' @examples
#' library(causaldata)
#'
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
#' model <- standardization(data = nhefs.nmv)
#' print(model)
#' summary(model)
#' print(model$ATE.summary)
#' print(model$ATE.summary$Estimate[[2]] -
#'       model$ATE.summary$Estimate[[3]]) # manually calculate risk difference

standardization <- function(data, f = NA, family = gaussian(), simple = pkg.env$simple, n.boot=50, ...) {

  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if no formula provided
  if(is.na(as.character(f))[1]) {
    # override simple
    if(simple != pkg.env$simple) {
      f <- build_formula(out = pkg.env$outcome, tr = pkg.env$treatment,
                         cov = pkg.env$covariates,
                         data = data, simple = simple)
    }
    # use default
    else {
      f <- formula(pkg.env$f_out)
    }
  }

  # make three copies of the dataset
  cp <- data
  cp$label <- "observed"
  tr0 <- cp
  tr0$label <- "cf_untreated"
  tr0[pkg.env$treatment] <- 0
  tr0[pkg.env$outcome] <- NA
  tr1 <- cp
  tr1$label <- "cf_treated"
  tr1[pkg.env$treatment] <- 1
  tr1[pkg.env$outcome] <- NA

  combined_data <- rbind(cp, tr0, tr1) # combine copies

  model_func <- function(data, indices, f, family, ...) {
    if(!anyNA(indices)) {
      data <- data[indices,]
    }

    # build model using all three copies
    model <- glm(formula = f, data = data, family = family, ...)
    data$Y_hat <- predict(model, data)

    # calculate means in each group
    means <- c(mean(data$Y_hat[data$label=="observed"]),     # estimated outcome of the observed
               mean(data$Y_hat[data$label=="cf_treated"]),   # estimated counterfactual of the treated
               mean(data$Y_hat[data$label=="cf_untreated"]), # estimated counterfactual of the untreated
               mean(data$Y_hat[data$label=="cf_treated"]) -  # estimated risk differnece
                 mean(data$Y_hat[data$label=="cf_untreated"]),
               mean(data$Y_hat[data$label=="cf_treated"]) /  # estimated risk ratio
                 mean(data$Y_hat[data$label=="cf_untreated"]))
    return(list("model" = model, "means" = means))
  }

  # build model and bootstrapped estimates
  result <- model_func(data=combined_data, indices=NA, f=f, family=family, ...)
  boot_result <- boot(data=combined_data, R=n.boot, f=f, family=family,
                      statistic = function(data, indices, f, family, ...) {
                        model_func(data, indices, f, family, ...)$means[[4]]
                      }, ...)

  # calculate 95% CI
  beta <- boot_result$t0
  SE <- sd(boot_result$t)
  ATE.summary <- data.frame(
    "Beta" = beta,
    "SE" = SE,
    conf_int(beta, SE),
    check.names=FALSE
  )

  model <- result$model
  ATE <- data.frame(Estimate = result$means)
  rownames(ATE) <- c("Observed effect", "Counterfactual (treated)",
                             "Counterfactual (untreated)", "Risk difference",
                             "Risk ratio")

  model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"
  output <- list("call" = model$call, "formula" = model$call$formula,
                 "model" = model, "ATE" = ATE, "ATE.summary" = ATE.summary)
  class(output) <- "standardization"
  return(output)
}

#' @export
print.standardization <- function(x, ...) {
  print(x$model, ...)
  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat("Estimate - ", x$ATE$Estimate[[4]], "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.standardization <- function(object, ...) {
  s <- summary(object$model, ...)
  s$ATE <- object$ATE.summary
  class(s) <- "summary.standardization"
  return(s)
}

#' @export
print.summary.standardization <- function(x, ...) {
  class(x) <- "summary.glm"
  print(x, ...)
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  print(x$ATE, row.names=FALSE)
  cat("\r\n")
}

#' @export
predict.standardization <- function(object, ...) {
    return(predict(object$model, ...))
}
