#' @exportClass stdm
setClass("stdm")

#' @exportClass summary.stdm
setClass("summary.stdm")

#' @title Parametric Standardization
#' @description `stdm` uses a standard \code{\link[stats:glm]{glm}} linear model to perform parametric standardization
#' by adjusting bias through including all confounders as covariates. The coefficient of treatment is the estimated
#' average causal effect.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{stdm} returns an object of \code{\link[base::class]{class} "stdm"}.
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm outcome model.
#'
#' An object of class \code{"stdm"} is a list containing the following:
#'
#' \tabular{ll}{
#'  \code{call} \tab the matched call. \cr
#'  \tab \cr
#'  \code{formula} \tab the formula used in the model. \cr
#'  \tab \cr
#'  \code{model} \tab the underlying glm model. \cr
#'  \tab \cr
#'  \code{ATE} \tab the estimated average treatment effect. \cr
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
#' model <- stdm(data = nhefs.nmv)
#' summary(model)

stdm <- function(data, f = NA, family = gaussian(), simple = pkg.env$simple,
                 p.f = NA, p.simple = pkg.env$simple, p.family = binomial(),
                 p.scores = NA, SW = T, ...) {

  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if user gives an outcome formula
  if(is.na(as.character(f))[1]) {
    f <- pkg.env$f_out
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

  combined_data <- rbind(cp, tr0, tr1)

  model <- glm(f, data=combined_data, family = family, ...)
  combined_data$Y_hat <- predict(model, combined_data)

  means <- c(mean(combined_data$Y_hat[combined_data$label=="observed"]),     # estimated outcome of the observed
                   mean(combined_data$Y_hat[combined_data$label=="cf_treated"]),   # estimated counterfactual of the treated
                   mean(combined_data$Y_hat[combined_data$label=="cf_untreated"]), # estimated counterfactual of the untreated
                   mean(combined_data$Y_hat[combined_data$label=="cf_treated"]) -  # estimated risk differnece
                     mean(combined_data$Y_hat[combined_data$label=="cf_untreated"]),
                   mean(combined_data$Y_hat[combined_data$label=="cf_treated"]) /  # estimated risk ratio
                     mean(combined_data$Y_hat[combined_data$label=="cf_untreated"]))
  ATE.summary <- data.frame(Estimate = means)
  rownames(ATE.summary) <- c("Observed effect", "Counterfactual (treated)",
                            "Counterfactual (untreated)", "Risk difference",
                            "Risk ratio")

  model$call$formula <- formula(f) # manually set model formula to prevent "formula = formula"
  output <- list("call" = model$call, "formula" = model$call$formula,
                 "model" = model, "ATE" = ATE.summary$Estimate[[4]], "ATE.summary" = ATE.summary)
  class(output) <- "stdm"
  return(output)
}

#' @export
print.stdm <- function(x) {
  print(x$model)
  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat(x$ATE, "\r\n")
}

#' @export
summary.stdm <- function(x) {
  s <- summary(x$model)
  s$ATE <- x$ATE.summary
  class(s) <- "summary.stdm"
  return(s)
}

#' @export
print.summary.stdm <- function(s) {
  class(s) <- "summary.glm"
  print(s)
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  print(s$ATE)
  cat("\r\n")
}

#' @export
predict.stdm <- function(x, newdata=NULL) {
  if(is.null(newdata)) {
    return(predict(x$model))
  }
  else{
    return(predict(x$model, newdata=newdata))
  }
}

estimate_ate <- function(data, indices) {
  # make three copies of the dataset
  cp <- data[indices,]
  cp$label <- "observed"
  tr0 <- cp
  tr0$label <- "cf_untreated"
  tr0[pkg.env$treatment] <- 0
  tr0[pkg.env$outcome] <- NA
  tr1 <- cp
  tr1$label <- "cf_treated"
  tr1[pkg.env$treatment] <- 1
  tr1[pkg.env$outcome] <- NA

  combined_data <- rbind(cp, tr0, tr1)

  model <- glm(f, data=combined_data, family = family, ...)
  combined_data$Y_hat <- predict(model, combined_data)


  return(c(mean(combined_data$Y_hat[combined_data$label=="observed"]),     # estimated outcome of the observed
           mean(combined_data$Y_hat[combined_data$label=="cf_treated"]),   # estimated counterfactual of the treated
           mean(combined_data$Y_hat[combined_data$label=="cf_untreated"]), # estimated counterfactual of the untreated
           mean(combined_data$Y_hat[combined_data$label=="cf_treated"]) -  # estimated risk differnece
             mean(combined_data$Y_hat[combined_data$label=="cf_untreated"]),
           mean(combined_data$Y_hat[combined_data$label=="cf_treated"]) /  # estimated risk ratio
             mean(combined_data$Y_hat[combined_data$label=="cf_untreated"])))
}


