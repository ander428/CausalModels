
#' @title One Parameter G-Estimation of Structural Nested Mean Models
#' @description `gestimation` uses the \code{\link[=propensity_scores]{propensity_scores}} function to generate inverse probability
#' weights. The weights can either be standardized weights or non-standardized weights. A grid search is done on \eqn{\alpha} for the
#' to construct the best \eqn{\Beta} coefficient in the structural nested mean model. Alternatively, a linear mean model can be used
#' for a closed form estimator.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param grid a list of possible \eqn{\Beta} values that will be used in the grid search.
#' @param ids (optional) see documentation for \code{\link[geepack:geeglm]{geeglm}}. By default rownames of the data will be used.
#' @param f (optional) an object of class "formula" that overrides the default parameter. NOTE: for g-estimation this should be
#' a propensity formula.
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' @param simple (optional) a boolean indicator to build default formula with interactions for the g-estimation model.
#' If true, interactions will be excluded. If false, interactions will be included. By default, simple is set to false.
#' NOTE: \eqn{\Beta}  will be appended to the end of the formula
#' @param p.f (optional) an object of class "formula" that overrides the default formula for the denominator of the IP
#' weighting function.
#' @param p.simple (optional) a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' NOTE: if this is changed, the coefficient for treatment may not accurately represent the average causal effect.
#' @param p.family the family to be used in the underlying propensity model.
#' By default, this is set to \code{\link[stats:binomial]{binomial}}.
#' @param p.scores (optional) use calculated propensity scores for the weights. If using standardized weights,
#' the numerator will still be modeled.
#' @param SW a boolean indicator to indicate the use of standardized weights. By default, this is set to true.
#' @param n.boot (optional) an integer value that indicates number of bootstrap iterations to calculate standard error.
#' If no value is given, the standard error from the underlying linear model will be used. NOTE: when type is 'one.grid'
#' bootstrapping is not performed. By default, this is set to 100.
#' @param type the type of g-estimation to perform. It must be one of \code{"one.grid"} or \code{"one.linear"} for a
#' one parameter grid and linear mean model estimation respectively.
#' @param ... additional arguments that may be passed to the underlying model.
#'
#' @returns \code{gestimation} returns an object of \code{\link[base:class]{class} "gestimation"}.
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glm} or \code{geeglm} model.
#'
#' An object of class \code{"gestimation"} is a list containing the following:
#'
#'  \item{call}{the matched call.}
#'  \item{formula}{the formula used in the model.}
#'  \item{model}{the underlying glm model. If the model performed a grid search, this will be renamed 'best.model'}
#'  \item{weights}{the estimated IP weights.}
#'  \item{type}{returns the value used for the 'type' parameter.}
#'  \item{ATE}{the estimated average treatment effect (risk difference).}
#'  \item{ATE.summary}{a data frame containing the ATE, SE, and 95\% CI of the ATE. }
#'
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
#' grid <- seq(from = 1,to = 5, by = 0.05)
#' gest.model <- CausalModels::gestimation(nhefs.nmv, grid = grid, type = "one.grid")
#' gest.model$ATE.summary

#' @export
gestimation <- function(data, grid, ids = list(), f = NA, family = binomial(), simple = pkg.env$simple,
                        p.f = NA, p.simple = pkg.env$simple, p.family = binomial(), p.scores = NA,
                        SW = TRUE, n.boot = 100, type = "one.grid",  ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if no formula provided
  if(is.na(as.character(f))[1]) {
    # override simple
    if(simple != pkg.env$simple) {
      f <- build_formula(out = pkg.env$treatment, cov = pkg.env$covariates,
                         data = data, simple = p.simple)
    }
    # use default
    else {
      f <- formula(pkg.env$f_tr)
    }
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

  # use standardized weights
  if(SW) {
    numer_scores <- propensity_scores(as.formula(paste(pkg.env$treatment,"~1")), family = binomial(), data = data)$p.scores

    data$weights <- numer_scores / p.scores
  }
  else {
    data$weights <- 1 / p.scores
  }

  # perform grid search
  if(type == "one.grid") {
    if(length(grid) == 0) {
      errorCondition("'grid' parameter must have be a list with at least one value.")
    }

    coefs <- cbind(rep(NA,length(grid)), rep(NA,length(grid)), rep(NA,length(grid)), rep(NA, length(grid)))
    colnames(coefs) <- list("Effect", "Estimate", "Std. Error", "Pr(>|W|)")
    data$weights <- rep(1, nrow(data))

    # validate rownames
    if(length(ids) != nrow(data)) {
      message("IDs list does not match number of rows in data. Using rownames by default.")
      data$ids <- rownames(data)
    }
    else {
      data$ids <- ids
    }

    f_terms <- as.character(f)
    f <- as.formula(paste(f_terms[[2]], f_terms[[1]], f_terms[[3]], "+", "beta"))

    i <- 0
    model <- NA
    models <- list()

    for (val in grid) {
      i <- i + 1
      data[[pkg.env$treatment]] <- as.numeric(as.character(data[[pkg.env$treatment]])) # geeglm requires numeric response
      data$beta <- data[[pkg.env$outcome]] - (val * data[[pkg.env$treatment]]-1)

      model <- geepack::geeglm(f, family=family, data=data, weights=weights, id=ids, corstr="independence")
      models[[i]] <- model
      estimate <- summary(model)$coefficients["beta", "Estimate"]

      coefs[i,1] <- val
      coefs[i,2] <- estimate
      coefs[i,3] <- summary(model)$coefficients["beta", "Std.err"]
      coefs[i,4] <- summary(model)$coefficients["beta", "Pr(>|W|)"]
    }

    result <- as.data.frame(coefs)
    rownames(result) <- paste("Effect =", result$Effect)

    conf_int_idx <- which(result$`Pr(>|W|)` > 0.05)

    ATE.summary <- data.frame(
      Beta = result$Effect[[which.min(abs(result$Estimate))]],
      SE = NA,
      `2.5 %` = result$Effect[[min(conf_int_idx)]],
      `97.5 %` = result$Effect[[max(conf_int_idx)]],
      check.names = FALSE
    )
    #colnames(ATE.summary) <- c("Beta", "SE", "2.5 %", "97.5 %")

    output <- list("call" = model$call, "formula" = f, "best.model" = models[[which.min(abs(result$Estimate))]],
                   "weights" = data$weights, "result" = result[-1], type = "one.grid",
                   "ATE" = ATE.summary$Beta, "ATE.summary" = ATE.summary)

    class(output) <- "gestimation"
    return(output)
  }

  # one parameter linear mean model
  else if(type == "one.linear") {
    data[[pkg.env$treatment]] <- as.numeric(as.character(data[[pkg.env$treatment]]))

    model_func <- function(data, indices, f, family, weights, ...) {
      if(!anyNA(indices)) {
        data <- data[indices,]
      }

      model <- glm(f, data = data, weights = weights, family = family)
      preds <- predict(model, data, type = "response")
      estimate <- sum(data$weights*data[[pkg.env$outcome]]*(data[[pkg.env$treatment]]-preds)) /
        sum(data$weights*data[[pkg.env$treatment]]*(data[[pkg.env$treatment]] - preds))
      return(list("model" = model, "ATE" = estimate))
    }

    # build model
    result <- model_func(data=data, indices=NA, f=f, family=family, weights = data$weights, ...)
    model <- result$model
    beta <- 0
    SE <- 0
    ATE <- result$ATE

    if(n.boot > 1) {
      # build bootstrapped estimates
      boot_result <- boot(data=data, R=n.boot, f=f, family=family, weights=data$weights,
                          statistic = function(data, indices, f, family, ...) {
                            model_func(data, indices, f, family, weights, ...)$ATE
                          }, ...)

      # calculate 95% CI
      beta <- boot_result$t0
      SE <- sd(boot_result$t)
    }

    # calculate causal stats
    ATE.summary <- data.frame(
      "Beta" = beta,
      "SE" = SE,
      conf_int(beta, SE),
      check.names=FALSE
    )

    output <- list("call" = model$call, "formula" = f, "model" = model, type = "one.linear",
                   "weights" = data$weights, "ATE" = ATE.summary$Beta, "ATE.summary" = ATE.summary)

    class(output) <- "gestimation"
    return(output)
  }

  else {
    errorCondition("Invalid model type. Must be one of ('one.grid', 'one.linear')")
  }
}

#' @export
print.gestimation <- function(x, ...) {
  if(x$type == "one.grid") {
    cat("Best Model:")
    cat("\r\n")
    print(x$best.model, ...)
  }
  else {
    print(x$model, ...)
  }

  cat("\r\n")
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  cat("Estimate - ", x$ATE, "\r\n")
  cat("SE       - ", x$ATE.summary$SE, "\r\n")
  cat("95% CI   - (", x$ATE.summary$`2.5 %`, ", ", x$ATE.summary$`97.5 %`, ")", "\r\n")
}

#' @export
summary.gestimation <- function(object, ...) {
  model <- ifelse(object$type == "one.grid", object$best.model, object$model)
  s <- summary(model, ...)
  s$ATE <- object$ATE.summary
  class(s) <- "summary.gestimation"
  return(s)
}

#' @export
print.summary.gestimation <- function(x, ...) {
  class(x) <- "summary.glm"
  print(x, ...)
  cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
  print(x$ATE, row.names = FALSE)
  cat("\r\n")
}

#' @export
predict.gestimation <- function(object, ...) {
  model <- ifelse(object$type == "one.grid", object$best.model, object$model)
  return(predict(model, ...))
}
