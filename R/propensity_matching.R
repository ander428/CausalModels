#' @title Propensity Matching
#' @description `propensity_matching` uses either stratification or standardization to model an outcome
#' conditional on the propensity scores. In stratification, the model will break the propensity scores
#' into groups and output a \code{\link[multcomp:glht]{glht}} model based off a contrast matrix which
#' estimates the change in average causal effect within groups of propensity scores. In standardization,
#' the model will output a \code{\link[=standardization]{standardization}} model that conditions on the
#' propensity strata rather than the covariates. The model can also predict the expected outcome.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param p.scores (optional) use calculated propensity scores for matching. Otherwise, propensity scores
#' will be automatically modeled.
#' @param p.simple a boolean indicator to build default formula with interactions for the propensity models.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param type a string representing the type of propensity model to be used. By default, the function will stratify. Standardization with
#' propensity scores may also be used. The value given for \code{type} must be in \code{c("strata", "stdm")}.
#' @param grp.width a decimal value to specify the range to stratify the propensity scores. If option \code{quant} is set to true,
#' this will represent the spread of percentiles. If false, it will represent the spread of raw values of propensity
#' scores. Must be a decimal between 0 and 1. By default, this is set to 0.1. This option is ignored for standardization.
#' @param quant a boolean indicator to specify the type of stratification. If true (default), the model will stratify by
#' percentiles. If false, the scores will be grouped by a range of their raw values. This option is ignored for standardization.
#' @param ... additional arguments that may be passed to the underlying \code{\link[=propensity_scores]{propensity_scores}} function.
#'
#' @returns \code{propensity_matching} returns an object of \code{\link[base:class]{class} "propensity_matching"}
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with the underlying \code{glht} or
#' \code{standardization} model.
#'
#' An object of class \code{"propensity_matching"} is a list containing the following:
#'
#'  \item{call}{the matched call.}
#'  \item{formula}{the formula used in the model.}
#'  \item{model}{either the underlying \code{glht} or \code{standardization} model.}
#'  \item{p.scores}{the estimated propensity scores}
#'  \item{ATE}{a data frame containing the ATE, SE, and 95\% CI of the ATE. }
#'  \item{ATE.summary}{either a data frame containing the \code{glht} or \code{standardization} summary. }
#'
#' @export
#'
#' @examples
#' library(causaldata)
#' library(multcomp)
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
#' pm.model <- propensity_matching(nhefs.nmv)
#' pm.model$ATE.summary
#' summary(pm.model)
#' head(data.frame(preds=predict(pm.model)))

propensity_matching <- function(data, f = NA, simple = pkg.env$simple, p.scores = NA, p.simple = pkg.env$simple,
                                type = "strata", grp.width = 0.1, quant = TRUE, ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # check valid model type
  if(!type %in% c("strata", "stdm")) {
    stop("Invalid model type! Must be on of the following values: 'strata', 'std', or 'std.boot'")
  }

  if(grp.width <= 0 || grp.width >= 1) {
    stop("Invalid parameter! You must set 0 < grp.width < 1")
  }

  # if no given propensity scores
  if(anyNA(p.scores)) {
    # if no formula provided
    if(is.na(as.character(f))[1]) {
      # override simple
      if(p.simple != pkg.env$simple) {
        f <- build_formula(out = pkg.env$treatment, cov = pkg.env$covariates,
                           data = data, simple = p.simple)
      }
      # use default
      else {
        f <- formula(pkg.env$f_tr)
      }
    }

    p.scores <- propensity_scores(f = f, data = data, ...)$p.scores
  }

  # initialize output objects
  model <- NA
  call <- NA
  ATE <- NA
  ATE.summary <- NA

  # using stratification
  if(type == "strata") {
    p.grp <- list()
    lookup <- NA
    if(quant) { # use percentiles
      # group by propensity percentiles w/ width = grp.width
      quants <- c(quantile(p.scores, probs=seq(0,1,grp.width)))
      p.grp <- cut(p.scores, breaks=quants, include.lowest = TRUE)

      # create lookup table for groups
      lookup <- data.frame("n" = table(p.grp), names(quants)[-1])
      colnames(lookup) <- c("breaks", "n", "percentile")
      lookup$grp.name <- paste("p.grp", seq(1:(length(quants)-1)), sep = "")

      levels(p.grp) <- 1:nrow(lookup) # rename levels
      lookup <- lookup[c("grp.name", "breaks", "n", "percentile")]
    }
    else { # use raw groupings
      # group by propensity breaks w/ width = grp.width
      quants <- seq(0,1,grp.width)
      p.grp <- cut(p.scores, breaks=quants, include.lowest = T)

      # create lookup table for groups
      lookup <- data.frame("n" = table(p.grp))
      colnames(lookup) <- c("breaks", "n")
      lookup$grp.name <- paste("p.grp", seq(1:(length(quants)-1)), sep = "")

      levels(p.grp) <- 1:nrow(lookup) # rename levels
      lookup <- lookup[c("grp.name", "breaks", "n")]
    }

    # build linear model to make estimates in the contrast matrix
    model.f <- as.formula(paste(pkg.env$outcome, "~", pkg.env$treatment, "* p.grp"))
    model <- glm(model.f, data = data)
    model$call$formula <- model.f

    # build contrast matrix of all propensity groups
    cont_mat <- contrast_matrix(model, nrow(lookup),
                                c(paste("Effect of", pkg.env$treatment, "for p.score in", lookup$breaks)))


    # set all treatment values to 1
    cont_mat[1:nrow(lookup), names(model$coefficients)[2]] <- 1

    # try to fill the diag of the matrix with 1s
    tryCatch(
      for(i in 2:nrow(lookup)) {
        cont_mat[i, paste(names(model$coefficients)[2],":",lookup$grp.name[[i]], sep = "")] <- 1
      },
      # if this fails, there are likely one or more groups with <= 1 samples
      error = function(e) {
        stop("Unable to stratify propensity scores. This is likely due to a lack of positivity in the groups.
           Try setting 'grp.width' to a larger value.")
      })

    # build model for contrast matrix
    model <- glht(model, cont_mat)
    call <- model$model$call

    # summarize model
    sum_model <- summary(model)$test

    results <- data.frame(
      "Estimate" = sum_model$coefficients,
      "Std. Error" = sum_model$sigma,
      conf_int(sum_model$coefficients, sum_model$sigma),
      check.names=FALSE
    )

    ATE.summary <- cbind(lookup, results)
    ATE <- results
  }

  # using standardization
  else if(type == "stdm") {
    model.f <- build_formula(out = pkg.env$outcome, tr=pkg.env$treatment,
                             cov=c("p.scores"), simple=simple, data = cbind(data, p.scores))
    model <- standardization(f = model.f, data = cbind(data, p.scores))

    call <- model$call
    call$formula <- model.f
    ATE <- model$ATE
    ATE.summary <- model$ATE.summary
  }

  output <- list("call" = call, "formula" = call$formula, "model" = model, "p.scores" = p.scores,
                 "ATE" = ATE, "ATE.summary" = ATE.summary, "type" = type)

  class(output) <- "propensity_matching"
  return(output)
}

#' @export
print.propensity_matching <- function(x, ...) {
  if(x$type == "strata") {
    cat("\r\n")
    cat("Average treatment effect of ", pkg.env$treatment, ":", "\r\n", sep = "")
    cat("\r\n")
    print(x$ATE)
  }

  else if(x$type == "stdm") {
    print(x$model, ...)
  }
}

#' @export
summary.propensity_matching <- function(object, ...) {
  if(object$type == "strata" || object$type == "stdm") {
    summary(object$model, ...)
  }
}

#' @export
predict.propensity_matching <- function(object, ...) {
  if(object$type == "strata") {
    return(predict(object$model$model, ...))
  }
  else if(object$type == "stdm") {
    return(predict(object$model, ...))
  }
}
