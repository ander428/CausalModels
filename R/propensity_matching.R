#' @exportClass propensity_matching
setClass("propensity_matching")

#' @title Propensity Scores
#' @description `propensity_matching` builds a logistic regression with the target as the treatment variable
#' and the covariates as the independent variables.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:binomial]{binomial}}
#' NOTE: if this is changed, the outcome of the model may not be the probabilities and the results will not be valid.
#' @param ... additional arguments that may be passed to the underlying \code{\link[stats:glm]{glm}} model.
#'
#' @returns \code{propensity_matching} returns an object of \code{\link[base::class]{class} "propensity_matching"}
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm model.
#'
#' An object of class \code{"propensity_matching"} is a list containing the following:
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
#' p.score <- propensity_matching(nhefs.nmv)
#' p.score

propensity_matching <- function(data, f = NA, simple = pkg.env$simple, p.scores = NA, p.simple = pkg.env$simple,
                                type = "strata", grp.width = 0.1, quant = T, ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # check valid model type
  if(!type %in% c("strata", "stdm", "stdm.boot")) {
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

  model <- NA
  call <- NA
  ATE <- NA
  ATE.summary <- NA

  # using stratification
  if(type == "strata") {
    p.grp <- list()
    lookup <- NA
    if(quant) {
      # group by propensity percentiles w/ width = grp.width
      quants <- c(quantile(p.scores, probs=seq(0,1,grp.width)))
      p.grp <- cut(p.scores, breaks=quants, include.lowest = T)

      # create lookup table for groups
      lookup <- data.frame("n" = table(p.grp), names(quants)[-1])
      colnames(lookup) <- c("breaks", "n", "percentile")
      lookup$grp.name <- paste("p.grp", seq(1:(length(quants)-1)), sep = "")

      levels(p.grp) <- 1:nrow(lookup) # rename levels
      lookup <- lookup[c("grp.name", "breaks", "n", "percentile")]
    }
    else {
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

    model.f <- as.formula(paste(pkg.env$outcome, "~", pkg.env$treatment, "* p.grp"))
    model <- glm(model.f, data = data)
    model$call$formula <- model.f

    cont_mat <- contrast_matrix(model, nrow(lookup),
                                c(paste("Effect of", pkg.env$treatment, "for p.score in", lookup$breaks)))


    cont_mat[1:nrow(lookup), names(model$coefficients)[2]] <- 1

    tryCatch(
      for(i in 2:nrow(lookup)) {
        cont_mat[i, paste(names(model$coefficients)[2],":",lookup$grp.name[[i]], sep = "")] <- 1
      },
      error = function(e) {
        stop("Unable to stratify propensity scores. This is likely due to a lack of positivity in the groups.
           Try setting 'grp.width' to a larger value.")
      })
    model <- glht(model, cont_mat)
    call <- model$model$call

    sum_model <- summary(model)$test

    results <- list(
      "Estimate" = sum_model$coefficients,
      "Std. Error" = sum_model$sigma,
      "z value" = sum_model$tstat,
      "P(>|z|)" = sum_model$pvalues
    )

    ATE.summary <- cbind(lookup, results)
    ATE <- ATE.summary["Estimate"]
  }

  # using standardization
  else if(type == "stdm") {
    model.f <- build_formula(out = pkg.env$outcome, tr=pkg.env$treatment,
                             cov=c("p.scores"), simple=simple, data = cbind(data, p.scores))
    model <- stdm(f = model.f, data = cbind(data, p.scores))

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
print.propensity_matching <- function(x) {
  if(x$type == "strata" || x$type == "stdm") {
    print(x$model)
  }
}

#' @export
summary.propensity_matching <- function(x) {
  if(x$type == "strata" || x$type == "stdm") {
    summary(x$model)
  }
}

#' @export
predict.propensity_matching <- function(x, newdata=NULL) {
  if(is.null(newdata)) {
    if(x$type == "strata") {
      return(predict(x$model$model))
    }
    else if(x$type == "stdm") {
      return(predict(x$model))
    }
  }
  else {
    if(x$type == "strata") {
      return(predict(x$model$model, newdata=newdata))
    }
    else if(x$type == "stdm") {
      return(predict(x$model, newdata=newdata))
    }
  }
}
