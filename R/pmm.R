#' @exportClass pmm
setClass("pmm")

#' @title Propensity Scores
#' @description `pmm` builds a logistic regression with the target as the treatment variable
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
#' @returns \code{pmm} returns an object of \code{\link[base::class]{class} "pmm"}
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm model.
#'
#' An object of class \code{"pmm"} is a list containing the following:
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
#' p.score <- pmm(nhefs.nmv)
#' p.score

pmm <- function(data, f = NA, p.scores = NA, type = "strata", simple = pkg.env$simple,
                grp.width = 0.1, quant = T, ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # check valid model type
  if(!type %in% c("strata", "std", "std.boot")) {
    stop("Invalid model type! Must be on of the following values: 'strata', 'std', or 'std.boot'")
  }

  # if no given propensity scores
  if(anyNA(p.scores)) {
    # if no formula provided
    if(is.na(as.character(f))[1]) {
      # override simple
      if(simple != pkg.env$simple) {
        f <- build_formula(out = pkg.env$treatment, cov = pkg.env$covariates,
                           data = data, simple = simple)
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
  if(type == "strata") {
    p.grp <- list()
    lookup <- NA
    if(quant) {
      # group by propensity percentiles w/ width = grp.width
      quants <- c(quantile(p.scores, probs=seq(0,1,grp.width)))
      p.grp <- cut(p.scores, breaks=quants)

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
      p.grp <- cut(p.scores,
                   breaks=quants)
      lookup <- data.frame("n" = table(p.grp))
      colnames(lookup) <- c("breaks", "n")
      table$grp.name <- paste("p.grp", seq(1:(length(quants)-1)), sep = "")
      levels(p.grp) <- 1:nrow(lookup) # rename levels
      lookup <- lookup[c("grp.name", "breaks", "n")]
    }

    model.f <- as.formula(paste(pkg.env$outcome, "~", pkg.env$treatment, "* p.grp"))
    model <- glm(model.f, data = data)

    cont_mat <- contrast_matrix(model, nrow(lookup),
                                c(paste("Effect of", pkg.env$treatment, "for p.score in", lookup$breaks)))

    cont_mat[1:nrow(lookup), names(model$coefficients)[2]] <- 1
    for(i in 2:nrow(lookup)) {
      cont_mat[i, paste(names(model$coefficients)[2],":",lookup$grp.name[[i]], sep = "")] <- 1
    }

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

  output <- list("call" = call, "formula" = call$formula, "model" = model, "p.scores" = p.scores,
                 "ATE" = ATE, "ATE.summary" = ATE.summary, "type" = type)

  class(output) <- "pmm"
  return(output)
}

#' @export
print.pmm <- function(x) {
  if(x$type == "strata") {
    print(x$model)
  }
}

#' @export
summary.pmm <- function(x) {
  if(x$type == "strata") {
    summary(x$model)
  }
}

#' @export
predict.pmm <- function(x, newdata=NULL) {
  if(x$type == "strata") {
    if(is.null(newdata)) {
      return(predict(x$model$model))
    }
    else{
      return(predict(x$model$model, newdata=newdata))
    }
  }
}


