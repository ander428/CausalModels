#' @exportClass outcome_regression
setClass("outcome_regression")

#' @title Propensity Scores
#' @description `outcome_regression` builds a logistic regression with the target as the treatment variable
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
#' @returns \code{outcome_regression} returns an object of \code{\link[base::class]{class} "outcome_regression"}
#'
#' The function \code{summary} can be used to obtain and print a summary of the underlying glm model.
#'
#' An object of class \code{"outcome_regression"} is a list containing the following:
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
#' p.score <- outcome_regression(nhefs.nmv)
#' p.score

outcome_regression <- function(data, f = NA, contrasts = list(), simple = pkg.env$simple, ...) {
  check_init()

  # grab function parameters
  params <- as.list(match.call()[-1])

  # if no formula provided
  if(is.na(as.character(f))[1]) {
    # override simple
    if(simple != pkg.env$simple) {
      f <- build_formula(out = pkg.env$outcome, tr=pkg.env$treatment,
                         cov = pkg.env$covariates, data = data, simple = simple)
    }
    # use default
    else {
      f <- formula(pkg.env$f_out)
    }
  }

  model <- NA
  call <- NA
  ATE <- NA
  ATE.summary <- NA

  model <- glm(f, data=data)
  model$call$formula <- f

  # grab and divide model matrix
  temp_map <- data.frame(model.matrix(model))
  num_vars <- sapply(temp_map, function(x) {length(unique(x)) > 2})
  disc <- temp_map[!num_vars][-1:-2]  # exclude intercept and treatment
  continuous <- temp_map[unlist(lapply(names(num_vars), function(x) {grepl(pkg.env$treatment, x, fixed = T)}))][-1]
  continuous_names <- names(continuous)


  cont_input <- list()
  for(i in 1:length(names(disc))) {
    cont_input[[i]] <- paste("Effect of", "qsmk", "at", names(disc[i]))
  }

  if(length(contrasts) > 0) {
    j <- length(cont_input) + 1
    for(i in 1:length(contrasts)) {
      cont_input[[j]] <- paste("Effect of", "qsmk", "at", names(contrasts[i]), "of", contrasts[[i]]) # issue is here
      j <- j + 1
    }
  }

  cont_input <- unlist(cont_input)

  cont_mat <- contrast_matrix(model, length(cont_input), cont_input)

  sum_cont <- 0
  for(cont in contrasts) {
    sum_cont <- sum_cont + length(cont)
  }

  j <- 1
  cont_mat[1:length(cont_input), "qsmk1"] <- 1

  for(i in 1:length(names(disc))) {
    cont_mat[j, names(disc[i])] <- 1
    j <- j + 1
  }

  if(length(contrasts) > 0) {
    for(i in 1:length(contrasts)) {
      name <- names(contrasts[i])
      for(val in contrasts[[i]]) {
        cont_mat[j, paste(pkg.env$treatment, levels(data[[pkg.env$treatment]])[2], ":", name, sep="")] <- val
        j <- j + 1
      }
    }
  }


  model <- glht(model, cont_mat)
  call <- model$model$call

  ATE.summary <- summary(model)
  ATE <- ATE.summary$test$Estimate

  output <- list("call" = call, "formula" = call$formula, "model" = model,
                 "ATE" = ATE, "ATE.summary" = ATE.summary)

  class(output) <- "outcome_regression"
  return(output)
}

#' @export
print.outcome_regression <- function(x) {
  print(x$model)
}

#' @export
summary.outcome_regression <- function(x) {
  summary(x$model)
}

#' @export
predict.outcome_regression <- function(x, newdata=NULL) {
  if(is.null(newdata)) {
    return(predict(x$model$model))
  }
  else {
    return(predict(x$model$model, newdata=newdata))
  }
}
