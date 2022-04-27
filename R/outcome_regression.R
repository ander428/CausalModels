#' @title Outcome Regression
#' @description `outcome_regression` builds a linear model using all covariates. The treatment effects are stratified
#' within the levels of the covariates. The model will automatically provide all discrete covariates in a contrast matrix.
#' To view estimated change in treatment effect from continuous variables, a list called \code{contrasts}, needs to be given
#' with specific values to estimate. A vector of values can be given for any particualar continuous variable.
#'
#' @param data a data frame containing the variables in the model.
#' This should be the same data used in \code{\link[=init_params]{init_params}}.
#' @param f (optional) an object of class "formula" that overrides the default parameter
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#' @param family the family to be used in the general linear model.
#' By default, this is set to \code{\link[stats:gaussian]{gaussian}}.
#' NOTE: if this is changed, the assumptions about the model output may be incorrect and may not provide
#' accurate treatment effects.
#' @param contrasts a list of continuous covariates and values in the model to be included in the contrast matrix
#' (e.g. \code{list(age = c(18, 25, 40), weight = c(90, 159))}).
#' @param ... additional arguments that may be passed to the underlying \code{\link[multcomp:glht]{glht}} model.
#'
#' @returns \code{outcome_regression} returns an object of \code{\link[base:class]{class} "outcome_regression"}
#'
#' The functions \code{print}, \code{summary}, and \code{predict} can be used to interact with
#' the underlying \code{glht} model.
#'
#' An object of class \code{"outcome_regression"} is a list containing the following:
#'
#'  \item{call}{the matched call.}
#'  \item{formula}{the formula used in the model.}
#'  \item{model}{the underlying glht model.}
#'  \item{ATE}{estimated change average treatment effects within each strata}
#'  \item{ATE.summary}{a more detailed summary of the ATE estimations from glht. }
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
#' out.mod <- outcome_regression(nhefs.nmv, contrasts = list(age = c(21, 55),
#'                               smokeintensity = c(5, 20, 40)))
#' print(out.mod)
#' summary(out.mod)
#' head(data.frame(preds=predict(out.mod)))

outcome_regression <- function(data, f = NA, simple = pkg.env$simple,
                               family = gaussian(), contrasts = list(),  ...) {
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

  model <- glm(f, data=data, family = family, ...)
  model$call$formula <- f

  # grab and divide model matrix
  temp_map <- data.frame(model.matrix(model))
  num_vars <- sapply(temp_map, function(x) {length(unique(x)) > 2})
  disc <- temp_map[!num_vars]
  disc <- disc[,-c(1,ncol(disc))]# exclude intercept and treatment
  has_interaction <- unlist(lapply(names(num_vars), function(x) {grepl(pkg.env$treatment,
                                                                       x, fixed = T)}))
  continuous <- temp_map[has_interaction][-1]


  # if there are no interactions, use original variables
  if(ncol(continuous) == 0) {
    continuous <- temp_map[num_vars]
  }
  # if only sum interactions, use original for those without one
  else {
    no_interaction <- list()
    for(i in 1:length(num_vars)) {
      no_interaction[i] <- (!has_interaction[[i]] && num_vars[[i]] &&
                              !any(grepl(names(num_vars)[[i]], names(continuous))))
    }
    continuous <- cbind(continuous, temp_map[unlist(no_interaction)])
  }

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
        if(paste(pkg.env$treatment, levels(data[[pkg.env$treatment]])[2],
                 ":", name, sep="") %in% colnames(cont_mat)) {
                   cont_mat[j, paste(pkg.env$treatment, levels(data[[pkg.env$treatment]])[2], ":", name, sep="")] <- val
        }
        else if(paste(name, ":", pkg.env$treatment,
                      levels(data[[pkg.env$treatment]])[2], sep="") %in% colnames(cont_mat)) {
          cont_mat[j, paste(name, ":", pkg.env$treatment, levels(data[[pkg.env$treatment]])[2], sep="")] <- val
        }
        else {
          cont_mat[j, name] <- val
        }
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
print.outcome_regression <- function(x, ...) {
  print(x$model, ...)
}

#' @export
summary.outcome_regression <- function(object, ...) {
  summary(object$model, ...)
}

#' @export
predict.outcome_regression <- function(object, ...) {
  return(predict(object$model$model, ...))
}
