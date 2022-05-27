#' @title Initialize CausalModels Package
#' @description This function is required to be run first before any other function can run.
#' This will set within the package the global outcome, treatment, and covariate functions for each model to use.
#'
#' @param outcome the outcome variable of interest (must be continuous).
#' @param treatment the treatment with the causal effect of interest on the outcome.
#' @param covariates  a list/vector of covariate names to be use for confounding adjustment.
#' @param data a data frame containing the variables in the model.
#' @param simple a boolean indicator to build default formula with interactions.
#' If true, interactions will be excluded. If false, interactions will be included. By
#' default, simple is set to false.
#'
#' @export
#'
init_params <- function(outcome, treatment, covariates, data, simple = FALSE) {
  params <- as.list(match.call()[-1])
  cov <- params$cov

  if(!is.numeric(data[[as.character(params$outcome)]])) {
    assign('init', FALSE, pkg.env)
    stop("Outcome variable must be numeric!")
  }

  if(!is.factor(data[[as.character(params$treatment)]])) {
    assign('init', FALSE, pkg.env)
    stop("Treatment must be of type factor!")
  }

  else if(length(levels(data[[as.character(params$treatment)]])) > 2 ) {
    assign('init', FALSE, pkg.env)
    stop(paste("Treatment must be binary!", as.character(params$treatment), "has more than two levels!"))
  }

  tryCatch(data[as.character(params$covariates)[-1]], error = function(e) {assign('init', TRUE, pkg.env); stop(e)})

  assign('outcome', as.character(params$outcome), pkg.env)
  assign('treatment', as.character(params$treatment), pkg.env)

  # if covariates are given as strings
  if(length(as.character(params$covariates)[-1]) == 0) {
    cov <- covariates
    assign('covariates', cov, pkg.env)
  }
  # if covariates are not given as strings
  else {
    assign('covariates', as.character(params$covariates)[-1], pkg.env)
  }

  assign('simple', simple, pkg.env)
  assign('p_simple', simple, pkg.env)

  f_out <- build_formula(out = as.character(params$outcome),
                           tr = as.character(params$treatment),
                           cov = pkg.env$covariates,
                           data = data, simple = simple)
  f_tr <- build_formula(out = as.character(params$treatment),
                         cov = pkg.env$covariates,
                         data = data, simple = simple)
  assign('f_out', f_out, pkg.env)
  assign('f_tr', f_tr, pkg.env)

  assign('init', TRUE, pkg.env)

  cat("Successfully initialized!\r\n\r\n")

  cat("Summary:\r\n\r\n")
  cat(paste("Outcome -", pkg.env$outcome, "\r\n"))
  cat(paste("Treatment -", pkg.env$treatment, "\r\n"))
  cat(paste("Covariates -", "[", paste(pkg.env$covariates,collapse = ', '), "]", "\r\n\r\n"))
  cat(paste("Size -", nrow(data), "x", ncol(data), "\r\n\r\n"))
  cat(paste("Default formula for outcome models:", "\r\n"))
  cat(paste(deparse(f_out, width.cutoff = 500), "\r\n\r\n", collapse=""))
  cat(paste("Default formula for propensity models:", "\r\n"))
  cat(paste(deparse(f_tr, width.cutoff = 500), "\r\n", collapse=""))
}
