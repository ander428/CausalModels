
build_formula <- function(out, cov = list(), tr = NA, data, simple = FALSE) {
  # set treatment variable value
  tr <- ifelse(!is.na(tr), as.character(tr), "")

  # separate continuous and distrete covariates
  num_vars <- unlist(lapply(data[cov], is.numeric))
  cont <- cov[num_vars]
  disc <- cov[!num_vars]
  glm
  # collapse discrete covariates with +
  disc_val <- paste(unlist(disc), collapse = " + ")
  disc_val <- ifelse(tr != "" && disc_val != "", paste("+", disc_val), disc_val) # add + before if prev term exists

  # return basic formula with no interactions
  if (simple) {
    cont_val <- paste(unlist(cont), collapse = " + ") # collapse all terms with +
    cont_val <- ifelse(disc_val != "" || tr != "", # add + before if prev term exists
      paste("+", cont_val), cont_val
    )

    # return formula that doesn't include a treatment
    if (tr == "") {
      return(as.formula(paste(out, "~", disc_val, cont_val, collapse = " ")))
    }
    # return a formula with a treatment
    else {
      return(as.formula(paste(out, "~", tr, disc_val, cont_val, collapse = " ")))
    }
  }
  # if treatment exists, add interations between treatment and continuous
  # add cubic terms for each continuous variable
  else {
    # helper function to write out interactions between two vars
    print_interact <- function(x, y, use.I = TRUE) {
      if (!use.I) {
        return(paste("(", x, "*", y, ")", collapse = ""))
      }

      return(paste("I(", x, "*", y, ")", collapse = ""))
    }

    cont_val <- ""
    # only execute if there are continuous covariates
    if (!anyNA(cont) && !identical(cont, character(0))) {
      # logic for a single continuous covariate
      if (length(cont) == 1) {
        tr_interact <- ifelse(tr != "", paste("+", print_interact(tr, cont, FALSE)), "")
        cont_val <- paste(cont, tr_interact, "+", print_interact(cont, cont))
      }
      # logic for a set of continuous covariates
      else {
        i <- 1
        for (var in cont) {
          tr_interact <- ifelse(tr != "", paste("+", print_interact(tr, var, FALSE)), "")
          plus_op <- ifelse(i == 1, "", "+") # if first value, don't add plus
          cont_val <- paste(cont_val, plus_op, var, tr_interact, "+", print_interact(var, var))
          i <- i + 1
        }
      }
      cont_val <- ifelse(tr != "" || disc_val != "", paste("+", cont_val), cont_val) # add + to front if prev terms exist
    }
  }
  if (tr == "") {
    return(as.formula(paste(out, "~", disc_val, cont_val, collapse = " ")))
  } else {
    return(as.formula(paste(out, "~", tr, disc_val, cont_val, collapse = " ")))
  }
}

check_init <- function() {
  if (!pkg.env$init) {
    stop("Parameters not initialized. Please run init_params(...) before using any other functions")
  }
}


conf_int <- function(x, se) {
  lb <- x - qnorm(0.975) * se
  ub <- x + qnorm(0.975) * se
  return(list("2.5 %" = lb, "97.5 %" = ub))
}


contrast_matrix <- function(model, rows, names) {
  mat <- matrix(0, nrow = rows, ncol = length(coef(model)))
  colnames(mat) <- names(coef(model))
  rownames(mat) <- names
  return(mat)
}
