# Main landpred function
# Take a formula parse then return landpred model object we can call predict, etc... on
#' Create a Landpred Object
#'
#' Parses the formula and data to create a landpred object used for landmark prediction.
#' Call `?landpred.pacakge` for more information on the legacy API.
#'
#' @param formula A formula object with a Surv object on the LHS and covariates on the RHS.
#' @param data The data frame.
#' @param discrete Logical, whether to use the discrete method (legacy).
#'
#' @importFrom survival survfit Surv coxph
#' @importFrom stats bw.nrd dnorm optimize glm glm.fit as.formula coef model.frame model.response nobs predict printCoefmat rexp sd setNames terms vcov
#' @import splines
#' @import sm
#' @import quantreg
#' @example inst/examples/new_api_workflow.R
#'
#' @return A landpred_object.
#' @export
landpred <- function(formula, data, discrete=FALSE) {

  tf <- terms(formula, specials=c("Surv"))
  surv_terms <- attr(tf, "specials")$Surv

  # Get survival terms of right by not including the one for response
  rhs_survival_terms <- surv_terms[surv_terms > attr(tf, "response")]

  # Throw error if more than one short covariate, and get the column of it
  short_cov <- NULL
  if(length(rhs_survival_terms) > 1) {
    stop("Only a singular short-term covariate can be included.")
  } else if(length(rhs_survival_terms) == 1) {
    short_cov = attr(tf, "variables")[[rhs_survival_terms[[1]] + 1]]
    short_cov <- deparse(short_cov)
  }

  mf <- model.frame(formula, data)

  x_l <- model.response(mf)
  response_expr <- attr(tf, "variables")[[attr(tf, "response") + 1]]
  x_l_name <- deparse(response_expr[[2]])

  # delta for long covariate
  d_l_name <- deparse(response_expr[[3]])

  x_s <- NULL
  x_s_name <- NULL
  d_s_name <- NULL
  if(length(rhs_survival_terms) != 0) {
    x_s <- if (!is.null(short_cov)) mf[[short_cov]] else NULL
    short_expr <- attr(tf, "variables")[[rhs_survival_terms[[1]] + 1]]
    x_s_name <- deparse(short_expr[[2]])
    d_s_name <- deparse(short_expr[[3]])
  }

  if(!inherits(x_l, "Surv") || (!is.null(x_s) && !inherits(x_s, "Surv"))) {
    stop("Response variable and short-term covariate must a survival object.")
  }

  covariates <- attr(tf, "term.labels")

  if(length(covariates) == 1 && !is.null(x_s)) {
    stop("Normal covariate must be provided with short term covariate")
  }

  if(!discrete && is.null(x_s)) {
    stop("Short term covariate must be provded with multivariate continuous models")
  }

  covariates <- covariates[!(covariates %in% c(short_cov))]


  Z <- mf[, covariates, drop=FALSE]

  if(discrete == TRUE && ncol(Z) > 1) {
    stop("Only a singular covariate is allowed if discrete=TRUE")
  }

  names <- list(
    x_l_name = x_l_name,
    d_l_name = d_l_name,
    x_s_name = x_s_name,
    d_s_name = d_s_name,
    covariates=covariates
  )

  new_landpred_object(
    x_l,
    x_s,
    Z,
    names=names,
    formula = formula,
    discrete = discrete
  )
}

# Create new landpred model
new_landpred_object <- function(x_l, x_s, Z, formula, names, discrete) {
  structure(
    list(
      X_L = x_l, X_S = x_s, Z = Z, formula = formula, discrete = discrete,
      names=names
    ),
    class="landpred_object"
  )
}

#' Print Method for Landpred Object
#'
#' Prints a summary of the landpred object.
#'
#' @param x A landpred_object.
#' @param ... Additional arguments.
#'
#' @export
print.landpred_object <- function(x, ...) {
  cat("\nCall:\n")
  cat("landpred(formula = ", deparse(x$formula), ")\n", sep="")
}

#' Summary Method for Landpred Object
#'
#' Prints a detailed summary of the landpred object.
#'
#' @param object A landpred_object.
#' @param ... Additional arguments.
#'
#' @export
summary.landpred_object <- function(object, ...) {
  cat("\nLandpred Object Summary\n")
  cat("Call get_model() to get time-specific model for t0 + tau\n\n")
  cat("Call:\n")
  cat("landpred(formula = ", deparse(object$formula), ")\n\n", sep="")
  cat(sprintf("Discrete: %-8s Short Covariate: %-8s N: %s\n", object$discrete, !is.null(object$X_S), length(object$X_L)))
}

#' Get Landpred Model (General)
#'
#' Creates a landpred model object for a specific landmark time and prediction window.
#' Dispatches to continuous or discrete model creation based on the landpred object type.
#'
#' @param landpred_obj A landpred object.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#' @param bw The bandwidth (for continuous models).
#' @param transform Transformation function (for continuous models).
#'
#' @return A landpred_model object (continuous or discrete).
#' @export
get_model <- function(landpred_obj, t0, tau, bw=NULL, transform=identity) {
  if(landpred_obj$discrete == FALSE) {
    glm_noinfo <- fit_glm_normal(landpred_obj, t0, tau)
    model <- new_landpred_model_continuous(
      landpred_obj, glm_noinfo, t0, tau, bw,
      transform
    )
  } else {
    model <- new_landpred_model_discrete(landpred_obj, t0, tau)
  }
  model
}

#' Optimize Bandwidth for Continuous Landpred Models
#'
#' Selects the optimal bandwidth by minimizing the Mean Squared Error (MSE) using cross-validation.
#'
#' @param landpred_obj A landpred object.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#' @param lower Lower bound for bandwidth search.
#' @param upper Upper bound for bandwidth search.
#' @param transform Transformation function for the short-term covariate (e.g., log). Default is identity.
#' @param reps Number of cross-validation repetitions. Default is 50.
#' @param train_prop Proportion of data used for training in each fold. Default is 0.66.
#'
#' @return The optimal bandwidth.
#' @export
optimize_bandwidth <- function(landpred_obj, t0, tau, lower = 0.05, upper = 5, transform = identity, reps = 50, train_prop = 0.66) {
  if (landpred_obj$discrete) {
    stop("Bandwidth optimization is only for continuous models.")
  }
  
  opt <- optimize(
    f = mse_cv,
    interval = c(lower, upper),
    landpred_obj = landpred_obj,
    t0 = t0,
    tau = tau,
    transform = transform,
    reps = reps,
    train_prop = train_prop
  )
  
  return(opt$minimum)
}

new_landpred_model_discrete <- function(landpred_obj, t0, tau) {
  structure(
    list(
      landpred_obj = landpred_obj,
      t0 = t0,
      tau = tau
    ),
    class = "landpred_model_discrete"
  )
}

#' Predict Method for Discrete Landpred Model
#'
#' Predicts probabilities using the discrete landpred model.
#'
#' @param object A landpred_model_discrete object.
#' @param newdata Optional new data.
#' @param ... Additional arguments.
#'
#' @return Predicted probabilities.
#' @export
predict.landpred_model_discrete <- function(object, newdata = NULL, ...) {
  handle_discrete_pred(object$landpred_obj, object$t0, object$tau, newdata)
}

#' Print Method for Discrete Landpred Model
#'
#' Prints the discrete landpred model results.
#'
#' @param x A landpred_model_discrete object.
#' @param ... Additional arguments.
#'
#' @export
print.landpred_model_discrete <- function(x, ...) {
  old_landpred_result <- get_old_landpred_results_discrete(x$landpred_obj, x$t0, x$tau)

  model_name <- ""
  model_prob_formula <- ""

  if(old_landpred_result$mode == "no-covariate") {
    model_name <- "No Discrete covariate + no short covariate"
    model_prob_formula <- "P(TL < t0 + tau)"
  } else if (old_landpred_result$mode == "single-covariate") {
    model_name <- "Discrete covariate + no short covariate"
    model_prob_formula <- "P(TL < t0 + tau|Z=z)"
  } else {
    model_name <- "Discrete covariate + short covariate"
    model_prob_formula <- "P(TL < t0 + tau|Z=z,T_s=t_s)"
  }

  cat(sprintf("\nDiscrete Landpred Model (%s):\n", model_name))
  print(x$landpred_obj)
  cat("\n")

  Probs <- old_landpred_result$Prob
  if(is.matrix(Probs)) {
    cat("Probs:\n")
    apply(Probs, 1, function(row) {
      cat(sprintf("P(TL < t0+tau|Z=%d): %.3f\n", row[1], row[2]))
    })
    cat("\n")
  } else if (!is.null(Probs)) {
    cat("Probs:\n")
    cat(sprintf("%s: %.3f", model_prob_formula, Probs))
    cat("\n\n")
  }

  cat(sprintf("t0: %-10.3f tau: %-10.3f", x$t0, x$tau))
}
handle_discrete_pred <- function(landpred_obj, t0, tau, newdata) {
  old_landpred_result <- get_old_landpred_results_discrete(landpred_obj, t0, tau, newdata)
  probs <- old_landpred_result$newdata[, "Probability", drop = TRUE]
  probs
}

# Wrapper around old landpred functions
# Given a landpred object we call the aproppriate old function
get_old_landpred_results_discrete <- function(landpred_obj, t0, tau, newdata = NULL) {
  x_l <- landpred_obj$X_L
  x_s <- landpred_obj$X_S
  Z   <- landpred_obj$Z
  names_list <- landpred_obj$names

  formatted_data <- data.frame(
    XL = as.numeric(x_l[, "time"]),
    DL = as.numeric(x_l[, "status"])
  )

  # Build these dataframes if we have Z and X_S
  # ts=xs, but i dont feel like changing the naming
  Z_df <- if (!is.null(Z)) data.frame(Z = Z) else NULL
  ts_df <- if (!is.null(x_s)) {
    data.frame(
      XS = as.numeric(x_s[, "time"]),
      DL = as.numeric(x_s[, "status"])
    )
  } else NULL

  # Format newdata if present
  newdata_formatted <- if (!is.null(newdata)) {
    data.frame(
      XL = newdata[, names_list[["x_l_name"]], drop = TRUE],
      DL = newdata[, names_list[["d_l_name"]], drop = TRUE]
    )
  } else NULL

  newdata_ts_df <- if (!is.null(newdata) && !is.null(x_s)) {
    data.frame(
      XS = newdata[, names_list[["x_s_name"]], drop = TRUE],
      DL = newdata[, names_list[["d_s_name"]], drop = TRUE]
    )
  } else NULL

  newdata_Z_df <- if (!is.null(newdata) && !is.null(Z)) {
    data.frame(Z = newdata[, names_list[["covariates"]], drop = TRUE])
  } else NULL

  # Get the result based if we have X_s or Z, etc...
  # Optionally pass in newdata if we do not have it.
  result <- if (is.null(x_s) && is.null(Z)) {
    do.call(Prob.Null, list(t0, tau, formatted_data,
                            if (!is.null(newdata_formatted)) list(newdata = newdata_formatted) else list()))
  } else if (is.null(x_s)) {
    do.call(Prob.Covariate, c(list(
      t0, tau, cbind(formatted_data, Z_df), short = FALSE),
      if (!is.null(newdata)) list(newdata = cbind(newdata_formatted, newdata_Z_df)) else list()))
  } else {
    do.call(Prob.Covariate.ShortEvent, c(list(
      t0, tau, cbind(formatted_data, ts_df, Z_df)),
      if (!is.null(newdata)) list(newdata = cbind(newdata_formatted, newdata_ts_df, newdata_Z_df)) else list()))
  }

  return(result)
}

