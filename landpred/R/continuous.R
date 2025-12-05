# ------ Legacy functions ----------------
Kern.FUN.CONT <- function(zz,zi,bw) ## returns an (n x nz) matrix
{
  out = (VTM.CONT(zz,length(zi))- zi)/bw
  norm.k = dnorm(out)/bw
  norm.k
}

g.logit <- function(xx){exp(xx)/(1+exp(xx))}

VTM.CONT <-
  function(vec, rown) {
    matrix(vec,
           nrow = rown,
           ncol = length(vec),
           byrow = T)
  }

Wi.FUN.CONT <- function(tt, data, t0, tau, weight.given = NULL)	{
  Xi <- data[,1]; Di <- data[,2]; wihat <- rep(0, length(Xi))
  tmpind1 <- (Xi > t0)&(Xi <= (t0+tau)); tt0 <- c(Xi[tmpind1], t0 + tau); Ghat.tt0 <- Ghat.FUN(tt0,data, weight.given=weight.given)
  wihat[Xi > (t0+tau)] <- 1/Ghat.tt0[length(tt0)]
  wihat[tmpind1] <- Di[tmpind1]/Ghat.tt0[-length(tt0)]
  wihat
}

loc.fun.ex.CONT <- function(t, data,tau, s, h, weight = NULL, transform=identity) {
  X1i = data[,1]; X2i = data[,2]; D1i = data[,3]; D2i = data[,4]; Zi = data[,-c(1:4)]
  if (is.null(weight))	{
    W2i <- Wi.FUN.CONT(X2i,
                       data = cbind(X2i, D2i, Zi),
                       t0 = tau,
                       tau = s)
  } else{
    W2i <- weight
  }
  index.sub = data[,2] > tau & (transform(data[,1]) < transform(tau)) & (data[,3] == 1)
  K = Kern.FUN.CONT(transform(X1i[index.sub]),t,h)
  est.mat = matrix(nrow=length(t), ncol = dim(as.matrix(Zi))[2] + 1)
  invinf.mat = matrix(nrow=length(t), ncol = (dim(as.matrix(Zi))[2] + 1)^2)
  tmpfm <- paste("1*(X2i< tau+s) ~ ", paste(names(Zi),collapse="+"))
  for(i in 1:length(t)) {
    w_final <- K[i,] * W2i[index.sub]
    
    # filter valid weights to avoid warnings
    valid_weights <- w_final > 0
    
    if(sum(valid_weights) > 0) {
      m = glm(as.formula(tmpfm), 
              data=data[index.sub, , drop=FALSE][valid_weights, , drop=FALSE],  
              family = "quasibinomial", 
              weights = w_final[valid_weights])
      est.mat[i,] = m$coeff
      invinf.mat[i,] = as.vector(vcov(m))
    } else {
      est.mat[i,] = NA
      invinf.mat[i,] = NA
    }
  }
  return(list("est.mat" = est.mat, "invinf.mat" = invinf.mat))
}

# --------------------------------------------

subset_and_format_df <- function(landpred_obj, indexes, transform=identity) {
  X_L <- landpred_obj$X_L[, "time"]
  X_S <- transform(landpred_obj$X_S[, "time"])
  X_D <- landpred_obj$X_L[, "status"]
  Z <- landpred_obj$Z

  subsetted_data <- cbind(data.frame(
    X_L=X_L[indexes],
    X_D=X_D[indexes],
    X_S=X_S[indexes]
  ), Z[indexes, , drop=FALSE])

  subsetted_data
}

# Subset of the data used to fit model with no info on short covariate
get_subset_indexes_noinfo <- function(landpred_obj, t0) {
  if (!is.null(landpred_obj$X_S)) {
    subset_indexes <- pmin(landpred_obj$X_L[, "time"], landpred_obj$X_S[, "time"]) > t0
  } else {
    subset_indexes <- landpred_obj$X_L[, "time"] > t0
  }
  subset_indexes
}

get_subset_indexes_noinfo <- function(landpred_obj, t0) {
  # Should match old logic: X_L > t0 (not pmin)
  subset_indexes <- landpred_obj$X_L[, "time"] > t0
  return(subset_indexes)
}

# Subset of the data used to fit model with information on short covariate
# maybe have to check if it was also censored?
get_subset_indexes_short <- function(landpred_obj, t0, indices=NULL, transform=identity) {
  subset_indices <- landpred_obj$X_L[, "time"] > t0 &
    (transform(landpred_obj$X_S[, "time"]) < transform(t0)) & (landpred_obj$X_S[, "status"] == 1)
  subset_indices
}

# Linear model for t0 + tau where our given t_s < t0
# no information on short covariate
#' Fit GLM with Normal Weights (No Short Covariate Info)
#'
#' Fits a GLM for the probability of the event occurring before \code{t0 + tau},
#' given survival up to \code{t0}, using only baseline covariates.
#'
#' @param landpred_obj A landpred object containing the data.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#'
#' @return A fitted glm object.
#' @return A fitted glm object.
#' @keywords internal
fit_glm_normal <-
  function(landpred_obj, t0, tau) {
    ts_formula <-
      paste("I(X_L < t0 + tau) ~", paste(names(landpred_obj$Z), collapse = "+"))

    weights <- Wi.FUN.CONT(
      landpred_obj$X_L[, "time"],
      data = cbind(
        landpred_obj$X_L[, "time"],
        landpred_obj$X_L[, "status"]
      ),
      t0 = t0,
      tau = tau
    )

    ts_gt_subset <- get_subset_indexes_noinfo(landpred_obj, t0)
    subset_gt_t0 <- subset_and_format_df(landpred_obj, ts_gt_subset)

    # fit as quasibinomial to suppress warnings about non-integer successes
    model <- glm(
      as.formula(ts_formula),
      data = subset_gt_t0,
      family = "quasibinomial",
      weights = weights[ts_gt_subset]
    )

    model
  }

# fit glm with kernel weights and information on short covariate
#' Fit GLM with Kernel Weights (Short Covariate Info)
#'
#' Fits a GLM for the probability of the event occurring before \code{t0 + tau},
#' given survival up to \code{t0} and information on a short-term covariate.
#' Uses kernel weighting based on the short-term covariate value.
#'
#' @param landpred_obj A landpred object containing the data.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#' @param t_s The time of the short-term covariate measurement.
#' @param bw The bandwidth for kernel weighting.
#' @param transform A transformation function for the time variable (e.g., log).
#' @param indices Optional indices to subset the data.
#'
#' @return A fitted glm object.
#' @return A fitted glm object.
#' @keywords internal
fit_short_glm <-
  function(landpred_obj, t0, tau, t_s, bw, transform, indices=NULL) {
    ts_formula <-
      paste("I(X_L < t0 + tau) ~", paste(names(landpred_obj$Z), collapse = "+"))

    weights <- Wi.FUN.CONT(
      landpred_obj$X_L[, "time"],
      data = cbind(
        landpred_obj$X_L[, "time"],
        landpred_obj$X_L[, "status"]
      ),
      t0 = t0,
      tau = tau
    )

    ts_lt_subset <- get_subset_indexes_short(landpred_obj, t0, transform=transform)

    if (!is.null(indices)) {
      ts_lt_subset <- intersect(ts_lt_subset, indices)
    }

    # get subset and apply transformation
    subset_lt_t0 <- subset_and_format_df(landpred_obj, ts_lt_subset, transform)
    K = Kern.FUN.CONT(subset_lt_t0$X_S, t_s, bw)[1, ]

    w_final <- weights[ts_lt_subset] * K
    valid_weights <- w_final > 0
    
    if(sum(valid_weights) > 0) {
      model <- glm(
        as.formula(ts_formula),
        data = subset_lt_t0[valid_weights, , drop=FALSE],
        family = "quasibinomial",
        weights = w_final[valid_weights]
      )
    } else {
      # when all weights are zero none of the observations will be used
      warning("All weights are zero in fit_short_glm. Returning NULL model.")
      return(NULL)
    }
    model
  }


handle_continuous_pred <- function(model, newdata) {
  landpred_obj <- model$landpred_obj
  t0 <- model$t0
  tau <- model$tau
  bw <- model$bw
  glm_noinfo <- model$glm_noinfo
  transform <- model$transform

  formatted_data <- newdata[, landpred_obj$names[["covariates"]], drop=FALSE]


  X_S <- newdata[, landpred_obj$name[["x_s_name"]], drop=TRUE]
  t_s_values <- unique(X_S)

  response <- numeric(0)

  for(s in t_s_values) {
    data_subset <- formatted_data[X_S==s, , drop=FALSE]

    model_specified <- NULL
    if(s < t0)  {
      model_specified <- glm_noinfo
    } else {
      glm_shortinfo <- fit_short_glm(
        landpred_obj, t0, tau, s, bw, transform
      )
      model_specified <- glm_shortinfo
    }

    preds <- predict(
      model_specified, newdata=data_subset,
      type="response"
    )

    response[X_S == s] <- preds
  }

 response
}

# confidence interval for short info model
# copied from legacy version but with support for any number of covariates
#' Estimate Variance of Coefficients
#'
#' Estimates the variance of the coefficients for the short-term GLM using
#' a perturbation resampling method.
#'
#' @param t The target time for the short-term covariate (usually t_s).
#' @param data.v The data frame used for estimation.
#' @param tau The landmark time.
#' @param s The prediction window.
#' @param h The bandwidth.
#' @param vmat A matrix of perturbation weights.
#' @param Ainv The inverse information matrix from the model fit.
#' @param weight Optional weights.
#' @param transform Transformation function.
#'
#' @return A list containing the estimated standard errors for the intercept and slopes.
#' @return A list containing the estimated standard errors for the intercept and slopes.
#' @keywords internal
var.fun <- function(t, data.v, tau, s, h, vmat, Ainv, weight = NULL, transform=identity) {
  length.t <- length(t)
  size <- nrow(data.v)

  X1i <- data.v[, 1]
  X2i <- data.v[, 2]
  D1i <- data.v[, 3]
  D2i <- data.v[, 4]
  Zi  <- data.v[, -c(1:4), drop = FALSE]

  index.sub <- which(
    data.v[, 2] > tau &
      transform(data.v[, 1]) < transform(tau) &
      data.v[, 3] == 1
  )

  # Compute weights
  if (is.null(weight)) {
    W2i <- Wi.FUN.CONT(
      X2i,
      data = cbind(X2i, D2i, Zi),
      t0 = tau,
      tau = s
    )
  } else {
    W2i <- weight
  }

  # Local estimation at data points for residuals
  loc.m <- loc.fun.ex.CONT(
    t   = transform(X1i[index.sub]),
    data = data.v,
    tau = tau,
    s   = s,
    h   = h,
    transform = transform
  )

  # Adjusted vmat
  vmat.c <- vmat - 1

  # Difference term
  # Calculate linear predictor row-wise
  linear_pred <- rowSums(loc.m$est.mat * cbind(1, Zi[index.sub, , drop = FALSE]))
  
  diff <- 1 * (X2i > tau & X2i < tau + s & D2i == 1)[index.sub] - g.logit(linear_pred)

  # Remove NA entries
  valid_idx <- which(!is.na(diff))
  index.sub <- index.sub[valid_idx]
  diff <- diff[valid_idx]

  # Kernel
  K <- Kern.FUN.CONT(transform(X1i[index.sub]), t, h)

  # Piece 1
  piece.1.int <- t(t(K) * (W2i[index.sub] * diff)) %*% t(vmat.c[, index.sub])

  # Piece 2
  W2i.star <- apply(
    vmat,
    1,
    Wi.FUN.CONT,
    tt = X2i,
    data = cbind(X2i, D2i, Zi),
    t0 = tau,
    tau = s
  )

  piece.2.int <- t(t(K) * diff) %*% (
    W2i.star[index.sub, ] -
      matrix(W2i[index.sub],
             ncol = ncol(W2i.star),
             nrow = length(index.sub),
             byrow = FALSE)
  )

  piece.3.int <- piece.1.int + piece.2.int

  # Number of covariates
  len.z <- ncol(Zi)
  total.params <- len.z + 1  # intercept + slopes

  # Build slope matrices
  if (len.z > 0) {
    slope.mat <- matrix(nrow = length.t * len.z, ncol = ncol(piece.3.int))
    for (j in seq_len(len.z)) {
      piece.1.slope <- t(t(K) * (W2i[index.sub] * Zi[index.sub, j] * diff)) %*% t(vmat.c[, index.sub])
      piece.2.slope <- t(t(K) * (Zi[index.sub, j] * diff)) %*% (
        W2i.star[index.sub, ] -
          matrix(W2i[index.sub],
                 ncol = ncol(W2i.star),
                 nrow = length(index.sub),
                 byrow = FALSE)
      )
      slope.mat[((j - 1) * length.t + 1):(j * length.t), ] <- piece.1.slope + piece.2.slope
    }
  } else {
    slope.mat <- matrix(0, nrow = 0, ncol = ncol(piece.3.int))
  }

  # Output matrices
  o.matrix <- matrix(NA_real_, nrow = length.t, ncol = total.params)
  o.entire <- vector("list", total.params)
  for (i in seq_len(total.params)) {
    o.entire[[i]] <- matrix(NA_real_, nrow = length.t, ncol = ncol(piece.3.int))
  }

  # Main computation loop
  for (l in seq_len(length.t)) {
    a <- matrix(Ainv[l, ], total.params, total.params)

    for (i in seq_len(total.params)) {
      o_i <- piece.3.int[l, ] * a[i, 1]
      if (len.z > 0) {
        for (j in seq_len(len.z)) {
          slope_row <- slope.mat[l + length.t * (j - 1), ]
          o_i <- o_i + slope_row * a[i, j + 1]
        }
      }
      o.entire[[i]][l, ] <- o_i
      o.matrix[l, i] <- sd(o_i)
    }
  }

  sl_mat <-
      o.matrix[, -1, drop = FALSE]

  out <- list(
    "int" = o.matrix[, 1],
    "sl"  = sl_mat
  )

  for (i in seq_len(total.params)) {
    out[[paste0("o", i)]] <- o.entire[[i]]
  }

  return(out)
}

#' Calculate Standard Errors for Coefficients
#'
#' Calculates standard errors for the coefficients of the landpred model.
#' If \code{t_s} is provided, it uses the perturbation resampling method.
#' Otherwise, it returns the standard errors from the GLM.
#'
#' @param model A landpred_model_continuous object.
#' @param t_s The time of the short-term covariate measurement.
#' @param samples The number of resampling iterations.
#'
#' @return A named vector of standard errors.
#' @return A named vector of standard errors.
#' @keywords internal
coefficient_se <-
  function(model,
           t_s=NULL,
           samples = 200
           ) {
    std_errors <- c()

    # If no t_s, use standard errors calculated by glm,
    # otherwise do special pertubation resampling method
    if (is.null(t_s) || t_s > model$t0) {
      std_errors <- summary(model$glm_noinfo)$coefficients[, "Std. Error"]
      std_errors <- as.numeric(std_errors)
    } else {
      # legacy formatted data
      df_formatted <-
        landpred_to_legacy_data(model$landpred_obj)

      fit = loc.fun.ex.CONT(
        t = t_s,
        data = df_formatted,
        tau = model$t0,
        s = model$tau,
        h = model$bw,
        transform = model$transform
      )

     size <- dim(df_formatted)[1]
     vmat <- matrix(rexp(samples * size, 1), nrow = samples, ncol = size)
     var_results <- var.fun(
       t = t_s,
       data.v = df_formatted,
       tau = model$t0,
       s = model$tau,
       h = model$bw,
       vmat = vmat,
       Ainv = fit$invinf,
       transform = model$transform
     )

      std_errors <- c(var_results$int, var_results$sl)
    }

    setNames(std_errors, names(coef(model$glm_noinfo)))
}

#' Get Landpred Model
#'
#' Fits the base GLM (no short covariate info) and creates a landpred model object.
#'
#' @param landpred_obj A landpred object.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#' @param bw The bandwidth.
#' @param transform Transformation function.
#'
#' @return A landpred_model_continuous object.
#' @export
get_model <- function(landpred_obj, t0, tau, bw, transform=identity) {
  glm_noinfo <- fit_glm_normal(landpred_obj, t0, tau)
  new_landpred_model_continuous(
    landpred_obj, glm_noinfo, t0, tau, bw,
    transform
  )
}

new_landpred_model_continuous <- function(landpred_obj, glm_noinfo, t0, tau, bw, transform) {
  structure(
    list(
      landpred_obj=landpred_obj,
      glm_noinfo=glm_noinfo,
      t0=t0,
      tau=tau,
      bw=bw,
      transform=transform
    ),
    class = "landpred_model_continuous"
  )
}

#' Predict Method for Landpred Continuous Model
#'
#' Predicts the probability of the event occurring given new data.
#'
#' @param object A landpred_model_continuous object.
#' @param newdata New data frame containing covariates and short-term event info.
#' @param type Type of prediction (default "response").
#' @param ... Additional arguments 
#'
#' @return A vector of predicted probabilities.
#' @export
predict.landpred_model_continuous <- function(object, newdata=NULL, type="response", ...) {
  handle_continuous_pred(object, newdata)
}

#' Extract Coefficients from Landpred Continuous Model
#'
#' Extracts coefficients. If \code{t_s} is provided, it fits the short-term GLM
#' and returns its coefficients.
#'
#' @param object A landpred_model_continuous object.
#' @param t_s Optional short-term covariate time.
#' @param ... Additional arguments.
#'
#' @return A named vector of coefficients.
#' @export
coef.landpred_model_continuous <- function(object, t_s=NULL, ...) {
  if(is.null(t_s) || t_s > object$t0) {
    return(coef(object$glm_noinfo))
  } else {
    glm_shortinfo <- fit_short_glm(
      object$landpred_obj, object$t0, object$tau, t_s, object$bw, object$transform
    )
    return(coef(glm_shortinfo))
  }
}

#' Print Method for Landpred Continuous Model
#'
#' Prints the continuous landpred model results.
#'
#' @param x A landpred_model_continuous object.
#' @param ... Additional arguments.
#'
#' @export
print.landpred_model_continuous <- function(x, ...) {
  cat("\nContinuous Landpred Model:\n")
  cat("Call coef() or summary() with t_s to see coefficients.\n")
  print(x$landpred_obj)
  cat("\n")
  cat(sprintf("t0: %-10.3f tau: %-10.3f", x$t0, x$tau))
}

#' Summary Method for Landpred Continuous Model
#'
#' Prints a summary of the model, including coefficients and standard errors.
#'
#' @param object A landpred_model_continuous object.
#' @param t_s Optional short-term covariate time.
#' @param ... Additional arguments.
#'
#' @export
summary.landpred_model_continuous <- function(object, t_s = NULL, ...) {
  cat("\nContinuous Landpred Model:\n\n")

  t_s_message <- if (is.null(t_s)) "No short covariate" else sprintf("t_s=%f", t_s)
  cat(sprintf("Coefficients (%s):\n", t_s_message))

  # Choose model
  glm_fit <-
    if (is.null(t_s) ||
        t_s > object$t0)
      object$glm_noinfo
  else
    fit_short_glm(object$landpred_obj, object$t0, object$tau, t_s, object$bw, object$transform)

  coef_est <- coef(glm_fit)

  se <- coefficient_se(object, t_s)

  coef_table <- data.frame(
    Estimate = coef_est,
    `Std. Error` = se,
    check.names = FALSE
  )

  printCoefmat(coef_table, P.values = FALSE, has.Pvalue = FALSE)
  cat("---\n")
  cat(sprintf("Fit on n=%d observations.\n", nobs(glm_fit)))

  cat(sprintf("\nt0: %-10.3f tau: %-10.3f\n", object$t0, object$tau))
}

landpred_to_legacy_data <- function(landpred_obj) {
  # Extract components
  X_L_time <- landpred_obj$X_L[, "time"]
  X_L_status <- landpred_obj$X_L[, "status"]
  X_S_time <- landpred_obj$X_S[, "time"]
  X_S_status <- landpred_obj$X_S[, "status"]
  Z <- landpred_obj$Z

  # Combine into legacy format: [X1i, X2i, D1i, D2i, Zi...]
  legacy_data <- cbind(
    X2i = X_S_time,     # Short covariate time
    X1i = X_L_time,     # Long covariate time
    D2i = X_S_status,   # Short covariate status
    D1i = X_L_status,   # Long covariate status
    Z                   # Additional covariates
  )

  return(legacy_data)
}

# Updated legacy bandwidth selection function using w_i instead of Wi.FUN
#' @keywords internal
min_BW_cens_ex <- function(data.cv, tau, s, h) {
  X1i = data.cv[, 1]
  X2i = data.cv[, 2]
  D1i = data.cv[, 3]
  D2i = data.cv[, 4]
  Zi = data.cv[, -c(1:4)]

  W2i <- Wi.FUN(data.cv[,2],data = cbind(data.cv[,2],data.cv[,4],data.cv[,6]),t0=tau,tau=s)

  observ = dim(data.cv)[1]
  nv = floor(observ / 3)
  BW.vec = vector(length = 3)
  replic = matrix(nrow = 100, ncol = 3)

  for (lay in 1:100) {
    for (k in 1:3) {
      ind.val = sample(1:observ, nv)
      ind.tra = setdiff(1:observ, ind.val)
      subset.val = ind.val[D1i[ind.val] == 1 &
                             X1i[ind.val] < log(tau) & X2i[ind.val] > tau]
      subset.tra = ind.tra[D1i[ind.tra] == 1 &
                             X1i[ind.tra] < log(tau) & X2i[ind.tra] > tau]

      # Create temporary landpred object for training data
      train_landpred <- list(
        X_L = cbind(time = X2i[subset.tra], status = D2i[subset.tra]),
        X_S = cbind(time = X1i[subset.tra], status = D1i[subset.tra]),
        Z = Zi[subset.tra, , drop = FALSE],
        names = list(covariates = colnames(Zi), x_s_name = "X_S_time")
      )
      class(train_landpred) <- "landpred_obj"

      # Use fit_short_glm instead of loc.fun.ex
      loc.model = fit_short_glm(
        landpred_obj = train_landpred,
        t0 = tau,
        tau = s,
        t_s = X1i[subset.val],  # validation time points
        bw = h,
        transform = identity
      )

      if (is.null(loc.model)) {
        BW.vec[k] = NA
        next
      }

      # Create prediction data for validation set
      # We need to construct val_data from data.cv
      # data.cv has X1i (Short), X2i (Long), D1i, D2i, Zi
      # predict expects newdata with Z and X_S_time (if needed)
      
      val_data <- as.data.frame(Zi[subset.val, , drop=FALSE])
      # Add X_S_time if needed (it is X1i)
      # The model expects the name "X_S_time" or whatever was used in training?
      # train_landpred names: x_s_name = "X_S_time"
      val_data$X_S_time <- X1i[subset.val]
      
      p.hat = predict(loc.model, newdata = val_data, type = "response")
      BW = sum((((X2i[subset.val] <= tau + s) - p.hat) ^ 2) * (X2i[subset.val] >
                                                                 tau) * (X1i[subset.val] < log(tau)) * (W2i[subset.val]) * (D1i[subset.val] ==
                                                                                                                              1))
      BW.vec[k] = BW
    }
    replic[lay, ] = cbind(BW.vec[1], BW.vec[2], BW.vec[3])
  }
  mean(replic, na.rm = TRUE)
}

# Adapter to use legacy bandwidth selection with new landpred API

# Convert landpred_obj to legacy data format
landpred_to_legacy_data <- function(landpred_obj) {
  # Extract components
  X_L_time <- landpred_obj$X_L[, "time"]
  X_L_status <- landpred_obj$X_L[, "status"]
  X_S_time <- landpred_obj$X_S[, "time"]
  X_S_status <- landpred_obj$X_S[, "status"]
  Z <- landpred_obj$Z

  # Combine into legacy format: [X1i, X2i, D1i, D2i, Zi...]
  legacy_data <- cbind(
    X1i = X_S_time,     # Short covariate time
    X2i = X_L_time,     # Long covariate time
    D1i = X_S_status,   # Short covariate status
    D2i = X_L_status,   # Long covariate status
    Z                   # Additional covariates
  )

  return(legacy_data)
}


#' Calculate MSE for Bandwidth Selection using Cross-Validation
#'
#' @param bw The bandwidth to test.
#' @param landpred_obj The landpred object.
#' @param t0 The landmark time.
#' @param tau The prediction window.
#' @param transform Transformation function for short-term covariate.
#' @param reps Number of repetitions.
#' @param train_prop Proportion of data to use for training.
#' @return The Mean Squared Error.
#' @keywords internal
mse_cv <- function(bw, landpred_obj, t0, tau, transform = identity, reps = 50, train_prop = 0.66) {
  X_L_time <- landpred_obj$X_L[, "time"]
  X_L_status <- landpred_obj$X_L[, "status"]
  X_S_time <- landpred_obj$X_S[, "time"]
  X_S_status <- landpred_obj$X_S[, "status"]
  Z <- landpred_obj$Z
  
  n <- nrow(landpred_obj$X_L)
  n_train <- floor(n * train_prop)
  
  mse_sum <- 0
  valid_reps <- 0
  
  W_all <- Wi.FUN(
    data = cbind(X_L_time, X_L_status),
    t0 = t0,
    tau = tau
  )
  
  for (i in 1:reps) {
    ind_train <- sample(1:n, n_train)
    ind_val <- setdiff(1:n, ind_train)
    
    subset_train <- ind_train[X_S_status[ind_train] == 1 & 
                              X_S_time[ind_train] < t0 & 
                              X_L_time[ind_train] > t0]
                              
    subset_val <- ind_val[X_S_status[ind_val] == 1 & 
                          X_S_time[ind_val] < t0 & 
                          X_L_time[ind_val] > t0]
    
    if (length(subset_train) < 10 || length(subset_val) < 10) next
    
    train_obj <- list(
      X_L = landpred_obj$X_L[subset_train, , drop=FALSE],
      X_S = landpred_obj$X_S[subset_train, , drop=FALSE],
      Z = landpred_obj$Z[subset_train, , drop=FALSE],
      names = landpred_obj$names
    )
    class(train_obj) <- "landpred_object"
    
    t_s_val <- transform(X_S_time[subset_val])
    
    tryCatch({
      model <- fit_short_glm(
        landpred_obj = train_obj,
        t0 = t0,
        tau = tau,
        t_s = t_s_val,
        bw = bw,
        transform = transform
      )
      
      val_data <- subset_and_format_df(landpred_obj, subset_val)
      
      preds <- predict(model, newdata = val_data, type = "response")
      
      truth <- (X_L_time[subset_val] <= t0 + tau)
      weights <- W_all[subset_val]
      
      # Filter out cases where weight is 0 (censored before t0+tau)
      # Actually Wi.FUN handles this by setting weight to 0?
      # Yes.
      
      mse <- sum(weights * (truth - preds)^2) / sum(weights)
      mse_sum <- mse_sum + mse
      valid_reps <- valid_reps + 1
      
    }, error = function(e) {
      # Ignore errors in optimization
    })
  }
  
  if (valid_reps == 0) return(Inf)
  return(mse_sum / valid_reps)
}
