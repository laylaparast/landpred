# Internal class to hold landpred information for discrete setting
# mode can be either no-covariate, single-covariate, or short-covariate
new_landpred_result_discrete <- function(
    mode, Prob=numeric(), data=numeric(), newdata=NA,
    t0=numeric(), tau=numeric()
    ) {

  # Prob can be null, as short-covariate mode does not return it currently
  stopifnot(is.numeric(t0))
  stopifnot(is.numeric(tau))

  structure(
    list(
      Prob=Prob,
      data=data,
      newdata=newdata,
      t0=t0,
      tau=tau,
      mode=mode
    ),
    class = "landpred_result_discrete"
  )
}
