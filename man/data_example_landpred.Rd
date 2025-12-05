\name{data_example_landpred}
\alias{data_example_landpred}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Hypothetical data to be used in examples.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Hypothetical data to be used in examples.
}
\usage{data(data_example_landpred)}
\format{
  A data frame with 4868 observations on the following 5 variables.
  \describe{
    \item{\code{XL}}{a numeric vector. XL = min(TL, C) where TL is the time of the long term event, C is the censoring time. }
    \item{\code{DL}}{a 0/1 vector. DL =1*(TL<C) where TL is the time of the long term event, C is the censoring time.}
    \item{\code{XS}}{a numeric vector. XS = min(TS, C) where TS is the time of the long term event, C is the censoring time. }
    \item{\code{DS}}{a 0/1 vector. DS =1*(TS<C) where TS is the time of the long term event, C is the censoring time.}
    \item{\code{Z}}{a 0/1 vector of discrete covariate values. }
  }
}
\examples{
data(data_example_landpred)
}
\keyword{datasets}
