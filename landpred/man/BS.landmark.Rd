\name{BS.landmark}
\alias{BS.landmark}
\title{
Estimates the Brier score.
}
\description{

This function calculates the Brier score given the data (truth) and corresponding estimated probabilities.
}
\usage{
BS.landmark(t0, tau, data, weight=NULL)
}
\arguments{
  \item{t0}{
The landmark time.
}
  \item{tau}{
The prediction window.}
  \item{data}{
A data matrix where the first column is XL = min(TL, C) where TL is the time of the long term event, C is the censoring time; the second column is DL =1*(TL<C); the third column is the estimated probability P(TL<t0+tau | TL>t0). }
  \item{weight}{
 an optional weight to be incorporated in all estimation.
}
 }
\value{

\item{Brier.score}{Estimated Brier score}
}
\references{Parast, Layla, Su-Chun Cheng, and Tianxi Cai. Incorporating short-term outcome information to predict long-term survival with discrete markers. Biometrical Journal 53.2 (2011): 294-307.
}

\examples{
data(data_example_landpred)
t0=2
tau = 8

out = Prob.Null(t0=t0,tau=tau,data=data_example_landpred)

#get data with predictions
data_pred = out$data

#calculate training BS
BS.landmark(t0=t0,tau=tau, data = data_pred[,c("XL","DL","prob_est")])
}
\keyword{ prediction accuracy}
\keyword{ survival }
