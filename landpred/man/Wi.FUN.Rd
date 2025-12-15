\name{Wi.FUN}
\alias{Wi.FUN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computes the inverse probability of censoring weights for a specific t0 and tau
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes the inverse probability of censoring weights for a specific t0 and tau i.e. this computes I(t0 < XL < t0+tau)*DL/G(XL) + I(XL>t0+tau)/G(t0+tau) where XL = min(TL, C), TL is the time of the long term event, C is the censoring time, DL =1*(TL<C) and G() is the estimate survival probability for censoring estimated using the Kaplan Meier estimator

}
\usage{
Wi.FUN(data, t0, tau, weight.given = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
n by k matrix, where k>= 2. A data matrix where the first column is XL = min(TL, C) where TL is the time of the long term event, C is the censoring time, and the second column is DL =1*(TL<C)
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
the landmark time..
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
the residual survival time for which probabilities are calculated. 
}
  \item{weight.given}{
%%     ~~Describe \code{weight.given} here~~
 an optional weight to be incorporated in estimation of this weight

}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Inverse probability of censoring weight.
}
\author{
%%  ~~who you are~~
Layla Parast
}

\examples{
data(data_example_landpred)
t0=2
tau = 8

W2i <- Wi.FUN(data_example_landpred[,1],data = data_example_landpred[,c(1:2)],t0=t0,tau=tau)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IPCW }
\keyword{ survival}% __ONLY ONE__ keyword per line
