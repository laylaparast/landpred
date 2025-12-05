\name{AUC.landmark}
\alias{AUC.landmark}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimates the area under the ROC curve (AUC).
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the AUC given the data (truth) and corresponding estimated probabilities; uses a continuity correction.
}
\usage{
AUC.landmark(t0, tau, data, short = TRUE, weight=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{
%%     ~~Describe \code{t0} here~~
The landmark time.
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
The prediction window.}
  \item{data}{
%%     ~~Describe \code{data} here~~
Matrix: n by k matrix, where k = 4 or 6. A data matrix where the first column is XL = min(TL, C) where TL is the time of the long term event, C is the censoring time, and the second column is DL =1*(TL<C), the second to last column is the covariate vector (can be NULL) and the last column is the estimated probability P(TL<t0+tau | TL>t0). }
  \item{short}{
%%     ~~Describe \code{short} here~~
Logical value indicating whether data includes short term event information. Should be TRUE if short term XS and DS are includes as third and fourth columns of data matrix, FALSE if not. Default is TRUE. 
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
 an optional weight to be incorporated in all estimation.
}
 }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{AUC.est}{Estimated AUC}
}
\references{Parast, Layla, Su-Chun Cheng, and Tianxi Cai. Incorporating short-term outcome information to predict long-term survival with discrete markers. Biometrical Journal 53.2 (2011): 294-307.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(data_example_landpred)
t0=2
tau = 8
Prob.Null(t0=t0,tau=tau,data=data_example_landpred)

out = Prob.Null(t0=t0,tau=tau,data=data_example_landpred)
out$Prob
out$data

AUC.landmark(t0=t0,tau=tau, data = out$data)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prediction accuracy}
\keyword{ survival }% __ONLY ONE__ keyword per line
