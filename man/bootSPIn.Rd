\name{bootSPIn}
\alias{bootSPIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
SPIn with Bootstrap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the shortest probability interval (spin) using SPIn with bootstrap.
}
\usage{
bootSPIn(x, n.boot = 50, conf = 0.95, bw = 0, lb = -Inf, ub = Inf, l = NA, u = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A vector of samples from the distribution.
}
  \item{n.boot}{
%%     ~~Describe \code{n.boot} here~~
Number of bootstraps.
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
Scalar, the confidence level desired.
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
Scalar, the bandwidth of the weighting kernel in terms of sample points. If not specified, sqrt(n) will be used, where n is the sample size.
}
  \item{lb,ub}{
%%     ~~Describe \code{lb} here~~
Scalars, the lower and upper bounds of the distribution. If specified, a pseudo-sample point equal to the corresponding bound will be added.
}
  \item{l,u}{
%%     ~~Describe \code{l} here~~
Scalars, weighting centers (if provided).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
spin.boot computes the shortest probability interval for a distribution using SPIn with bootstrap.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
spin.boot returns an object of \code{\link{class}} 'SPIn'. An object of class 'SPIn' is a list containing the following components:
\item{spin }{A vector of length 2 with the lower and upper endpoints of the interval.}
\item{conf }{The confidence level.}
\item{x }{A vector of samples from the distribution.}
\item{w.l,w.u }{Vectors of the computed weights.}
\item{l.l,l.u,u.l,u.u }{Endpoints of the weights.}
}
\references{
%% ~put references to the literature/web site here ~
Simulation efficient shortest probability intervals. (arXiv:1302.2142)
}
\author{
%%  ~~who you are~~
Ying Liu yliu@stat.columbia.edu
}
\note{
%%  ~~further notes~~
This function assumes that the distribution is unimodal, and computes only 1 interval, not the set of intervals that are appropriate for multimodal distributions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{plot.SPIn}},\code{\link{SPIn}}
}
\examples{
x <- rgamma(100,3)
bootSPIn(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
