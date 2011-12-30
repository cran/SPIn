\name{plot.SPIn}
\alias{plot.SPIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot the Results by SPIn
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot the histogram, the kernel density estimate, the shortest probability interval and the central interval.
}
\usage{
\method{plot}{SPIn}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
SPIn object, result of \code{\link{SPIn}}.
}
  \item{...}{
See \code{plot}.
}
}

\references{
%% ~put references to the literature/web site here ~
An Optimal Weighting Strategy for Better Shortest Probability Intervals. (In preparation)
}
\author{
%%  ~~who you are~~
Ying Liu yliu@stat.columbia.edu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{SPIn}}
}
\examples{
x <- rgamma(100,3)
r <- SPIn(x)
plot(r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
