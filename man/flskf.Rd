\name{flskf}
\Rdversion{1.1}
\alias{flskf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
flskf(A, b, mu = 1, ncap = length(b))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{ncap}{
%%     ~~Describe \code{ncap} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, b, mu=1, ncap=length(b))
{
  m <- nrow (A)
  n <- ncol (A)
  X <- array (0,c(n,ncap))
  R <- array (0,c(n,n))
  diag (R) <- 1/mu
  a <- t(A[1,,drop=FALSE])
  phi <- c((1/mu) * crossprod(a) + 1)
  w <- (1/phi) * R \%*\% a
  rho <- b[1]
  X[,1] <- rho * w
  for (j in 2:ncap) {
    R = R - phi * tcrossprod(w)
    diag (R) <- diag (R) + 1/mu
    a <- t(A[j,,drop=FALSE])
    phi <- c(t(a) \%*\% R \%*\% a)
    phi <- phi + 1
    rho <- c(b[j] - crossprod (a, X[,j-1]))
    w <- (1/phi) * R \%*\% a
    X[,j] <- X[,j-1] + rho*w
  }
  X
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
