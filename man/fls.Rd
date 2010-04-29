\name{fls}
\Rdversion{1.1}
\alias{fls}
\title{Flexible Least Squares}
\description{
The Kalaba-Tesfatsion Flexible Least Squares method. 
Numbered equations refer to the paper:
R. Kalaba and L. Tesfatsion, Time-Varying Linear Regression Via 
Flexible Least Squares, Computers and Mathematics with Applications, 
Vol. 17 (1989), pp. 1215-1245.
}
\usage{
fls(A, b, mu = 1, ncap = length(b))
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
  M <- array (0,c(n,n,ncap))
  E <- array (0,c(n,ncap))
  X <- array (0,c(n,ncap))
  R <- matrix(0,n,n)
  diag(R) <- diag(R) + mu
  for (j in 1:ncap) {
    Z <- solve(qr(R + tcrossprod(A[j,]),LAPACK=TRUE),diag(1.0,n));
    M[,,j] <- mu*Z             # (5.7b)
    v <- b[j]*A[j,]
    if(j==1) p <- rep(0,n)
    else p <- mu*E[,j-1]
    w <- p + v
    E[,j] <- Z \%*\% w           # (5.7c)
    R <- -mu*mu*Z
    diag(R) <- diag(R) + 2*mu
  }
# Calculate eqn (5.15) FLS estimate at ncap
  Q <- -mu*M[,,ncap-1]
  diag(Q) <- diag(Q) + mu
  Ancap <- A[ncap,,drop=FALSE]
  C <- Q + t(Ancap) \%*\% Ancap
  d <- mu*E[,ncap-1,drop=FALSE] + b[ncap]*t(Ancap)
  X[,ncap] <- C \%*\% d
  X[,ncap] <- solve(qr(C,LAPACK=TRUE),d)
# Use eqn (5.16) to obtain smoothed FLS estimates for 
# X[,1], X[,2], ..., X[,ncap-1]
  for (j in 1:(ncap-1)) {
    l <- ncap - j
    X[,l] <- E[,l] + M[,,l] \%*\% X[,l+1]
  }
  X
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
