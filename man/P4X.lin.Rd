\name{P4X.lin}
\alias{P4X.lin}
\title{Find Ray Parameter for Distance}
\description{Calculates ray parameter and azimuth of a ray that reaches
  some given point, using a 2-D iterative optimization search.
}
\usage{
P4X.lin(x, y, zs, zr, ATM = CheckAtm.lin(list()), maxerror = 3) 
}

\arguments{
  \item{x}{Abscissa of receiver relative to source (m)}
  \item{y}{Ordinate of receiver relative to source (m)}
  \item{zs}{Elevation of source (m)}
  \item{zr}{Elevation of receiver (m)}
  \item{ATM}{Linear atmosphere}
  \item{maxerror}{Maximum permitted error for solution (m)}
}
\details{
  'maxerror' should be set with the necessary arrival time precision in
  mind.  Smaller values of maxerror will produce more accurate arrival
  locations (and arrival times), but will require longer calculation times.
}
\value{List with following elements:
  \item{p}{ray parameter}
  \item{az}{azimuth}
  \item{error}{distance between inputs (x,y) and final position of
  result ray}
  }
\author{Jake Anderson}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
ATM = CheckAtm.lin()
P4X.lin(100, 100, 100, 00, ATM)
P4X.lin(100, 100, 100, 00, ATM, 0.01)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

