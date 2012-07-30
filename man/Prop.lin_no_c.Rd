\name{Prop.lin_no_c}
\alias{Prop.lin_no_c}
\title{Acoustic Ray Tracing}
\description{Given a ray parameter, azimuth, source/receiver
  elevations, calculates where ray lands at receiver elevation.
  Does not use compiled C code, so is slower, but in case of
  some bug in the C code, this function will still work.
}
\usage{
Prop.lin_no_c(p, az, zs, zr, ATM = CheckAtm.lin(list())) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ray parameter (s/m)}
  \item{az}{azimuth (degrees)}
  \item{zs}{source elevation (m)}
  \item{zr}{receiver elevation (m)}
  \item{ATM}{Linear atmosphere}
}

\value{List including the following elements:
  \item{x, y}{ending position of ray}
  \item{t}{arrival time of ray}
  \item{A}{arrival amplitude}
  \item{p}{ray parameter}
  }

\author{Jake Anderson}
\references{
Garces, M.A., Hansen, R.A., Lindquist, K.G., 1998.  Traveltimes for infrasonic waves propagating in a stratified atmosphere.  Geophysical Journal International 135, 255-263.
}
\examples{
ATM = CheckAtm.lin()
Prop.lin_no_c(0.001, 45, 100, 0, ATM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

