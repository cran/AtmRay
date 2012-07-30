\name{A2P.lin}
\alias{A2P.lin}
\alias{P2A.lin}
\title{Convert Between Incidence Angle and Ray Parameter}
\description{Given an atmosphere and elevation, converts incidence angle
to ray parameter, or vice-versa.}
\usage{
A2P.lin(angle, z, az, ATM)
P2A.lin(p, z, az, ATM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Ray parameter (s/m)}
  \item{angle}{Incidence angle (from vertical) in degrees}
  \item{z}{Elevation (m)}
  \item{az}{Compass azimuth (degrees)}
  \item{ATM}{Linear atmosphere}
}
\value{Ray parameter (s/m) or incidence angle (degrees).
  }
\author{Jake Anderson}

\examples{
ATM = CheckAtm.lin()
A2P.lin(50, 50, 0, ATM)
P2A.lin(0.002321347, 50, 0, ATM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

