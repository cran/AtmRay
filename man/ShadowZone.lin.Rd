\name{ShadowZone.lin}
\alias{ShadowZone.lin}
\title{Calculate Distance to Acoustic Shadow Zone}
\description{In atmospheres where effective sound speed varies with
  elevation, rays will eventually turn and reverse their vertical
  direction.  As a result, for every source/receiver elevation pair,
  there exists some horizontal distance beyond which waves fail to
  propagate, referred to as a shadow zone.}
\usage{
ShadowZone.lin(zr, zs, az = 0, ATM)	
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zr}{receiver elevation (m)}
  \item{zs}{source elevation (m)}
  \item{az}{wave propagation azimuth (degrees)}
  \item{ATM}{linear atmosphere}
}
\details{
  zr and zs are interchangeable and not strictly the end and beginning
  of the ray.  That is, if the values of zr and zs are switched, this
  function will return the same result, even though a shadow zone only
  exists for a single configuration (when waves propagate in the
  direction of increasing sound speed).  
}
\value{List including the following elements:
  \item{X}{distance to edge of shadow zone (m)}
  \item{p}{ray parameter of ray that reaches distance X}
  }
\author{Jake Anderson}
\examples{
ATM = CheckAtm.lin(list(gc = -0.005)) # typical sound speed lapse rate
ShadowZone.lin(0, 100, 0, ATM) 
}

\keyword{ misc }

