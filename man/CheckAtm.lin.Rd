\name{CheckAtm.lin}
\alias{CheckAtm.lin}
\title{Check Linear Atmosphere}
\description{Verifies that an atmosphere contains all required elements
  and fills in missing elements with default values.}
\usage{
CheckAtm.lin(ATM = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{h}{Focal depth (km)}
  \item{ATM}{List that may include characteristics of an atmosphere}
}
\details{
  A linear atmosphere variable must contain the elements described in
  the "Value" section.  Default values for these are as follows:
  z0: 0
  c0: 330
  wx0: 0
  wy0: 0
  rho0: 1.2929 * exp(-ATM$z0/6800)
  gc: -10^9
  gwx: 0
  gwy: 0
  grho: -0.0001901058 * exp(-ATM$z0/6800)

  Other functions may encounter problems when working with an effective
  sound speed gradient of zero.
  
}
\value{List including the following elements:
  \item{z0}{Elevation of intercept layer (m)}
  \item{c0}{Intrinsic sound speed at intercept layer (m/s)}
  \item{wx0}{Zonal (east-west) wind at intercept layer (m/s)}
  \item{wy0}{Meridional (north-south) wind at intercept layer (m/s)}
  \item{rho0}{Density at intercept layer (kg/m^3)}
  \item{gc}{Vertical intrinsic sound speed gradient (1/s)}
  \item{gwx}{Vertical zonal wind gradient (1/s)}
  \item{gwy}{Vertical meridional wind gradient (1/s)}
  \item{grho}{Vertical density gradient (kg/m^4)}
  }

\author{Jake Anderson}

\examples{
# quickly make a new default atmosphere
ATM = CheckAtm.lin()

# fill in missing values for an existing atmosphere
ATM = list(c0 = 343, gc = -0.006)
ATM = CheckAtm.lin(ATM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

