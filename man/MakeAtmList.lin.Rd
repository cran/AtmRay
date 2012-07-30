\name{MakeAtmList.lin}
\alias{MakeAtmList.lin}
\title{Make List of Linear Atmospheres}
\description{Given a set of values that must be taken by each field in a
linear atmosphere, returns a list of atmospheres that span all
combinations of those values (similar to meshgrid).}
\usage{
MakeAtmList.lin(z0 = 0, c0 = 343, gc = 0, wx0 = 0, gwx = 0, wy0 = 0, gwy = 0, rho0 = 1.2929 * exp(-z0/6800), grho = -0.0001901323 * exp(-z0/6800))
}
\arguments{
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
\details{
  This is useful for modeling wave propagation in a range of atmospheres
  when actual atmospheric characteristics are poorly constrained.
}
\value{List of linear atmospheres.
  }
\author{Jake Anderson}
\examples{
# make atmospheres spanning a range of base sound speeds and
# sound speed gradients
c0 = seq(330, 336, 0.1) 
gc = seq(-0.006, -0.004, 0.0005)

ATM_list = MakeAtmList.lin(c0 = c0, gc = gc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

