\name{MakeArrivals.lin}
\alias{MakeArrivals.lin}
\title{Calculate Arrival Function}
\description{Calculates arrival time and arrival amplitude of many
  sources to produce an arrival function.}
\usage{
MakeArrivals.lin(xs, ys, zs, xr, yr, zr, dt, nt, timing, ATM = CheckAtm.lin(list()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xs}{source location abscissae (m)}
  \item{ys}{source location ordinates (m)}
  \item{zs}{source location elevations (m)}
  \item{xr}{receiver location abscissae (m)}
  \item{yr}{receiver location ordinates (m)}
  \item{zr}{receiver location elevations (m)}
  \item{dt}{time interval (s)}
  \item{nt}{number of time steps to evaluate}
  \item{timing}{source times (s)}
  \item{ATM}{linear atmosphere}
}
\details{
  This is suitable for calculating arrival functions for spatially
  distributed acoustic sources.  Each of sx, sy, and sz should be the
  same length, and each of rx, ry, and rz should be the same length.
}
\value{
  Returns a matrix of dimension nt x length(rx), with rows corresponding
  to time steps and columns to receivers.
}

\author{Jake Anderson}

\examples{
MakeArrivals.lin(xs = c(100, 150), ys = c(100, 150), zs = c(100, 150), xr = 0, yr = 0, zr = 0,
dt = 0.01, nt = 100, timing = c(0, 0), ATM = CheckAtm.lin())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

