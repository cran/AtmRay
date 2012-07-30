\name{PlotAtm.lin}
\alias{PlotAtm.lin}
\title{Plot Linear Atmosphere}
\description{Makes a plot showing effective sound speed structure.}
\usage{PlotAtm.lin(ATM, zlim = c(0, 100), winddir = 90, col =
sky.colors(500),
TOPO = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ATM}{Linear Atmosphere}
  \item{zlim}{Height limits in plot}
  \item{winddir}{Direction of wind to plot}
  \item{col}{color vector, such as the output of heat.colors}
  \item{TOPO}{list containing vectors x and y, and matrix z, or NULL}
}
\details{
  When plotting topography, an east-west cross-section where
  TOPO$y is zero is plotted.  TOPO should be formatted as a potential
  input to functions like 'contour' or 'image'.
}
\value{
  None; plot side effects only.
}

\author{Jake Anderson}

\examples{
ATM = CheckAtm.lin()
PlotAtm.lin(ATM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

