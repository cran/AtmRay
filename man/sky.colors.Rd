\name{sky.colors}
\alias{sky.colors}
\title{Sky Colormap}
\description{Returns a vector of colors one would see in a sunrise,
  ranging from orange to blue.
}
\usage{
sky.colors(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{length of output color vector}
}

\value{
  Vector of colors of length n.
}

\author{Jake Anderson}

\examples{
col = sky.colors(12)
M = matrix(1:12, 12, 1)
image(M, col = col)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

