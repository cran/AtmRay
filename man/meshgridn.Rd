\name{meshgridn}
\alias{meshgridn}
\title{Multi-Dimensional Meshgrid}
\description{Inputs a list of vectors, and returns a list of vectors
  such that every possible combination of input vector values occurs
  once in the output.
}
\usage{
meshgridn(L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{list of vectors}
}
\value{
  List of vectors of equal length such that each combination of elements
  drawn from the input vectors occurs exactly once in the output list.
  }
\author{Jake Anderson}

\examples{
meshgridn(list(1:5, 10:12, 100:104))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

