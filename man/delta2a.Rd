\name{delta2a}
\alias{delta2a}
\alias{a2delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tresholds Mapping
}
\description{
Functions to transform (and antitrasforms) the thresholds into an unconstrained space.
}
\usage{
delta2a(delta)
a2delta(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{ a set of increaing thresholds}
  \item{a}{ a vector (of size \code{length(delta)-1}) of threshold log-differences}
}
\details{
delta_k = log(a_k - a_{k+1})
}
\value{
A vector of tresholds or of tresholds log-differences
}
\author{
Euloge Clovis Kenne Pagui and Antonio Canale
}
\examples{
thresh <- c(0,2,4)
logdiffs <- a2delta(thresh)
}
