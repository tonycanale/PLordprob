\name{PLordprob-package}
\alias{PLordprob-package}
\alias{PLordprob}
\docType{package}
\title{
Multivariate Ordered Probit Model via Pairwise Likelihood
}
\description{
Fit multivariate ordered probit model, i.e. the extension of
the scalar ordered probit model where the observed variables have
dimension greater than one. Estimation of the parameters is done via
maximization of the pairwise likelihood, a special case of the
composite likehood obtained as product of bivariate marginal
distributions. The package uses the Fortran 77 subroutine SADMVN 
by Alan Genz, with minor adaptations made by Adelchi Azzalini 
in his "mvnormt" package for evaluating the two-dimensional Gaussian 
integrals involved in the pairwise log-likelihood. Optimization 
of the latter objective function  is performed via quasi-Newton 
box-constrained optimization algorithm, as implemented in nlminb.
}
\details{
\tabular{ll}{
Package: \tab PLordprob\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-29\cr
License: GPL-2
}
}
\author{
Euloge Clovis Kenne Pagui [aut,cre], Antonio Canale [aut], Alan Genz [ctb],  Adelchi Azzalini [ctb] \cr
Maintainer: Euloge Clovis Kenne Pagui <kenne@stat.unipd.it>
}
\references{
Cox D. R. , Reid N. (2004) A note on pseudolikelihood constructed from marginal densities. \emph{Biometrika}, \bold{91}, 729--737.

Genz A. (1992) Numerical computation of multivariate normal probabilities. \emph{Journal of computational and graphical statistics}, \bold{2}, 141--149.

Kenne Pagui, E. C. and Canale, A. (2014) Pairwise likelihood inference for multivariate categorical responses, Technical Report, Department of Statistics, University of Padua. 

Lindsay B. (1988) Composite likelihood methods. \emph{Comtemporary Mathematics}, \bold{80}, 221--240.
}

