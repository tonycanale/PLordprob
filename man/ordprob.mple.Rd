\name{ordprob.mple}
\alias{ordprob.mple}
\title{
Maximum Pairwise Likelihood Estimation
}
\description{
Estimate the parameters of a multivariate ordered probit model maximizing the pairwise log-likelihood function.
}
\usage{
ordprob.mple(y, x = NULL, K, start.par = list(type = "default"), 
       same.means = FALSE, eval.max = 1000, iter.max = 600, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
matrix with n rows and q columns containing the categorical responses. Each line is a vector of size q representing the responses for a single statistical unit.
}
  \item{x}{
additional n times p matrix of subject specific covariates. 
}
  \item{K}{
Number of levels of the categorical responses.
}
  \item{start.par}{
list containing parameters for the maximization algorithm:
\code{type} is character string. If "default" the initialization is as default. Otherwise each value should be passed, \code{cor} is a vector with the initial values for the polycoric correlations, \code{beta} is a vector with the initial values for the regression coefficients, \code{xi} is a vector with the initial values for the column specific means, and \code{thres} is a vector with the initial values for the thresholds.
}
  \item{same.means}{
logical. If code{TRUE} all the q variable are assumed to have the same mean.
}
  \item{eval.max}{
see \code{help(nlminb)}
}
  \item{iter.max}{
see \code{help(nlminb)}
}
  \item{\dots}{
additional arguments to be passed.
}
}
\details{
The code is implemented in R software with call to C functions for the most demanding operations. To evaluate the Gaussian integrals, the package uses the Fortran 77 subroutine SADMVN. The default choice of initialization is the first threshold equal to zero and the remaining thresholds equally spaced with distance one. As for the covariance components, we consider as  starting values  the sample covariances of the observed categorical variables treated as continuous. Optimization of the pairwise log-likelihood function is performed via quasi-Newton box-constrained optimization algorithm, as implemented in \code{nlminb}.}
\value{
A list with components:
\item{par}{The best set of parameters found.}
\item{objective}{The value of the negative pairwise likelihood  corresponding to \code{par}}
\item{convergence}{An integer code. 0 indicates successful convergence.}
\item{message}{A character string giving any additional information returned by the optimizer, or \code{NULL}. For details, see \code{nlminb} documentation.}
\item{iterations}{Number of iterations performed.}
\item{evaluations}{Number of objective function and gradient function evaluations.}
\item{thresh}{The set of thresholds partitioning the latent sample space.}
\item{xi}{Vector of the item means.}
\item{cor}{Estimated polychoric correlation matrix.}
}
\references{
Cox D. R. , Reid N. (2004) A note on pseudolikelihood constructed from marginal densities. \emph{Biometrika}, \bold{91}, 729--737.

Genz A. (1992) Numerical computation of multivariate normal probabilities. \emph{Journal of computational and graphical statistics}, \bold{2}, 141--149.

Kenne Pagui, E. C. and Canale, A. (2014) Pairwise likelihood inference for multivariate categorical responses, Technical Report, Department of Statistics, University of Padua. 

Lindsay B. (1988) Composite likelihood methods. \emph{Comtemporary Mathematics}, \bold{80}, 221--240.
}
\examples{
library(mnormt)

\dontrun{# constrained polychoric correlation matrix
sigma <- matrix(0.6,5,5)
diag(sigma) <- 1
set.seed(123)
y <- rmnorm(50, varcov=0.5*sigma)
y <- floor(y)+3
fit <- ordprob.mple(y, K=4, same.mean=TRUE)
#polychoric correlations
fit$cor

# unconstrained polychoric correlation matrix
set.seed(1)
S <- matrix(rWishart(1, df=5, Sigma=diag(1,5)), 5, 5)
prec <- diag(1/sqrt(diag(S)))
cor <- prec\%*\%S\%*\%prec
sigma <- cor
set.seed(123)
y=rmnorm(70, c(1,1,0,0,0), varcov=sigma)
y=floor(y)+4
fit <- ordprob.mple(y, K=7, same.mean=FALSE)
#polychoric correlations
fit$cor
}
}
