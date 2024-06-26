\name{BF_dist}
\alias{BF_dist}
\title{Generate bootstrap samples for the bootstrapping Fieller's confidence interval of the penalized inverse-variance weighted (pIVW) method}
\usage{
BF_dist(
  object,
  beta_hat = 0,
  tau2 = 0,
  lambda = 1,
  n_boot = 1000,
  seed_boot = 1
)
}
\arguments{
\item{object}{An data frame containing Bx, By, Bxse and Byse.}

\item{beta_hat}{The causal effect estimate.}

\item{tau2}{The estimated variance of the horizontal pleiotropy.}

\item{lambda}{The penalty parameter in the pIVW estimator. By default, \code{lambda=1}.}

\item{n_boot}{The number of the bootstrap samples. By default, \code{n_boot=1000}.}

\item{seed_boot}{The seed for random sampling in the bootstrap method. By default, \code{seed_boot=1}.}
}
\value{
\item{z_b}{A vector containing the bootstrap samples for the bootstrapping Fieller's confidence interval.}
}
\description{
Internal function of the penalized inverse-variance weighted (pIVW) method, which generates bootstrap samples for the bootstrapping Fieller's confidence interval.
}
\keyword{internal}
