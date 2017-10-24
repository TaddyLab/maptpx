\name{topicVar}
\alias{topicVar}
\alias{expit}
\alias{logit}
\title{
topic variance
}
\description{ Tools for looking at the variance of document-topic weights. }
\usage{
topicVar(counts, theta, omega) 
logit(prob)
expit(eta)
}
\arguments{
  \item{counts}{A matrix of multinomial response counts, as inputed to the \code{topics} or \code{predict.topics} functions.}
  \item{theta}{A fitted topic matrix, as ouput from the \code{topics} or \code{predict.topics} functions.}
  \item{omega}{A fitted document topic-weight matrix, as ouput from the \code{topics} or \code{predict.topics} functions.}
  \item{prob}{A probability vector (positive and sums to one) or a matrix with probability vector rows.}
  \item{eta}{A vector of the natural exponential family parameterization for a probability vector (with first category taken as null) or a matrix with each row the NEF parameters for a single observation.}
}
\details{ These function use the natural exponential family (NEF) parametrization of a probability vector \eqn{q_0 ... q_{K-1}} with the first element corresponding to a 'null' category; that is, with 
	  \eqn{NEF(q) = e_1 ... e_{K-1}} and setting \eqn{e_0 = 0}, the probabilities are
\deqn{q_k = \frac{exp[e_k]}{1 + \sum exp[e_j]}.} 
Refer to Taddy (2012) for details. }
\value{ \code{topicVar} returns an array with dimensions \eqn{(K-1,K-1,n)}, where \code{K=ncol(omega)=ncol(theta)} and \code{n = nrow(counts) = nrow(omega)}, filled with the posterior covariance matrix for the NEF parametrization of each row of \code{omega}.  Utility \code{logit} performs the NEF transformation and \code{expit} reverses it.  }
\references{
Taddy (2012), \emph{On Estimation and Selection for Topic Models}.
\url{http://arxiv.org/abs/1109.4518}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}

\seealso{
  topics, predict.topics
}
