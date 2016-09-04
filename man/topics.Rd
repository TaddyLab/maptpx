\name{topics}
\alias{topics}
\title{
Estimation for Topic Models
}
\description{ MAP estimation of Topic models  }
\usage{
topics(counts, K, shape=NULL, initopics=NULL, 
  tol=0.1, bf=FALSE, kill=2, ord=TRUE, verb=1,
  admix=TRUE, nbundles=1, use_squarem=FALSE,
  init.adapt=FALSE, light=1, method_admix=1, ...)
}
\arguments{
  \item{counts}{
    A matrix of multinomial response counts in \code{ncol(counts)} phrases/categories 
    for \code{nrow(counts)} documents/observations. Can be either a simple \code{matrix} or a \code{simple_triplet_matrix}.}
  \item{K}{The number of latent topics.  If \code{length(K)>1}, 
  \code{topics} will find the Bayes factor (vs a null single topic model) for each element and return parameter
  estimates for the highest probability K.}
  \item{shape}{ Optional argument to specify the Dirichlet prior concentration parameter as \code{shape} for topic-phrase probabilities.  Defaults to \code{1/(K*ncol(counts))}.  For fixed single \code{K}, this can also be a \code{ncol(counts)} 
  		by \code{K} matrix of unique shapes for each topic element.}
  \item{initopics}{ Optional start-location for \eqn{[\theta_1 ... \theta_K]}, the topic-phrase probabilities.  
  		    Dimensions must accord with the smallest element of \code{K}.  
  		    If \code{NULL}, the initial estimates are built by incrementally adding topics. }
  \item{tol}{ Convergence tolerance: optimization stops, conditional on some extra checks, when the \emph{relative} posterior increase over a full paramater set update is less than \code{tol}.}
  \item{bf}{An indicator for whether or not to calculate the Bayes factor for univariate \code{K}.  
  If \code{length(K)>1}, this is ignored and Bayes factors are always calculated. }
  \item{kill}{For choosing from multiple \code{K} numbers of topics (evaluated in increasing order),
   the search will stop after \code{kill} consecutive drops in the corresponding Bayes factor.  Specify \code{kill=0}
   if you want Bayes factors for all elements of \code{K}. }
 \item{ord}{If \code{TRUE}, the returned topics (columns of \code{theta}) will be ordered by decreasing usage (i.e., 
 by decreasing \code{colSums(omega)}).}
  \item{verb}{A switch for controlling printed output.  \code{verb > 0} will print something, with the level of detail
  increasing with \code{verb}.}
  \item{admix}{if TRUE, admixture model used, else a mixture model is used. Defaults to TRUE}
  \item{nbundles}{apply the active set method update to the topic proportions \code{omega} after these many steps if \code{light=1}. Defaults to 1, which essentially means the update is used at each step.}
  \item{use_squarem}{if TRUE, uses the SQUAREM implementation due to Ravi Varadhan. Else uses the usual EM + Quasi Newton Acceleration update.}
  \item{init.adapt}{if TRUE, adpatively initializes the topic distributions \code{theta}. Deafults to FALSE.}
  \item{light}{May tale values 0,1,2. If \code{light=2}, we select a set of genes after a burn-in period for each cluster and iterate over these selected genes. If \code{light=1}, the active set method update occurs after every \code{nbundles} steps. If \code{light=0}, no active set method is used for updating the topic proportions \code{omega}. Defaults to 1.}
  \item{method_admix}{May take values 1 or 2. If method_admix=1, R code is used for EM update. If \code{emthod_admix=2}, C code is used for EM update. Usually C code is preferredfor small datasets, but R code is performs more stablly for large data and is also faster. Defaults to 1.}
  
  \item{...}{Additional arguments to the undocumented internal \code{tpx*} functions. }
 }
\details{A latent topic model represents each i'th document's term-count vector \eqn{X_i} 
(with \eqn{\sum_{j} x_{ij} = m_i} total phrase count)
as having been drawn from a mixture of \code{K} multinomials, each parameterized by topic-phrase
probabilities \eqn{\theta_i}, such that 
\deqn{X_i \sim MN(m_i, \omega_1 \theta_1 + ... + \omega_K\theta_K).}
We assign a K-dimensional Dirichlet(1/K) prior to each document's topic weights 
\eqn{[\omega_{i1}...\omega_{iK}]}, and the prior on each \eqn{\theta_k} is Dirichlet with concentration \eqn{\alpha}.
The \code{topics} function uses quasi-newton accelerated EM, augmented with sequential quadratic programming 
for conditional \eqn{\Omega | \Theta} updates, to obtain MAP estimates for the topic model parameters. 
We also provide Bayes factor estimation, from marginal likelihood
calculations based on a Laplace approximation around the converged MAP parameter estimates.  If  input \code{length(K)>1}, these
Bayes factors are used for model selection. Full details are in Taddy (2011). }
\note{ Estimates are actually functions of the MAP (K-1 or p-1)-dimensional 
logit transformed natural exponential family parameters.  }
\value{
 An \code{topics} object list with entries
\item{K}{The number of latent topics estimated.  If input \code{length(K)>1}, 
on output this is a single value corresponding to the model with the highest Bayes factor. }
\item{theta}{The \code{ncol{counts}} by \code{K} matrix of estimated topic-phrase probabilities.}
\item{omega}{The \code{nrow{counts}} by \code{K} matrix of estimated document-topic weights. }
\item{BF}{The log Bayes factor for each number of topics in the input \code{K}, against a null single topic model.}
\item{D}{Residual dispersion: for each element of \code{K}, estimated dispersion parameter
 (which should be near one for the multinomial), degrees of freedom, and p-value for a test of whether the true dispersion is \eqn{>1}.}
\item{X}{The input count matrix, in \code{dgTMatrix} format.}
}
\references{
Taddy (2012), \emph{On Estimation and Selection for Topic Models}.
\url{http://arxiv.org/abs/1109.4518}
}
\author{
  Matt Taddy \email{taddy@chicagobooth.edu}
}
\seealso{
  plot.topics, summary.topics, predict.topics, wsjibm, congress109, we8there
}
\examples{
## see also http://faculty.chicagobooth.edu/matt.taddy/teaching/we8there.R

## Simulation Parameters
K <- 10
n <- 100
p <- 100
omega <- t(rdir(n, rep(1/K,K)))
theta <- rdir(K, rep(1/p,p))

## Simulated counts
Q <- omega\%*\%t(theta)
counts <- matrix(ncol=p, nrow=n)
totals <- rpois(n, 100)
for(i in 1:n){ counts[i,] <- rmultinom(1, size=totals[i], prob=Q[i,]) }

## Bayes Factor model selection (should choose K or nearby)
summary(simselect <- topics(counts, K=K+c(-5:5)), nwrd=0)

## MAP fit for given K
summary( simfit <- topics(counts,  K=K, verb=2), n=0 )

## Adjust for label switching and plot the fit (color by topic)
toplab <- rep(0,K)
for(k in 1:K){ toplab[k] <- which.min(colSums(abs(simfit$theta-theta[,k]))) }
par(mfrow=c(1,2))
tpxcols <- matrix(rainbow(K), ncol=ncol(theta), byrow=TRUE)
plot(theta,simfit$theta[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
plot(omega,simfit$omega[,toplab], ylab="fitted values", pch=21, bg=tpxcols)
title("True vs Fitted Values (color by topic)", outer=TRUE, line=-2)

## The S3 method plot functions
par(mfrow=c(1,2))
plot(simfit, lgd.K=2)
plot(simfit, type="resid")

}
