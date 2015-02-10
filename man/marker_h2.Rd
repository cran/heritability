\name{marker_h2}
\alias{marker_h2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a marker-based estimate of heritability, given phenotypic observations
at individual plant or plot level. % (Kruijer \emph{et al.})
}
\description{
Given a genetic relatedness matrix and phenotypic observations at individual
plant or plot level, this function computes REML-estimates of the genetic and
residual variance and their standard errors, using the AI-algorithm (Gilmour \emph{et al.} 1995).
Based on this, heritability estimates and confidence intervals are given
(the estimator \eqn{h_r^2} in Kruijer \emph{et al.}).
}
\usage{
marker_h2(data.vector, geno.vector, covariates = NULL, K, alpha = 0.05,
          eps = 1e-06, max.iter = 100, fix.h2 = FALSE, h2 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.vector}{
A vector of phenotypic observations. Needs to be of type numeric. May contain
missing values.
}
  \item{geno.vector}{
A vector of genotype labels, either a factor or character. This vector should
correspond to \code{data.vector}, and hence needs to be of the same length.
}
  \item{covariates}{
A data-frame or matrix with optional covariates, the rows corresponding to
the phenotypic observations in \code{data.vector} and \code{geno.vector}.
May contain missing values. Factors are not allowed, and need to be encoded by
columns of type numeric or integer. The data-frame or matrix should not contain an intercept,
which is included by default.
}
  \item{K}{
A genetic relatedness or kinship matrix, typically marker-based.
Must have row- and column-names corresponding to the levels of \code{geno.vector}
}
  \item{alpha}{
Confidence level, for the 1-alpha confidence intervals.
}
  \item{eps}{
Numerical precision, used as convergence criterion in the AI-algorithm.
}
  \item{max.iter}{
Maximal number of iterations in the AI-algorithm.
}
  \item{fix.h2}{
Compute the log-likelihood and inverse AI-matrix for a fixed heritability value. Default is \code{FALSE}.
}
  \item{h2}{
When \code{fix.h2} is \code{TRUE}, the value of the heritability. Must be of type numeric, between 0 and 1.
}
}
\details{
\itemize{
%
\item Given phenotypic observations \eqn{Y_{ij}} for genotypes \eqn{i=1,...,n} and replicates
\eqn{j = 1,...,n_i}, the mixed model
\eqn{Y_{ij} = \mu + G_i + E_{ij}}
is assumed. The vector of additive genetic effects \eqn{(G_1,...,G_n)'} follows a
multivariate normal distribution with mean zero and covariance \eqn{\sigma_A^2 K},
where \eqn{\sigma_A^2} is the additive genetic variance, and \eqn{K} is a genetic relatedness matrix derived from a dense set of markers.
The errors \eqn{E_{ij}} are independent and normally distributed with variance \eqn{\sigma_E^2}.
Under certain assumptions (see Speed \emph{et al.} 2012) the marker- or chip-heritability \eqn{h^2 = \sigma_A^2 / (\sigma_A^2 + \sigma_E^2)}
equals the narrow-sense heritability.
%
\item It is assumed that the genetic relatedness matrix \eqn{K} is scaled such that \eqn{trace(P K P) = n - 1}, where
\eqn{P} is the projection matrix \eqn{I_n - 1_n 1_n' / n}, for the identity matrix \eqn{I_n} and \eqn{1_n} being a column vector of ones.
If this is not the case, \eqn{K} is automatically scaled prior to fitting the mixed model.
%
\item The model can optionally include a term \eqn{X_{ij} \beta}, where \eqn{X_{ij}}
is the row vector with observations on \eqn{k} extra covariates and the vector \eqn{\beta} contains their effects.
In this case the argument \code{covariates} should be the (N x k) matrix or
data-frame with rows \eqn{X_{ij}} (N being the total number of observations).
Observations where either \eqn{Y_{ij}} or any of the covariates is missing are discarded.
%
\item Confidence intervals for heritability are constructed using the delta-method and the inverse AI-matrix.
The delta-method can be applied either directly to the
function \eqn{(\sigma_A^2,\sigma_E^2) -> \sigma_A^2 / (\sigma_A^2 + \sigma_E^2)} or to the function
\eqn{(\sigma_A^2,\sigma_E^2) -> log(\sigma_A^2 / \sigma_E^2)}. In the latter case,
a confidence interval for \eqn{log(\sigma_A^2 / \sigma_E^2)} is obtained, which is back-transformed to a confidence interval for heritability.
This approach (proposed in Kruijer \emph{et al.}) has the advantage that intervals are always contained in the unit interval.
%
\item The AI-algorithm is run for \code{max.iter} iterations. If by then there is no convergence a warning is printed and the current estimates are returned.
}
}
\value{
A list with the following components:
\itemize{
  \item{va}{:  REML-estimate of the (additive) genetic variance.}
  \item{ve}{:  REML-estimate of the residual variance.}
  \item{h2}{:  Plug-in estimate of heritability: \eqn{va / (va + ve)}.}
  \item{conf.int1}{:  1-alpha confidence interval for heritability.}
  \item{conf.int2}{:  1-alpha confidence interval for heritability, obtained by
  application of the delta method on a logarithmic scale.}
  \item{inv.ai}{:  The inverse of the average information (AI) matrix.}
  \item{loglik}{:  The log-likelihood.}
}
}
\references{
\itemize{
\item Gilmour \emph{et al.} Gilmour, A.R., R. Thompson and B.R. Cullis (1995) Average Information {REML}: An Efficient Algorithm for Variance Parameter
	Estimation in Linear Mixed Models. \emph{Biometrics}, volume 51, number 4, 1440-1450.
\item Kruijer, W. \emph{et al.} (2015) Marker-based estimation of heritability in immortal
populations. \emph{Genetics}, Vol. 199(2), p. 1-20.
\item Speed, D., G. Hemani, M. R. Johnson, and D.J. Balding (2012) Improved heritability
estimation from genome-wide snps. \emph{the American journal of human genetics} 91: 1011-1021.
}
}
\author{
Willem Kruijer.
}
%\note{
%Works well for small data-sets (i.e. at most 1000 observations).
%}
\seealso{
For marker-based estimation of heritability using genotypic means, see
\code{\link{marker_h2_means}}.
}
\examples{
data(LD)
data(K_atwell)
# Heritability estimation for all observations:
#out <- marker_h2(data.vector=LD$LD,geno.vector=LD$genotype,
#                 covariates=LD[,4:8],K=K_atwell)
# Heritability estimation for a randomly chosen subset of 20 accessions:
set.seed(123)
sub.set <- which(LD$genotype \%in\% sample(levels(LD$genotype),20))
out <- marker_h2(data.vector=LD$LD[sub.set],geno.vector=LD$genotype[sub.set],
                 covariates=LD[sub.set,4:8],K=K_atwell)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ heritability }
\keyword{ narrow-sense }
\keyword{ marker-based }