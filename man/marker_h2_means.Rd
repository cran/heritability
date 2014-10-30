\name{marker_h2_means}
\alias{marker_h2_means}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a marker-based estimate of heritability, given genotypic means.% (Kruijer \emph{et al.})
}
\description{
Given a genetic relatedness matrix and genotypic means,
this function computes REML-estimates of the genetic and
residual variance and their standard errors, using the AI-algorithm (Gilmour \emph{et al.} 1995).
Based on this, heritability estimates and confidence intervals are given
(the estimator \eqn{h_m^2} in Kruijer \emph{et al.}).
}
\usage{
marker_h2_means(data.vector, geno.vector, K, Dm=NULL, alpha = 0.05, eps = 1e-06,
       max.iter = 100, fix.h2 = FALSE, h2 = 0.5, grid.size=99)
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
  \item{K}{
A genetic relatedness or kinship matrix, typically marker-based.
Must have row- and column-names corresponding to the levels of \code{geno.vector}
}
\item{Dm}{
See details.
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
  \item{grid.size}{If the AI-algorithm has not converged after \code{max.iter} iterations,
the likelihood is computed on the grid of heritability values 1/(\code{grid.size}+1),...,\code{grid.size}/(\code{grid.size}+1); see details.
}
}
\details{
\itemize{
%
\item Given phenotypic observations \eqn{Y_{i}} for genotypes \eqn{i=1,...,n}, the mixed model
\eqn{Y_{i} = \mu + G_i + E_{i}}
is assumed. Typically, the \eqn{Y_{i}} are genotypic means or BLUEs obtained from fitting a linear (mixed) model to
the raw data, containing several plants or plots for each genotype.
The vector of additive genetic effects \eqn{(G_1,...,G_n)'} follows a
multivariate normal distribution with mean zero and covariance \eqn{\sigma_A^2 K},
where \eqn{\sigma_A^2} is the additive genetic variance, and \eqn{K} is a genetic relatedness matrix derived from a dense set of markers.
The vector of errors \eqn{(E_1,...,E_n)'} follows a
multivariate normal distribution with mean zero and covariance \eqn{\sigma_E^2 R},
where \eqn{R} is the covariance of the means obtained from the initial analysis.
In case of a completely randomized design with \eqn{r_i} replicates for genotypes
\eqn{i=1,...,n}, \eqn{R} is diagonal with elements \eqn{1 / r_i}.
Under certain assumptions (see Speed \emph{et al.} 2012) the marker- or chip-heritability \eqn{h^2 = \sigma_A^2 / (\sigma_A^2 + \sigma_E^2)}
equals the narrow-sense heritability.
%
\item As in the \code{\link{marker_h2}} function, it is assumed that the genetic relatedness matrix \eqn{K} is scaled such that \eqn{trace(P K P) = n - 1}, where
\eqn{P} is the projection matrix \eqn{I_n - 1_n 1_n' / n}, for the identity matrix \eqn{I_n} and \eqn{1_n} being a column vector of ones.
If this is not the case, \eqn{K} is automatically scaled prior to fitting the mixed model.
%
\item No covariates can be included, as it is assumed that these are available at plant- or plot level, and accounted for in the genotypic means.
%
\item The resulting heritability estatimes are less accurate than those obtained from individual plant or plot data, and the likelihood can be monotone
in \eqn{h^2 = \sigma_A^2 / (\sigma_A^2 + \sigma_E^2)}. If the AI-algorithm has not converged after \code{max.iter} iterations,
the likelihood is computed on the grid of heritability values 1/(\code{grid.size}+1),...,\code{grid.size}/(\code{grid.size}+1)
%
\item As in the \code{\link{marker_h2}} function, confidence intervals for heritability are constructed using the delta-method and the inverse AI-matrix.
The delta-method can be applied either directly to the
function \eqn{(\sigma_A^2,\sigma_E^2) -> \sigma_A^2 / (\sigma_A^2 + \sigma_E^2)} or to the function
\eqn{(\sigma_A^2,\sigma_E^2) -> log(\sigma_A^2 / \sigma_E^2)}. In the latter case,
a confidence interval for \eqn{log(\sigma_A^2 / \sigma_E^2)} is obtained, which is back-transformed to a confidence interval for heritability.
This approach (proposed in Kruijer \emph{et al.}) has the advantage that intervals are always contained in the unit interval.
%
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
  \item{loglik.vector}{:  Empty numeric vector if the AI-algorthm converged within \eqn{max.iter} iterations. Otherwise it contains the log-likelihood on a grid.}
}
}
%   \item{inv.ai.varcomp}{To do.}
\references{
\itemize{
\item Gilmour \emph{et al.} Gilmour, A.R., R. Thompson and B.R. Cullis (1995) Average Information {REML}: An Efficient Algorithm for Variance Parameter
	Estimation in Linear Mixed Models. \emph{Biometrics}, volume 51, number 4, 1440-1450.
\item Kruijer, W. \emph{et al.} (submitted) Marker-based estimation of heritability in immortal
populations.
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
For marker-based estimation of heritability using individual plant or plot data, see
\code{\link{marker_h2}}.
}
\examples{
data(means_LDV)
data(R_matrix_LDV)
data(K_atwell)
out <- marker_h2_means(data.vector=means_LDV$LDV,geno.vector=means_LDV$genotype,
                       K=K_atwell,Dm=R_matrix_LDV)
# Takes about a minute:
#data(means_LD)
#data(R_matrix_LD)
#out <- marker_h2_means(data.vector=means_LD$LD,geno.vector=means_LD$genotype,
#                       K=K_atwell,Dm=R_matrix_LD)
# The likelihood is monotone increasing:
#plot(x=(1:99)/100,y=out$loglik.vector,type="l",ylab="log-likelihood",lwd=2,
#     main='',xlab='h2',cex.lab=2,cex.axis=2.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ heritability }
\keyword{ narrow-sense }
\keyword{ marker-based }