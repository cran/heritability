\name{repeatability}
\alias{repeatability}
\title{ANOVA-based estimates of repeatability}
\description{
Given a population where each genotype is phenotyped for a number of genetically
identical replicates (either individual plants or plots in a field trial), the repeatability
or intra-class correlation (Kempthorne, 1957) can be estimated by
\eqn{V_g / (V_g + V_e)}, where \eqn{V_g = (MS(G) - MS(E)) / r} and \eqn{V_e = MS(E)}. In these expressions,
\eqn{r} is the number of replicates per genotype, and \eqn{MS(G)} and \eqn{MS(E)} are
the mean sums of squares for genotype and residual error obtained from analysis
of variance. In case \eqn{MS(G) < MS(E)}, \eqn{V_g} is set to zero.
See Singh \emph{et al.} (1993) or Lynch and Walsh (1998), p.563.
When the genotypes have differing numbers of replicates, \eqn{r} is replaced by
\eqn{\bar r = (n-1)^{-1} (R_1 - R_2 / R_1)}, where \eqn{R_1 = \sum r_i} and \eqn{R_2 = \sum r_i^2}.
Under the assumption that all differences between genotypes are genetic,
repeatability equals broad-sense heritability; otherwise it only provides an upper-bound for broad-sense heritability.
}
\usage{
repeatability(data.vector, geno.vector, line.repeatability = FALSE,
              covariates.frame = data.frame())
}
\arguments{
\item{data.vector}{A vector of phenotypic observations. Needs to be of type numeric. May contain
missing values.}
\item{geno.vector}{A vector of genotype labels, either a factor or character. This vector should
correspond to \code{data.vector}, and hence needs to be of the same length.}
\item{line.repeatability}{If \code{TRUE}, the line-repeatability or
line-heritability \eqn{\sigma_G^2 / (\sigma_G^2 + \sigma_E^2 / r)} is estimated,
otherwise (the default) the repeatability at plot- or plant level, which is \eqn{\sigma_G^2 / (\sigma_G^2 + \sigma_E^2)}.
}
\item{covariates.frame}{A data-frame with additional covariates, the rows corresponding to
\code{geno.vector} and the phenotypic observations in \code{data.vector}.
May contain missing values. Each column can be numeric or a factors.}
}
\value{
A list with the following components:
\itemize{
  \item{repeatability}{:  the estimated repeatability.}
  \item{gen.variance}{:  the estimated genetic variance.}
  \item{res.variance}{:  the estimated residual variance.}
  \item{line.repeatability}{:  whether repeatability was estimated at the
  individual plant or plot level (the default), or at the level of
  genotypic means (in the latter case, \code{line.repeatability=TRUE}) }
  \item{average.number.of.replicates}{:  The average number of replicates. See the description above.}
  \item{conf.int}{:  Confidence interval for repeatability. See Singh \emph{et al.} (1993) or Lynch and Walsh (1998)}
}
}
\references{
\itemize{
\item Kempthorne, Oscar (1957) An introduction to genetic statistics. \emph{New York: John Wiley and Sons}
\item Kruijer, W. \emph{et al.} (submitted) Marker-based estimation of heritability in immortal
populations.
\item Lynch, M., and B. Walsh (1998) Genetics and Analysis of Quantitative Traits. Sinauer As-
sociates, 1st edition.
\item Singh, M., S. Ceccarelli, and J. Hamblin (1993) Estimation of heritability from varietal
trials data. \emph{Theoretical and Applied Genetics} 86: 437-441.
}
}
\examples{
repeatability(data.vector=rep(rnorm(26),each=5) + rnorm(5*26),
              geno.vector=rep(letters,each=5))
}
\author{
Willem Kruijer \email{willem.kruijer@wur.nl}
}

