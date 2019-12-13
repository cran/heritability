\name{means_floweringTime}
\alias{means_floweringTime}
\alias{means_LD}
\alias{means_LDV}
\docType{data}
\title{
Flowering time from Atwell \emph{et al.} (2010): accession means.
}
\description{
Accession means for the flowering time data contained in \code{\link{LD}} and \code{\link{LDV}}.
}
\usage{data(means_LD); data(means_LDV)}
\format{
  Data-frames with flowering time means:
  \describe{
    \item{\code{genotype}}{a factor, the levels being the accession or ecotype identifiers}
    \item{\code{LD}}{Flowering time under long days, in number of days}
    \item{\code{LDV}}{Flowering time under long days and vernalization, in number of days}
  }
}
\details{
Following Kruijer \emph{et al.} (appendix A)
these means were defined as the least-squares estimate for the factor accession, in a linear model
containing both accession and replicate effects.
Consequently there are differences compared to Atwell \emph{et al.} (2010), where just the
arithmetic averages are considered.
}
%\source{
%Atwell \emph{et al.} (2010).
%}
\references{
\itemize{
\item Atwell, S., Y. S. Huang, B. J. Vilhjalmsson, G. Willems, M. Horton, \emph{et al.} (2010)
Genome-wide association study of 107 phenotypes in Arabidopsis thaliana inbred lines. Nature
465: 627-631.
\item Kruijer, W. \emph{et al.} (2015) Marker-based estimation of heritability in immortal
populations. \emph{Genetics}, Vol. 199(2), p. 1-20.
}
}
\seealso{
Together with the covariance matrices contained in \code{\link{R_matrix_LD}} and
\code{\link{R_matrix_LDV}}, the means contained in \code{means_LD} and \code{means_LDV}
can be used to estimate heritability, using the function \code{\link{marker_h2_means}}.
For the corresponding genetic relatedness matrix, see \code{\link{K_atwell}}.
For the individual plant data, see \code{\link{floweringTime}}.
}
\examples{
data(means_LD)
data(means_LDV)
}
\keyword{datasets}
