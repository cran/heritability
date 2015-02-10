\name{floweringTime}
\alias{floweringTime}
\alias{LD}
\alias{LDV}
\docType{data}
\title{
Flowering time data taken from Atwell \emph{et al.} (2010).
}
\description{
Two data-frames containing individual plant data on flowering time under different conditions:
\code{LDV} (Flowering time under long days and vernalization) and \code{LD} (Flowering time under long days, without vernalization).
}
\usage{data(LD); data(LDV)}
\format{
  Data-frames with flowering time observations, genotype and design information:
  \describe{
    \item{\code{genotype}}{a factor, the levels being the accession or ecotype identifiers}
    \item{\code{LD}}{Flowering time under long days, in number of days}
    \item{\code{LDV}}{Flowering time under long days and vernalization, in number of days}
    \item{\code{replicate}}{The replicate (or block) each plant is contained in (factor with levels \code{1} to \code{6})}
    \item{\code{rep1}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 1 and 0 otherwise}
\item{\code{rep2}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 2 and 0 otherwise}
\item{\code{rep3}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 3 and 0 otherwise}
\item{\code{rep4}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 4 and 0 otherwise}
\item{\code{rep5}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 5 and 0 otherwise}
  }
}
\details{
All plants that had not flowered by the end of the experiment were given a phenotypic value of 200.
Only accessions for which SNP-data are available are included here: 167 accessions in case of LD and 168 accessions in case of LDV.
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
For the corresponding genetic relatedness matrix, see \code{\link{K_atwell}}.
}
\examples{
data(LD); data(LDV)
}
\keyword{datasets}
