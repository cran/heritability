\name{BT_LW_H}
\alias{BT_LW_H}
\docType{data}
\title{
Bolting time and leaf width for the Arabidopsis hapmap population.
}
\description{
Bolting time and leaf width for the Arabidopsis hapmap population
}
\usage{data(BT_LW_H)}
\format{
  A data frame with phenotypic observations on bolting time and leaf width:
  \describe{
    \item{\code{genotype}}{a factor, the levels being the accession or ecotype identifiers}
    \item{\code{BT}}{Bolting time, in number of days}
    \item{\code{LW}}{Leaf width}
    \item{\code{replicate}}{The replicate (or block) each plant is contained in (factor with levels \code{1} to \code{3})}
    \item{\code{rep1}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 1 and 0 otherwise}
\item{\code{rep2}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 2 and 0 otherwise}
  }
}
%\details{
%Randomized complete block design. Missing values.
%}
%\source{
%Kruijer \emph{et al.} (2014).
%}
\author{
Willem Kruijer <willlem.kruijer@wur.nl>;
experiments conducted by Rik Kooke <rik.kooke@gmail.com>
}
\references{
\itemize{
\item Kruijer, W. \emph{et al.} (2015) Marker-based estimation of heritability in immortal
populations. \emph{Genetics}, Vol. 199(2), p. 1-20.
}
}
\seealso{
For the corresponding genetic relatedness matrix, see \code{\link{K_hapmap}}.
}
\examples{
data(BT_LW_H)
}
\keyword{datasets}
