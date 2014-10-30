\name{leafArea}
\alias{leafArea}
\alias{LA_H}
\alias{LA_S}
\docType{data}
\title{
\emph{Arabidopsis} leaf area data for the hapmap and Swedish regmap population.
}
\description{
\emph{Arabidopsis} leaf area data for the hapmap and Swedish regmap population.
}
\usage{data(LA_H); data(LA_S)}

\format{
  Data frame with leaf area observations:
  \describe{
    \item{\code{genotype}}{a factor, the levels being the accession identifiers}
    \item{\code{LA13_H}}{Leaf area 13 days after sowing, in numbers of pixels (hapmap)}
    \item{\code{LA13_S}}{Leaf area 13 days after sowing, in numbers of pixels (Swedish regmap)}
    \item{\code{replicate}}{The replicate (or block) each plant is contained in (factor with levels \code{1} to \code{4})}
    \item{\code{rep1}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 1 and 0 otherwise}
    \item{\code{rep2}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 2 and 0 otherwise}
    \item{\code{rep3}}{numeric encoding of the factor replicate: equals 1 if the plant is in replicate 3 and 0 otherwise}
    \item{\code{x}}{The within image x-coordinate of the plant. A factor with levels \code{1} \code{2} \code{3}}
    \item{\code{y}}{The within image y-coordinate of the plant. A factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{x1}}{numeric encoding of the factor \code{x}: equals 1 if the plant is in position 1 and 0 otherwise}
    \item{\code{x2}}{numeric encoding of the factor \code{x}: equals 1 if the plant is in position 2 and 0 otherwise}
    \item{\code{y1}}{numeric encoding of the factor \code{y}: equals 1 if the plant is in position 1 and 0 otherwise}
    \item{\code{y2}}{numeric encoding of the factor \code{y}: equals 1 if the plant is in position 2 and 0 otherwise}
    \item{\code{y3}}{numeric encoding of the factor \code{y}: equals 1 if the plant is in position 3 and 0 otherwise}
}
}
%\details{
%add.
%}
\author{
Willem Kruijer <willlem.kruijer@wur.nl>;
experiments conducted by Padraic Flood <flood@mpipz.mpg.de>
}
%\source{
%Kruijer, W. \emph{et al.}.
%}
\references{
\itemize{
\item Kruijer, W. \emph{et al.} (submitted) Marker-based estimation of heritability in immortal
populations.
}
}
\seealso{
For the corresponding genetic relatedness matrices, see \code{\link{K_hapmap}} and \code{\link{K_swedish}}.
}
\examples{
data(LA_H); data(LA_S)
}
\keyword{datasets}
