\name{K_arabidopsis}
\alias{K_atwell}
\alias{K_hapmap}
\alias{K_swedish}
\docType{data}
\title{
Marker-based relatedness matrices for 3 populations of Arabidopsis thaliana.
}
\description{
Marker-based relatedness matrices based on the SNP-data from Horton \emph{et al.} (2012).
Three matrices are provided:
(a) \code{K_atwell}, for the 199 accessions studied in Atwell \emph{et al.} (2010).
(b) \code{K_hapmap}, for a subset of 350 accessions taken from the Arabidopsis hapmap (Li \emph{et al.}, 2010).
(c) \code{K_swedish}, for 304 Swedish accessions.
All of these are part of the world-wide regmap of 1307 accessions, described in Horton \emph{et al.} (2012).
}
\usage{data(K_atwell); data(K_hapmap); data(K_swedish)}
\format{
Matrices whose row- and column names are the ecotype or seed-stock IDs of the accessions.
}
\details{
The matrices were computed using equation (2.2) in Astle and Balding (2009); see also
Goddard \emph{et al.} (2009). The \code{heritability}-package does not contain functions
to construct relatedness matrices from genotypic data, but such functions can be found in many other software packages.
For example, the \code{kin} function in the R-package \code{synbreed} (Wimmer \emph{et al.}, 2012), or
(outside R) GCTA (Yang \emph{et al.}, 2011), LDAK (Speed \emph{et al.}, 2012), Fast-LMM (Lippert, 2011) and GEMMA (Zhou and Stephens, 2012).
}
\references{
\itemize{
\item W. Astle and D.J. Balding (2009) Population Structure and Cryptic Relatedness in Genetic Association Studies.
\emph{Statistical Science}, Vol. 24, No. 4, 451-471.
\item Atwell, S., Y. S. Huang, B. J. Vilhjalmsson, G. Willems, M. Horton, \emph{et al.} (2010)
Genome-wide association study of 107 phenotypes in Arabidopsis thaliana inbred lines. \emph{Nature}
465: 627-631.
\item Goddard, M.E., Naomi R. Wray, Klara Verbyla and Peter M. Visscher (2009) Estimating Effects and Making Predictions from Genome-Wide Marker Data.
\emph{Statistical Science}, Vol. 24, No. 4, 517-529.
\item Horton, M. W., A. M. Hancock, Y. S. Huang, C. Toomajian, S. Atwell, \emph{et al.}
(2012) Genome-wide patterns of genetic variation in worldwide Arabidopsis thaliana accessions
from the RegMap panel. \emph{Nature Genetics} 44: 212-216.
\item Li, Y., Y. Huang, J. Bergelson, M. Nordborg, and J. O. Borevitz (2010) Association
mapping of local climate-sensitive quantitative trait loci in arabidopsis thaliana. \emph{PNAS} vol. 107, number 49.
\item Lippert, C., J. Listgarten, Y. Liu, C.M. Kadie, R.I. Davidson, \emph{et al.} (2011) FaST
linear mixed models for genome-wide association studies. \emph{Naure methods} 8: 833-835.
\item Speed, D., G. Hemani, M. R. Johnson, and D.J. Balding (2012) Improved heritability
estimation from genome-wide snps. \emph{the American journal of human genetics} 91: 1011-1021.
\item Wimmer, V., T. Albrecht, H-J. Auinger, and C-C. Schoen (2012) synbreed: a framework
for the analysis of genomic prediction data using r. \emph{bioinformatics} 28: 2086-2087.
\item Yang, J., S.H. Lee, M.E. Goddard, and P.M. Visscher (2011) GCTA: a tool for genomewide
complex trait analysis. \emph{the American journal of human genetics} 88: 76-82.
\item Zhou, X., and M. Stephens, (2012) Genome-wide efficient mixed-model analysis for association
studies. \emph{Nature genetics} 44: 821-824.
}
}
\seealso{
For phenotypic data for the population described in Atwell \emph{et al.} (2010), see \code{\link{LD}} and \code{\link{LDV}}.
For phenotypic data for the hapmap, see \code{\link{BT_LW_H}} and \code{\link{LA_H}}.
For phenotypic data for the Swedish regmap, see \code{\link{LA_S}}.
}
\examples{
data(K_atwell)
data(K_hapmap)
data(K_swedish)
}
\keyword{datasets}
