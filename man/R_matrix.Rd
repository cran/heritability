\name{R_matrix}
\alias{R_matrix}
\alias{R_matrix_LDV}
\alias{R_matrix_LD}
\docType{data}
\title{
Covariance matrix of the accession means for flowering time.%Marker-based kinship matrix based on the Atwell \emph{et al.} (2010) data.
}
\description{
Covariance matrices of the accession means for flowering time contained in \code{\link{means_LD}}
and \code{\link{means_LDV}}, derived from the Atwell \emph{et al.} (2010) data.
}
\usage{data(R_matrix_LDV);data(R_matrix_LD)}
\format{
Matrix whose row- and column names are the ecotype-IDs of the accessions contained in
\code{\link{LD}} and \code{\link{LDV}}.
}
\details{
The matrix was computed as in Kruijer \emph{et al.}, Appendix A.
}
%\source{
%Atwell \emph{et al.} (2010).
%}
\references{
\itemize{
\item Atwell, S., Y. S. Huang, B. J. Vilhjalmsson, G. Willems, M. Horton, \emph{et al.} (2010)
Genome-wide association study of 107 phenotypes in Arabidopsis thaliana inbred lines. Nature
465: 627-631.
\item Kruijer, W. \emph{et al.} (submitted) Marker-based estimation of heritability in immortal
populations.
}
}
\seealso{
Together with the corresponding means contained in \code{\link{means_LD}} and \code{\link{means_LDV}},
these matrices can be used to estimate heritability, using the function \code{\link{marker_h2_means}}.
}
%For the corresponding accession means, see \code{\link{means_LD}} and \code{\link{means_LDV}}.
\examples{
data(R_matrix_LD); data(R_matrix_LDV)
}
\keyword{datasets}
