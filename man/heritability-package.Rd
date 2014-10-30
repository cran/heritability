\name{heritability-package}
\alias{heritability-package}
\alias{heritability}
\docType{package}
\title{
Marker-based estimation of heritability using individual plant or plot data.
}
\description{
The package implements marker-based estimation
of heritability when observations on genetically identical replicates are available.
These can be either observations on individual plants (e.g. in a growth chamber) or plot-level data in a field trial.
The function \code{marker_h2} estimates heritability using a mixed model for the
individual plant or plot data, as proposed in Kruijer \emph{et al.}
For comparison, also mixed-model based estimation using genotypic means (\code{marker_h2_means})
and estimation of repeatability with ANOVA (\code{repeatability}) are implemented.
For illustration the package contains several datasets for the model species \emph{Arabidopsis thaliana}.
}
\details{
\tabular{ll}{
Package: \tab heritability\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-30\cr
License: \tab GPL-3\cr
}
}
\author{
Willem Kruijer
Maintainer: Willem Kruijer <willlem.kruijer@wur.nl>
}
\references{
Kruijer, W. \emph{et al.} (submitted) Marker-based estimation of heritability in immortal
populations.
}
\keyword{ heritability }
\keyword{ mixed model }
\examples{
# A) marker-based estimation of heritability, given individual plant-data
# and a marker-based relatedness matrix:
data(LDV)
data(K_atwell)
# This may take up to 30 sec.
#out1 <- marker_h2(data.vector=LDV$LDV,geno.vector=LDV$genotype,
#                  covariates=LDV[,4:8],K=K_atwell)
#
# B) marker-based estimation of heritability, given genotypic means
# and a marker-based relatedness matrix:
data(means_LDV)
data(R_matrix_LDV)
data(K_atwell)
out2 <- marker_h2_means(data.vector=means_LDV$LDV,geno.vector=means_LDV$genotype,
                        K=K_atwell,Dm=R_matrix_LDV)
#
# C) estimation of repeatability using ANOVA:
data(LDV)
out3 <- repeatability(data.vector=LDV$LDV,geno.vector= LDV$genotype,
                      covariates.frame=as.data.frame(LDV[,3]))
}