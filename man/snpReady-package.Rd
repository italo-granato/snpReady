\name{snpReady-package}
\alias{snpReady-package}
\alias{snpReady}
\docType{package}
\title{
\packageTitle{snpReady}
}
\description{
\packageDescription{snpReady}
}
\details{
The first two function, \code{\link{raw.data}} and \code{\link{G.matrix}}, was created to clean, summarize and prepare genomic datasets to Genome Selection and Genome Association analysis. \code{\link{popgen}}, the third function, it is possible to estimate population genetic parameters for whole population or by subgroups, per genotypes or per loci. Furthermore, fixed alleles or even exclusive markers can be identify using popgen().

}
\author{
Italo S C Granato and Roberto Fritsche-Neto

Maintainer: Italo S C Granato \email{italo.granato@gmail.com}
}
\references{
VanRaden, P.M. (2008) Efficient Methods to Compute Genomic Predictions. J Dairy Sci 91: 4414-4423.

Yang, J., B. Benyamin, B.P. McEvoy, S. Gordon, A.K. Henders, D.R. Nyholt, et al. (2010) Common SNPs explain a large proportion of the heritability for human height. Nature Genet. 569: 565--569.

}

\keyword{ package; genetics }

\examples{
#raw.data
data(maize.line)
M <- as.matrix(maize.line)
raw.data(M, frame="long", base=TRUE, sweep.sample= 0.8, 
          call.rate=0.95, maf=0.05, input=TRUE, outfile="-101")

#G.matrix
data(maize.hyb)
G.matrix(maize.hyb, method="VanRaden", format = "wide")

#popgen
data(maize.hyb)
x <- popgen(maize.hyb) 
}
