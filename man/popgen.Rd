\name{popgen}
\alias{popgen}
\title{Population genetics from genomic data}
\description{
Allows for estimating parameters of population genetics from genomic data. Besides, it also allows the estimate of same parameters considering subpopulations.
}
\usage{
popgen(M, subgroups, plot = FALSE)
}
\arguments{
\item{M}{Object of class \code{matrix}. A (non-empty) matrix of molecular markers, considering the count of reference alleles per loci (0, 1 or 2). Markers must be in columns and individuals in rows. Missing data should be assigned as \code{NA}}

\item{subgroups}{A \code{vector} with information for subgroups or subpopulations.}
\item{plot}{If \code{TRUE}, a graphical output is produced. See \code{details} }
}
\details{
The number of subgroups is defined by the user and accepts any data type (\code{character}, \code{integer} ...) to distinguish subpopulations.
These two inputs must have the same sort for rows (genotypes).

}
\value{
  Two-level lists are returned (\code{whole} and \code{bygroup}), one with general information for markers and individuals and another by subgroups (if applicable).

For \code{whole}, a list containing estimates parameters for
\describe{
\item{$Markers}{For each marker it presents the allelic frequency (\eqn{p} and \eqn{q}),
Minor Allele Frequency (\eqn{MAF}), expected heterozygosity (\eqn{H_e}), observed
heterozygosity (\eqn{H_o}), Nei's Genetic Diversity (\eqn{DG}) and Polymorphism Informative Content(\eqn{PIC}), proportion of missing (\eqn{Miss}), \eqn{\chi^2} statistic for the Hardy-Weinberg equilibrium test and its pvalue} 
\item{$Genotypes}{It presents observed heterozygosity (\eqn{H_o}) and coefficient of inbreeding (\eqn{F_i}) estimated as excess of homozygous relative to the expected (Keller et al. (2011))}
\item{$Population}{The same parameters as those for markers except PIC are returned for general population along with lower and upper boundaries}
\item{$Variability}{shows estimates of effective population size (\eqn{Ne}), additive (\eqn{Va}) and dominance (\eqn{Vd}) variances components, and a summary of number of groups, genotypes and markers}

In the presence of \code{subgroups}, the same populational parameters are estimated considering each subpopulation accompanied by its exclusive and fixed alleles. Moreover, a list with the F-statistics (\eqn{F_IT}, \eqn{F_IS} and \eqn{F_ST}) for genotypes and markers are exhibited. For genotypes, it shows the statistics considering all subpopulations and a pairwise framework, and for markers loci, the parameters are presented only considering all subpopulations.

The \code{plot} produces a histogram for the estimates of \eqn{MAF}, \eqn{GD}, \eqn{PIC} and \eqn{He} considering the whole population and subpopulations, when it is available. Also, a heat map of the pairwise \eqn{F_ST} between populations is produced.

}
}
\examples{
# hybrid maize data
data(maize.hyb)
x <- popgen(maize.hyb) 

# using subpopulations
PS<-c(rep(1,25), rep(2,25))
x <- popgen(maize.hyb, subgroups=PS)
}

\references{
Weir, B.S. and C.C. Cockerham. (1984). Estimating F-Statistics for the Analysis of Population Structure. Evolution 38: 1358-1370. doi:10.2307/2408641.

Keller M.C., Visscher P.M., Goddard M.E. (2011) Quantification of inbreeding due to distant ancestors and its detection using dense single nucleotide polymorphism data. Genetics 189:237-249. doi: 10.1534/genetics.111.130922
}
