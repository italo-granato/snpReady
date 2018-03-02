\name{popgen}
\alias{popgen}
\title{Population genetics from genomic data}
\description{
Allows for estimating parameters of population genetics from genomic data. Besides, it also allows the estimate of same parameters considering subpopulations.
}
\usage{
popgen(M, subgroups)
}
\arguments{
\item{M}{Object of class \code{matrix}. A (non-empty) matrix of molecular markers, considering the count of reference alleles per loci (0, 1 or 2). Markers must be in columns and individuals in rows. Missing data should be assigned as \code{NA}}

\item{subgroups}{A \code{vector} with information for subgroups or subpopulations.}
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
heterozygosity (\eqn{H_o}), Nei's Genetic Diversity (\eqn{DG}) and Polymorphism Informative Content(\eqn{PIC})} 
\item{$Genotypes}{It presents observed heterozygosity (\eqn{H_o}) and coefficient of inbreeding (\eqn{F_i})}
\item{$Population}{The same parameters as those for markers except PIC are returned for general population along with lower and upper boundaries}
\item{$Variability}{shows estimates of effective population size (\eqn{Ne}), additive (\eqn{Va}) and dominance (\eqn{Vd}) variances components, and a summary of number of groups, genotypes and markers}

In the presence of \code{subgroups}, the same populational parameters are estimated considering each subpopulation accompanied by its exclusive and fixed alleles. Moreover, a list with the F-statistics for genotypes and markers are exhibited. For genotypes, it shows the statistics considering all subpopulations and a pairwise framework, and for markers loci, the parameters are presented only considering all subpopulations.
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
}
