\name{popgen}
\alias{popgen}
\title{Population genetics from genomic data}
\description{
This function allows for estimating parameters of population genetics from genomic data. In addition,
it also allows estimations considering subpopulations.
}
\usage{
popgen(M, subgroups)
}
\arguments{
\item{M}{object of class \code{matrix}. A (non-empty) matrix of molecular markers, considering the number favorable alleles per loci (0, 1 or 2). Markers must be in columns and individuals in rows. Missing data should be assign as \code{NA}}

\item{subgroups}{a \code{vector} with information for subgroups or subpopulations.}
}
\details{
The matrix of makers is of dimension \eqn{n} x \eqn{p}, in which individuals are in rows and markers in columns.
The number of subgroups is user defined  and accepts any data type (\code{character}, \code{integer}, \code{numeric}...) to distinguish subpopulations.
These two dataset must have the same sort for rows (genotypes).
}
\value{
  Two-level lists are returned (\code{general} and \code{bygroup}), one with general information for markers and individuals and another by subgroups (if applicable).

For \code{general}, a list containing estimates parameters for
\describe{
\item{$Markers}{For each marker it presents the allelic frequency (\eqn{p} and \eqn{q}),
Minor Allele Frequency (\eqn{MAF}), expected heterozygosity (\eqn{H_e}), observed
heterozygosity (\eqn{H_o}), Nei's Genetic Diversity (\eqn{DG}) and Polymorphism Informative Content(\eqn{PIC})} 
\item{$Genotypes}{It presents observed heterozygosity (\eqn{H_o}), coefficient of inbreeding (\eqn{F_i}) and selfing index (\eqn{S_i})}
\item{$Population}{Treating whole individuals as one population, the same parameters produced for markers are returned for general population with its mean, lower and upper limits}
\item{$Variability}{shows estimates of effective population size (\eqn{Ne}), additive (\eqn{Va}) and dominance (\eqn{Vd}) variances components, and a summary of number of groups, genotypes and markers}

if \code{bygroups} is provided, same output produced for general is created for subgroups. Moreover, two more list are presented each with number of exclusive and fixed alleles for each group
}
\examples{
# hybrid maize data
data(maize.hyb)
x <- popgen(maize.hyb) 

# using subpopulations
PS<-c(rep(1,25), rep(2,25))
x <- popgen(maize.hyb, subgroups=PS)
}
