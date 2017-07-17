\name{raw.data}
\alias{raw.data}
\title{Preparation of genomic data}
\description{
This function gets genomic data ready to be used in packages or softwares that
perform genomic predictions
}
\usage{
raw.data(data, frame = c("long","wide"), hapmap=NULL, base=TRUE, sweep.sample=1,
                     call.rate=0.95, maf=0.05, imput=TRUE,
                     imput.type=c("wright", "mean"), outfile=c("012","-101","structure"))
}
\arguments{
\item{data}{object of class \code{matrix}}

\item{frame}{Format of genomic data to be imputed. Two formats are currently supported. \code{"long"} is used for objects
with sample ID (1st column), marker ID (2nd column), fist allele (3rd column) and second allele (4th column). While \code{"wide"} is for objects with dimension \eqn{n \times m} where markers must be in columns and individuals in rows}

\item{hapmap}{\code{matrix}. Object with information on SNPs, chromosome and position}

\item{base}{\code{logical}. Are genotypes coded as nitrogenous bases? if \code{TRUE}, data are converted to numeric. If \code{FALSE}, it follows to clean up}

\item{sweep.sample}{\code{numeric}. Threshold for removing samples from data by missing rate. Samples with missing rate above the defined threshold are
removed from dataset}

\item{call.rate}{\code{numeric}. Threshold for removing marker by missing genotype rate. SNP with \code{call rate} below threshold are removed from dataset. Must be between {0, 1}}

\item{maf}{Threshold for removing SNP by minor allele frequency. Must be between {0, 1}}

\item{imput}{Should imputation of missing data be performed?. Default is \code{TRUE} }

\item{imput.type}{Type of imputation. It can be by "wright" or "mean". See \code{details} }

\item{outfile}{\code{character}. Type of output to be produced. \code{"012"} outputs matrix coded as 0 to \code{AA}, 1 to \code{Aa} and 2 to \code{aa}. \code{"-101"}
presents marker matrix coded as -1, 0 and 1 to \code{aa}, \code{Aa} and \code{AA}, respectively. \code{"structure"} returns a matrix suitable for STRUCTURE Software.
For this, each individual is splited in two rows, one for each allele. Nitrogenous bases are then recoded to a specific number, so A is 1, C is 2, G is 3 and T is 4. This format is only acceptable if base are \code{TRUE}}
}
\details{
The function allows flexible input of genomic data. Data might be in long format with 4 columns or in wide 
format where markers are in columns and individuals in rows. Both numeric and nitrogenous bases are accepted. Samples and markers can be eliminated based on missing data rate. Markers can also be eliminated based on
the frequency of the minor allele.  
Two methods of imputation are currently implemented. One is carried out through combination of allele frequency and individual inbreeding coefficient estimated from markers. 
\deqn{ 
p(x_{ij}) = \left \{
          \begin{array}{ll}
          0 = (1 - p_j)^2 + p_j (1 - p_j) F_i \\
          1 = 2 p_j (1 - p_j) - 2 p_j (1 - p_j) F_i\\
          2 = p_j^2 + p_j (1 - p_j) F_i        
                \end{array}
              \right.
}

Hence, for missing values, genotypes are imputed based on their probability of occurrence. This probability depends both on genotype frequency and inbreeding of the individual a specific locus.
The second method is based on mean of SNP. Thus, each missing point in a SNP \eqn{j} is replaced by mean of SNP \eqn{j}  
\deqn{x_{ij} = 2p_j}
  
   
}
\value{
Returns a properly coded marker matrix output and a report specifying which individuals are removed by \code{sweep.sample} and which markers are removed by \code{"call.rate"}
and \code{maf}.
}
\examples{
data(maize.line)
M <- as.matrix(maize.line)
raw.data(M, frame="long", base=TRUE, sweep.sample= 0.8, 
         call.rate=0.95, maf=0.05, imput=TRUE, outfile="-101")

}
\references{
# missing
}
\seealso{
# # missing
}
