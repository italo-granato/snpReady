\name{G.matrix}
\alias{G.matrix}
\encoding{UTF-8}
\title{Generate Genetic Relationship Matrix (GRM)}
\usage{
G.matrix(M, method=c("VanRaden", "UAR", "UARadj", "GK"), format=c("wide", "long"))
}
\arguments{
\item{M}{\code{matrix}. Matrix of markers in which \eqn{n} individuals are in rows and \eqn{p} markers in columns.}

\item{method}{Method for constructing the GRM. Four methods are currently supported. \code{"VanRaden"} indicates the method proposed by Vanraden (2008) for additive
genomic relationship and its counterpart for dominance genomic relationship. \code{"UAR"} and \code{"UARadj"} represent methods proposed by Yang et al. (2010) for additive genomic relationship. \code{GK} represents the Gaussian kernel for additive genomic. See \code{Details}}

\item{format}{Type of object to be returned. \code{wide} returns a \eqn{n \times n} matrix. \code{long} returns the low diagonal from GRM as a table with 3 columns. See \code{Details}}
}
\value{
The GRM is returned in a matrix or table format
}
\description{
This function allows to create four different types of Genomic Relationship Matrix (GRM)
}
\details{

G.matrix provides four diferent types of relationship matrix
\describe{
\item{\code{VanRaden}}{ Relationship matrix is estimated as proposed by Vanraden (2008)
 \deqn{G = \frac{XX'}{2\sum p_i(1-p_i)}} where \eqn{X} is the centered marker matrix}

\item{\code{UAR} and \code{UARadj}}{The genomic relationship matrices proposed by Yang et al. (2010) named by Powell et al. (2010) as Unified Additive Relationship (UAR)
and Adjusted UAR, respectively. GRM is then obtained by
\deqn{A_{jk} = \frac{1}{N}\sum_i{A_{ijk}} = \begin{cases}
\dfrac{1}{N} \sum_i{\dfrac{(x_{ij} - 2p_{i})(x_{ik} - 2p_i)}{2p_i(1-p_i)}}, j \neq k \cr
1 +  \dfrac{1}{N} \sum_i{\dfrac{x_{ij}^{2}(1 + 2p_{i})x_{ij} + 2p_i^{2}}{2p_i(1-p_i)}}, j = k
\end{cases}}  
where \eqn{p_i} is the allele frequency at SNP \eqn{i} and \eqn{x_{ij}} is the SNP genotype that takes a value of 0, 1 or 2 for the genotype of the \eqn{j^{th}} individual at SNP \eqn{i}}  

\item{\code{GK}}{the Gaussian kernel is obtained by
\deqn{ K (x_i, x_{i'}) = \dfrac{exp(-d_{ii'}^2)}{quantile(d^2, 0.05)}}
where \eqn{d_{ii'}^2} is the square of euclidian distance between two individuals}
}
The \code{format} argument is the desired output format. For \code{"wide"}, the relationship output produced is in matrix format, with \eqn{n \times n} dimension. 
If \code{"long"} is the chosen format, the inverse of the relationship matrix is computed and converted to a table. In this case, the low triangular part of the relationship matrix
is changed to a table with 3 columns representing the respective rows, columns and values (Used mainly by ASReml)

If the relationship matrix is not positive definite, a near positive definite matrix is created and solved, followed by a warning message.
}
\examples{
#(1) Additive and dominance relationship matrix 
data(maize.hyb)
x <- G.matrix(maize.hyb, method = "VanRaden", format = "wide")
A <- x$Ga
D <- x$Gd
}
\references{
VanRaden, P.M. (2008) Efficient Methods to Compute Genomic Predictions. Journal of Dairy Science, 91:4414-4423 
Yang, J., Benyamin, B., McEvoy, B.P., et al (2010) Common SNPs explain a large proportion of the heritability for human height. Nature Genetics 42:565-569
}