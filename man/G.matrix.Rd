\name{G.matrix}
\alias{G.matrix}
\encoding{UTF-8}
\title{Estimation of Genomic Relationship Matrix}
\usage{
G.matrix(M, method=c("VanRaden", "UAR", "UARadj", "GK"), format=c("wide", "long"), 
         plot = FALSE)
}
\arguments{
\item{M}{\code{matrix}. Matrix of markers in which \eqn{n} individuals are in rows and \eqn{p} markers in columns. This matrix do not need to be centered.}

\item{method}{Method to built the GRM. Four methods are currently supported. \code{"VanRaden"} indicates the method proposed by Vanraden (2008) for additive
genomic relationship and its counterpart for dominance genomic relationship. \code{"UAR"} and \code{"UARadj"} are methods proposed by Yang et al. (2010) for additive genomic relationship. \code{"GK"} represents the Gaussian kernel for additive genomic. See \code{Details}}

\item{format}{Type of object to be returned. \code{wide} returns a \eqn{n \times n} matrix. \code{long} returns the low diagonal from GRM as a table with 3 columns. See \code{Details}}
\item{plot}{If \code{TRUE}, a graphical output is produced. See \code{Details}}
}
\value{
It returns the GRM. If the method is \code{VanRaden}, additive and dominance matrices are produced. Otherwise, only the additive form. 
If \code{plot} is \code{TRUE} a heat map of the pairwise relationship is save as pdf into the working directory . Also, a 3D plot with the three first principal components is generated. 
}
\description{
It generates four different types of Genomic Relationship Matrix (GRM)
}
\details{
G.matrix provides four different types of relationship matrix. The \code{VanRaden} represents the relationship matrix estimated as proposed by Vanraden (2008):
 \deqn{G = \frac{XX'}{trace(XX')/n}}
\eqn{X} is the centered marker matrix. For any marker locus \eqn{i}, \eqn{x_i = m_i - 2p_{i}} where \eqn{m_i} is the vector of SNP genotypes coded as allele couting (0, 1 and 2).

\code{UAR} is the genomic relationship matrices proposed by Yang et al. (2010) and named as Unified Additive Relationship (UAR). This matrix is then obtained by
\deqn{
G_{UAR} = A_{jk} = \frac{1}{N} \sum_i{A_{ijk}} = \left \{
\begin{array}{ll}
\frac{1}{N} \sum_i{\frac{(x_{ij} - 2p_{i})(x_{ik} - 2p_i)}{2p_i(1-p_i)}}, j \neq k \\
1 +  \frac{1}{N} \sum_i{\frac{x_{ij}^{2}(1 + 2p_{i})x_{ij} + 2p_i^{2}}{2p_i(1-p_i)}}, j = k
\end{array}
\right.
}  
where \eqn{p_i} is the allele frequency at SNP \eqn{i} and \eqn{x_{ij}} is the SNP genotype that takes a value of 0, 1 or 2 for the genotype of the \eqn{j^{th}} individual at SNP \eqn{i}.
The same authors proposed an adjustment in the original UAR matrix (\code{UARadj}) to reduce the bias in estimation of variance in the relationship in causal loci. Thus:
\deqn{
G_{UARadj} = \left \{
\begin{array}{ll}
 \beta A_{jk}, j \neq k \\
1 + \beta(A_{jk} - 1), j = k
\end{array}
\right.
}  
where \eqn{\beta = 1 - frac{c + 1/N}{var(A_{jk}}}, c is a constant dependent on MAF of causal variants assumed as 0. 
\code{GK} represents the Gaussian kernel, obtained by
\deqn{ K (x_i, x_{i'}) = \frac{exp(-d_{ii'}^2)}{quantile(d^2, 0.5)}}
where \eqn{d_{ii'}^2} is the square of euclidian distance between two individuals
The \code{format} argument is the desired output format. For \code{"wide"}, the relationship output produced is in matrix format, with \eqn{n \times n} dimension. 
If \code{"long"} is the chosen format, the inverse of the relationship matrix is computed and converted to a table. In this case, the low triangular part of the relationship matrix
is changed to a table with three columns representing the respective rows, columns, and values (Used mainly by ASReml)

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
Pérez-Elizalde, S.,Cuevas, J.; Pérez-Rodríguez, P.; Crossa, J. (2015) Selection of The Bandwidth Parameter in a Bayesian Kernel Regression Model for Genomic-Enabled Prediction. J Agr Biol Envir S, 20-4:512-532
  
Yang, J., Benyamin, B., McEvoy, B.P., et al (2010) Common SNPs explain a large proportion of the heritability for human height. Nature Genetics 42:565-569

VanRaden, P.M. (2008) Efficient Methods to Compute Genomic Predictions. Journal of Dairy Science, 91:4414-4423

}