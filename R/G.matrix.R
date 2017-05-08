#' Generate Genetic Relationship Matrix
#'
#' @description This function allows to create four different types of Genomic Relationship Matrix (GRM)
#' 
#' @usage G.matrix(M, method=c("VanRaden", "UAR", "UARadj", "GK"), format=c("wide", "long"))
#'        
#' @param M \code{matrix}. Matrix of markers in which \eqn{n} individuals are in rows and \eqn{p} markers in columns.
#' @param method Method for constructing the GRM. Four methods are currently supported. \code{"VanRaden"} indicates the method proposed by Vanraden (2008) for additive
#' genomic relationship and its counterpart for dominance genomic relationship. \code{"UAR"} and \code{"UARadj"} represent methods proposed by Yang et al. (2010) for additive genomic relationship. \code{GK} represents the Gaussian kernel for additive genomic. See 'Detais'
#' @param format \code{character}. Class of object to be returned. \code{"wide"} returns a \eqn{n} x \eqn{n} matrix.
#' \code{"long"} returns the GRM as a table with 3 columns. See 'Details'
#' @details
#' Function G.matrix provides three diferent types of relationship matrix. For \code{VanRaden} method, the relationship matrix is estimated as proposed by Vanraden (2008):
#'  \deqn{G = \frac{XX'}{2\sum p_i(1-p_i)}} where \eqn{X} is the centered marker matrix.  
#' \code{"UAR"} and \code{"UARadj"} represents the genomic relationship matrices proposed by Yang et al. (2010) named by Powell et al. (2010) as Unified Additive Relationship (UAR)
#' and Adjusted UAR, respectively. GRM is then obtained by:
#' A_{jk} = \frac{1}{N}\sum_i{A_{ijk}} = \begin{cases}
#' \dfrac{1}{N} \sum_i{\dfrac{(x_{ij} - 2p_{i})(x_{ik} - 2p_i)}{2p_i(1-p_i)}}, j \neq k \cr
#' 1 +  \dfrac{1}{N} \sum_i{\dfrac{x_{ij}^{2}(1 + 2p_{i})x_{ij} + 2p_i^{2}}{2p_i(1-p_i)}}, j = k
#'\end{cases}  
#' where: p_i is the allele frequency at SNP $i$, x_{ij} is the SNP genotype that takes a value of 0, 1 or 2 if the genotype of the \eqn{j}th
#' individual at SNP \eqn{i} is \eqn{aa}, \eqn{Aa} or \eqn{AA}, respectively.
#' For \code{GK}, the Gaussian kernel is obtained by:
#' \eqn{ K (x_i, x_{i'}) = exp(-d_{ii'}^2)}  
#' Where \eqn{d_{ii'}^2} is the euclidian distance square between two individuals.  
#' The \code{format} argument is the desired output format. For \code{"wide"}, the relationship output produced is in matrix format, with \eqn{n x n} dimension. 
#' If \code{"long"} is the chosen format, the inverse of the relationship matrix is produced and converted to a table. In this case, the upper triangular part of the relationship matrix
#' is changed to a table with 3 columns representing the respective rows, columns and values (Used mainly by ASReml).
#' If the relationship matrix is not positive definite, a near positive definite matrix is created and solved, followed by a warning message.
#' @return The GRM is returned in a matrix or table format.
#' @examples
#' #(1) Additive and dominance relationship matrix 
#' data(maize.hyb)
#' x <- G.matrix(maize.hyb, method="VanRaden", format = "wide")
#' A <- x$Ga
#' D <- x$Gd
#' @references
#' VanRaden, P.M. (2008) Efficient Methods to Compute Genomic Predictions. Journal of Dairy Science, 91:4414-4423 
#' 
#' Yang, J., Benyamin, B., McEvoy, B.P., et al (2010) Common SNPs explain a large proportion of the heritability for human height. Nature Genetics 42:565–569
#'
#'  

#' @export
G.matrix <- function(M, method=c("VanRaden", "UAR", "UARadj", "GK"), format=c("wide", "long")){
  
  if (any(is.na(M)))
    stop("Matrix should not have missing values")
  
  match.arg(method)
  match.arg(format)
  
  if(missing(method))
    stop("Method argument is missing")
     
  if(missing(format))
    stop("Format argument is missing")
  
  N <- nrow(M) 
  m <- ncol(M) 
  p <- colSums(M)/(2*N)

  WWG <- function(M, p){
    w <- scale(x = M, center = T, scale = F)
    
    S <- ((M==2)*1) * -rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
    
    WWl <- w %*% t(w)
    Ga <- WWl/(sum(diag(WWl))/N) + diag(1e-6, nrow(WWl))
    
    SSl <- S %*% t(S)
    Gd <- SSl/(sum(diag(SSl))/N)
    
    return(list(Ga=Ga,Gd=Gd))
  }
  
  UAR <- function(M, p, metho = c("UAR", "UARadj")){
    
    mrep <- scale(x = M, center = T, scale = F)
    X <- (1/m)*(mrep %*% (t(mrep) * (1/(2*p*(1-p)))))
    X[lower.tri(X, diag = T)] <- 0
    X <- X + t(X)
    
    numerator <- M^2 - t(t(M) * (1+2*p)) + matrix(rep(2*p^2, each=N), ncol=m) 
    diag(X) <- 1+(1/m)*colSums(t(numerator) *  (1/(2*p*(1-p))))
    
    if (metho == "UARadj"){
      B <- 1-((6.2*10^-6 + (1/m))/var(c(X)))
      X[lower.tri(X, diag = FALSE)] <- B*X[lower.tri(X, diag = FALSE)]
      X[upper.tri(X, diag = FALSE)] <- B*X[upper.tri(X, diag = FALSE)]  
      diag(X) <- 1 + (B*(diag(X)-1))
    }
      return(X)
  }
  
  toSparse <- function(m){
    comb <- data.frame(row = rep(1:nrow(m), each=nrow(m)),
                       column = rep.int(1:nrow(m), nrow(m)))
    x <- comb[comb$row >= comb$column,]
    x$value <- m[cbind(x$row, x$column)]
    attr(x, "rowNames") <- rownames(m)
    return(x)}
  
  posdefmat <- function(mat){
    #' @importFrom matrixcalc is.positive.definite
    if(is.positive.definite(round(mat, 18))){
      g = solve(mat)
    }else{
      #' @importFrom Matrix nearPD
      g <- solve(nearPD(mat)$mat)
      warning("The matrix was adjusted for the nearest positive definite matrix")
    }
    return(g)
  }
  
  if (method == "VanRaden"){
    Gmat <- WWG(M, p)
    namesG <- names(Gmat)
    if (format == "long"){
      Gmat <- lapply(Gmat, function(x) toSparse(posdefmat(x)))
      
    }
    return(Gmat)
  }
  
  if (method %in% c("UAR", "UARadj")){
    uar <- UAR(M, p, metho = method)
    if (format == "long")
      {uar <- toSparse(posdefmat(uar))}
    return(Ga=uar)
  }
  
  if(method == "GK"){
    w <- scale(x = M, center = T, scale = F)
    D <- as.matrix(dist(w)) ^ 2
    if(quantile(D, 0.05) == 0)
      stop("Was not possible to compute the 5% quantile for the distance matrix")
    GK <- exp(-D / quantile(D, 0.05))
    
    if(format == "long")
      {GK <- toSparse(posdefmat(GK))}
    return(GK)
  }

}
