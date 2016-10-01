#' Generate Genetic Relationship Matrix
#'
#' @description Generate Genomic Relationship Matrix (GRM)
#' 
#' @usage G.matrix(Z, method=c("WW", "UAR", "UARadj"), format=c("long","wide"))
#'        
#' @param \code{Z} \code{matrix}. Matrix of markers in which \eqn{n} individuals are in rows and \eqn{p} markers in columns.
#' @param \code{method} \code{character}. Method for constructing the GRM. Three methods are currently supported. \code{"WW"} indicates the method proposed by Vanraden (2008) for additive
#' and dominant genomic relationship. \code{"UAR"} and \code{"UARadj"} represent methods proposed by Yang et al. (2010) for additive genomic relationship. See 'Detais'
#' @param \code{format} \code{character}. Class of object to be returned. \code{"wide"} returns a \eqn{n} x \eqn{n} matrix.
#' \code{"long"} returns the GRM as a table with 3 columns. See 'Details'
#' @details
#' Function G.matrix provides three diferent types of relationship matrix. For \code{"WW"} method, the relationship matrix is estimated as proposed by Vanraden (2008):
#'  \deqn{G = \frac{ZZ'}{2\sum p_i(1-p_i)}} where: \eqn{Z} is the marker matrix. The SNP genotype takes a value of 0, 1 or 2 if the genotype of the \eqn{j}th
#' individual at SNP \eqn{i} is \eqn{aa}, \eqn{Aa} or \eqn{AA}, respectively.
#' \code{"UAR"} and \code{"UARadj"} represents the genomic relationship matrix proposed by Yang et al. (2010) named by Powell et al. (2010) as Unified Additive Relationship (UAR)
#' and Adjusted UAR, respectively. GRM is then obtained by:
#' \deqn{A_{jk} = \frac{1}{N}\sum_i{A_{ijk}} = \begin{cases}
#' \dfrac{1}{N} \sum_i{\dfrac{(x_{ij} - 2p_{i})(x_{ik} - 2p_i)}{2p_i(1-p_i)}}, j \neq k \cr
#' 1 +  \dfrac{1}{N} \sum_i{\dfrac{x_{ij}^{2}(1 + 2p_{i})x_{ij} + 2p_i^{2}}{2p_i(1-p_i)}}, j = k
#'\end{cases}}
#' where: \eqn{p_i} is the allele frequency at SNP \eqn{i}, \eqn{x_{ij}} is the SNP genotype that takes a value of 0, 1 or 2 if the genotype of the \eqn{j}th
#' individual at SNP \eqn{i} is \eqn{aa}, \eqn{Aa} or \eqn{AA}, respectively.
#' The \code{format} argument is the desired output format. For \code{"wide"}, the relationship output produced is in matrix format, with \eqn{n x n} dimension. 
#' If \code{"long"} is the chosen format, the inverse of the relationship matrix is produced and converted to a table. In this case, the upper triangular part of the relationship matrix
#' is changed to a table with 3 columns representing the respective rows, columns and values (Used mainly by ASReml).
#' If the relationship matrix is not positive definite, a near positive definite matrix is created and solved, followed by a warning message.
#' @return The GRM is returned in a matrix or table format.
#' @examples
#'
#' #(1) Additive and dominance relationship matrix 
#' Z <- data(maize)
#' x <- G.matrix(Z, method="WW", format = "wide")
#' A <- x$Ga
#' D <- x$Gd
#' 

#' @export
G.matrix <- function(Z, method=c("WW", "UAR", "UARadj"), format=c("wide", "long")){
  coded <- unique(as.vector(Z))
  if (any(is.na(match(coded, c(0,1,2)))))
    stop("SNPs must be coded as 0, 1, 2")
  
  if (any(is.na(Z)))
    stop("Matrix must not have missing values")
  
  if(missing(method))
    stop("Method argument is missing")
     
     if(missing(format))
    stop("Format argument is missing")
  
  N <- nrow(Z) 
  n <- ncol(Z) 
  p <- colSums(Z)/(2*nrow(Z))

    WWG <- function(Z, p){
    w <- Z - matrix(rep(2*p, each=nrow(Z)), ncol = ncol(Z))
    
    S <- ((Z==2)*1) * -rep(2*p^2, each=nrow(Z)) + ((Z==1)*1) * rep(2*p*(1-p), each=nrow(Z)) + ((Z==0)*1) * (-rep(2*(1-p)^2, each=nrow(Z)))
    
    WWl <- w %*% t(w)
    Ga <- WWl/sum(2*p*(1-p)) + diag(1e-6, nrow(WWl))
    
    SSl <- S %*% t(S)
    Gd <- SSl/sum((2*p*(1-p))^2)
    
    return(list(Ga=Ga,Gd=Gd))
  }
  
  UAR <- function(Z, p, adj=FALSE){
    
    mrep <- Z - matrix(rep(2*p, each=nrow(Z)), ncol = ncol(Z))
    
    X <- (1/n)*(mrep %*% (t(mrep) * (1/(2*p*(1-p)))))
    numerator <- Z^2 - t(t(Z) * (1+2*p)) + matrix(rep(2*p^2, each=nrow(Z)), ncol=ncol(Z)) 
    diag(X) <- 1+(1/n)*colSums(t(numerator) *  (1/(2*p*(1-p))))
    
    if (adj==TRUE){
      B <- 1-((6.2*10^-6 + (1/n))/var(c(X)))
      X[lower.tri(X, diag = FALSE)] <- B*X[lower.tri(X, diag = FALSE)]
      X[upper.tri(X, diag = FALSE)] <- B*X[upper.tri(X, diag = FALSE)]  
      diag(X) <- 1 + (B*(diag(X)-1))
      return(X)
    }
    else{
      return(X)
    }
  }
  
  
  toSparse <- function(m){
    comb <- data.frame(row = rep(seq(nrow(m)), each=nrow(m)),
                       column = rep.int(seq(nrow(m)), nrow(m)))
    x <- comb[comb$row >= comb$column,]
    x$value <- m[cbind(x$row, x$column)]
    attr(x, "rowNames") <- rownames(m)
    return(x)}
  
  posdefmat <- function(mat){
    #' @importFrom matrixcalc is.positive.definite
    if(is.positive.definite(mat)){
      g = solve(mat)
    }else{
      #' @importFrom Matrix nearPD
      g <- solve(nearPD(mat)$mat)
      warning("The matrix was adjusted for the nearest positive definite matrix")
    }
    return(g)
  }
  
  if (method == "WW" & format == "wide"){
    Gww <- WWG(Z, p)
    return(Gww)
  } 
  if (method == "WW" & format == "long"){
    Gmat <- WWG(Z, p)
    Aww <- toSparse(posdefmat(Gmat$Ga))
    Dww <- toSparse(posdefmat(Gmat$Gd))
    return(list(Ga=Aww, Gd=Dww))
  }
  if (method=="UAR" & format == "wide"){
    uar <- UAR(Z, p)
    return(Ga=uar)
  }
  if (method=="UAR" & format == "long"){
    Gmat <- UAR(Z, p)
    uarsp <- toSparse(posdefmat(Gmat))
    return(Ga=uarsp)
  }
  if (method=="UARadj" & format == "wide"){
    uaradj <- UAR(Z, p, adj = TRUE)
    return(Ga=uaradj)
  }
  if (method=="UARadj" & format == "long"){
    uaradj <- UAR(Z, p, adj = TRUE)
    uaradjsp <- toSparse(posdefmat(uaradj))
    return(Ga=uaradjsp)
  }
}
