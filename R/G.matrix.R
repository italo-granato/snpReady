#' Generate realized relationship matrix
#'
#' @param \code{Z} \code{matrix}. Matrix of markers where individuals are in rows and SNP markers are in the columns.
#' @param \code{method} \code{character} method to construct the genomic relationship matrix. Three methods are currently supported. \code{"WW"} indicates method propposed by Vanraden (2008) for additive
#' and dominant genomic relationship. \code{"UAR"} and \code{"UARadj"} represents methods propposed by yang et al. (2010) for additive genomic relationship. See 'Detais'
#' @param \code{frame} \code{character} format of matrix to be returned. \code{"matrix"} is used for value returned as matrix \eqn{n} x \eqn{n}.
#' \code{"column"} is used for return in a table format with 3 columns. See 'Detais'
#' @details
#' Function g.matrix provides three diferent types of relationship matrix. For \code{"ww"} method, the relationship matrix was built as propposed by Vanraden (2008) where
#'  \deqn{G = \frac{ZZ'}{2\sum p_i(1-p_i)}} where \eqn{Z} is the marker matrix coded with classical parameterization. The SNP genotype takes a value of 0, 1 or 2 if the genotype of the \eqn{j}th
#' individual at SNP i \eqn{A_1 A_1}, \eqn{A_1 A_2} or \eqn{A_2 A_2}, respectively was centered for allele frequency.
#' \code{"UAR"} and \code{"UARadj"} represents the genomic relationship matrix propposed by Yang et al. (2010) named by Powell et al. (2010) as unified additive relationship (UAR)
#' and adjusted UAR, respectively. To obtain a genome-wide relationship, they use:
#' \deqn{A_{jk} = \frac{1}{N}\sum_i{A_{ijk}} = \begin{cases}
#' \dfrac{1}{N} \sum_i{\dfrac{(x_{ij} - 2p_{i})(x_{ik} - 2p_i)}{2p_i(1-p_i)}}, j \neq k \cr
#' 1 +  \dfrac{1}{N} \sum_i{\dfrac{x_{ij}^{2}(1 + 2p_{i})x_{ij} + 2p_i^{2}}{2p_i(1-p_i)}}, j = k
#'\end{cases}}
#' Where, \eqn{p_i} is the allele frequency at SNP \eqn{i}, \eqn{x_{ij}} is SNP genotype that takes a value of 0, 1 or 2 if the genotype of the \eqn{j}th
#' individual at SNP \eqn{i} is \eqn{A_1 A_1}, \eqn{A_1 A_2} or \eqn{A_2 A_2}, respectively.
#' The \code{frame} argument is the desired format output. For \code{"matrix"}, the relationship output produced is in matrix format, with dimension \eqn{n x n}. 
#' If \code{"column"} is the format chosen, the inverse of relationship matrix is produced and converted to a table where the upper triangular part of a matrix
#' is converted to a table with columns pointing respectives row, column and value. This format is used mainly by asreml. To solve the inverse of a relationship matrix
#' that should be a positive definite. In case it's not, a near positive definite matrix is created and then solved and a warning is produced.
#' @return The selected matrix is returned in a matrix or table format.
#' @examples
#'
#' #(1) Additive and dominance relationship matrix 
#' Z <- data(maize)
#' x <- G.matrix(Z, method="WW", frame = "matrix")
#' A <- x$Ga
#' D <- x$Gd
#' 


G.matrix <- function(Z, method=c("WW, UAR, UARadj"), frame=c("matrix, column")){
  coded <- unique(as.vector(Z))
  if (any(is.na(match(coded, c(0,1,2)))))
    stop("SNPs must be coded as 0, 1, 2")
  if (any(is.na(Z)))
    stop("matrix must not have missing values")
  N <- nrow(Z) 
  n <- ncol(Z) 
  n.hom <- colSums(Z==2) 
  nhete <- colSums(Z==1)
  n.ind <- nrow(Z) - colSums(is.na(Z))
  p <- (2*n.hom + nhete)/(2*n.ind)
  q <- 1-p
  
  WWG <- function(Z){
    repW <- function(i){
      x <- Z[,i]
      x[x==2] <- 2*q[i]
      x[x==1] <- q[i]-p[i]
      x[x==0] <- (-2)*p[i]
      return(x)}
    w <- sapply(1:ncol(Z), function (i) repW(i))
    repS <- function(j){
      s <- Z[,j]
      s[s==2] <- (-2)*q[j]^2
      s[s==1] <- 2*p[j]*q[j]
      s[s==0] <- (-2)*p[j]^2
      return(s)}
    s <- sapply(1:ncol(Z), function (x) repS(x))
    WWl <- w%*%t(w)
    I <- diag(1e-6, nrow=dim(WWl)[1], ncol=dim(WWl)[2])
    Ga <- WWl/(sum(diag(WWl))/nrow(Z)) + I
    
    SSl <- s%*%t(s)
    Gd <- SSl/(sum(diag(SSl))/nrow(Z))
    
    return(list(Ga=Ga,Gd=Gd))
  }
  
  UAR <- function(Z, adj=FALSE){
    id <- expand.grid(id1 = seq(N), id2 = seq(N))
    Rel <- function (z){
      i=id[z,1]
      j=id[z,2]
      if (j==i){
        y <- 1 +(1/n)*sum((Z[i,]^2-(1+2*p)*Z[i,]+2*p^2)/(2*p*(1-p)))
      }else
      {
        y <- (1/n)*sum(((Z[i,]-2*p)*(Z[j,]-2*p))/(2*p*(1-p)))
      }
      return(y)
    }
    
    UARel <- sapply(1:nrow(id), function(x) Rel(x))
    
    if (adj==TRUE){
      vid <- id[,1]==id[,2]
      (Beta <- 1-((6.2*10^-6 + (1/n))/var(UARel)))
      UARel[vid] <- 1 + (Beta*(UARel[vid]-1))
      UARel[!vid] <- Beta*(UARel[!vid])
    }
    
    Ad <- matrix(UARel, N, N, byrow=TRUE)
    colnames(Ad) <- rownames(Ad) <- rownames(Z)
    return(Ad)
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
  
  
  if (method == "WW" & frame == "matrix"){
    Gww <- WWG(Z)
    return(Gww)
  }
  
  if (method == "WW" & frame == "column"){
    Gmat <- WWG(Z)
    Aww <- toSparse(posdefmat(Gmat$Ga))
    Dww <- toSparse(posdefmat(Gmat$Gd))
    return(list(Ga=Aww, Gd=Dww))
  }
  if (method=="UAR" & frame == "matrix"){
    uar <- UAR(Z)
    return(Ga=uar)
  }
  if (method=="UAR" & frame == "column"){
    Gmat <- UAR(Z)
    uarsp <- toSparse(posdefmat(Gmat))
    return(Ga=uarsp)
  }
  if (method=="UARadj" & frame == "matrix"){
    uaradj <- UAR(Z, adj = TRUE)
    return(Ga=uaradj)
  }
  if (method=="UARadj" & frame == "column"){
    uaradj <- UAR(Z, adj = TRUE)
    uaradjsp <- toSparse(posdefmat(uaradj))
    return(Ga=uaradjsp)
  }
}
