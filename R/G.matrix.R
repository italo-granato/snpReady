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
  
  if(any(p == 0 | p == 1))
    stop("Monomorphic markers are no accepted")

  WWG <- function(M, p){
    w <- scale(x = M, center = T, scale = F)
    
    S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
    
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
    w <- scale(x = M, center = T, scale = T)
    D <- as.matrix(dist(w)) ^ 2
    if(quantile(D, 0.5) == 0)
      stop("Was not possible to compute the 50% quantile for the distance matrix")
    GK <- exp(-D / quantile(D, 0.05))
    
    if(format == "long")
      {GK <- toSparse(posdefmat(GK))}
    return(GK)
  }

}
