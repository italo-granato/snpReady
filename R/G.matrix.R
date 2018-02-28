G.matrix <- function(M, method=c("VanRaden", "UAR", "UARadj", "GK"), format=c("wide", "long"),
                     plot = FALSE){
  
  if (any(is.na(M)))
    stop("Matrix should not have missing values")
  
  if(missing(method))
    stop("Method argument is missing")
     
  N <- nrow(M)
  m <- ncol(M)
  p <- colMeans(M)/2
  
  if(any(p == 0 | p == 1))
    stop("Monomorphic markers are no allowed")

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
  
   switch(method,
         "VanRaden" = {
           Ga <- WWG(M, p)
           },
         "UAR" = {
           Ga <- UAR(M, p, metho = method)
         },
         "UARadj" = {
           Ga <- UAR(M, p, metho = method)
         },
         "GK" = {
           w <- scale(x = M, center = T, scale = T)
           D <- as.matrix(dist(w)) ^ 2
           if(quantile(D, 0.5) == 0)
             stop("Was not possible to compute the 5% quantile for the distance matrix")
           Ga <- exp(-D / quantile(D, 0.05))
         },
         {
           stop("Method selected is not available")
         })
if(format == "long"){
    if(is.list(Ga)){
      G <- lapply(Gmat, function(x) toSparse(posdefmat(x)))
      names(G) <- names(Ga)
    }else{
      G <- toSparse(posdefmat(Ga))
    }
    return(G)
  }else{
    if(plot){
      if(method == "VanRaden"){
        for(i in 1:length(Ga)){
          pdf = pdf(paste("heatmap_", method, names(Ga)[i],".pdf", sep = ""), width = 10, height = 4)
          heatmap(Ga[[i]], scale = "none", Rowv = NA, Colv = NA, cexRow = 0.5, 
                  cexCol = 0.5,  main=paste(method, names(Ga)[i], sep = " "))
          dev.off()
        }
      }else{
        pdf = pdf(paste("heatmap_", method, ".pdf", sep = ""), width = 10, height = 4)
        heatmap(Ga, scale = "none", Rowv = NA, Colv = NA, cexRow = 0.5, 
                cexCol = 0.5,  main=method)
        dev.off()		  
      }
      
      tmp <- Ga
      if(method == "VanRaden")
        tmp <- tmp$Ga
							        
      pckin <- princomp(tmp, cor=TRUE, scores=TRUE)
      #' @importFrom rgl plot3d
      plot3d(pckin$scores[,1:3])
    }
    return(Ga)
  }
}