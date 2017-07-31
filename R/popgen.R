popgen <- function(M, subgroups=NULL)
  {
  if(is.null(colnames(M)))
  stop("Colnames is missing")
  
  hasAllMiss <- colSums(is.na(M)) == nrow(M)
  
  if(any(hasAllMiss))
    warning("There are some markers with all data points missing. They were removed from dataset")
  
  Z <- as.matrix(M[, !hasAllMiss])
  
  if(is.null(subgroups))
    subgroups <- 1
  
  labelSG <- unique(subgroups)
  nSG <- length(labelSG)

  g.of.p<-function(M){

      m<-ncol(M)
      g<-nrow(M)
      
      p <- colMeans(M, na.rm = T)/2
      fs <- cbind(p, 1-p)
      MAF <- apply(fs, 1, min)
      q <- 1-p
      Hesp <- 2*p*q
      Hobs <- colMeans(M==1, na.rm = T)
      Dg <- 1-p^2-q^2
      PIC <- 1-(p^2 + q^2) - (2*p^2*q^2)
      propMiss <- colSums(is.na(M))/g
      markers <- round(cbind(p, q, MAF, "He"=Hesp, "Ho"=Hobs, "DG"=Dg, PIC, "Miss" = propMiss), 2)
      markers[is.nan(markers)] <- NA
      markers <- as.data.frame(markers)
      
      
      Hg.obs <- rowMeans(M == 1, na.rm = T)
      Fi <- 1- Hg.obs/mean(Hesp, na.rm = TRUE)
      Si <- (2*Fi)/(1+Fi)
      
      genotypes <- round(cbind("Ho" = Hg.obs, Fi, Si),2)
      
      meanMrk <- colMeans(markers, na.rm = TRUE)
      rangeMrk <- t(apply(X = markers, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
      
      meanGen <- colMeans(genotypes, na.rm = TRUE)
      rangeGen <- t(apply(X = genotypes, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
      
      population <- round(rbind(cbind(meanMrk, rangeMrk)[c(6,7,3),], cbind(meanGen, rangeGen)), 2)
      rownames(population) <- c(rownames(population)[1:4], "F", "S")
      colnames(population) <- c("mean", "lower", "upper")
      
      Ne <- 1/(2*mean(Fi)*g)
      Va <- sum(2*p*q)
      Vd <- sum((2*p*q)^2)
      variance <- t(round(data.frame(Ne, Va, Vd, "number of genotypes" = g, "number of markers" = m),2))
      colnames(variance) <- ("estimate")
      
      average <- list("Markers" = markers, "Genotypes" = genotypes, "Population" = population, "Variability" = variance)
      return(average)
  }
  
  general <- g.of.p(Z)
  
  bygroup <- c("There are no subgroups")
  
  if(nSG > 1){
    bygroup <- lapply(labelSG, function(i) g.of.p(Z[subgroups == i, ]) )
    names(bygroup) <- labelSG
    
    pbyg <- sapply(X = as.vector(labelSG), FUN = function(x) bygroup[[x]]$Markers$p)
    
    for(i in 1:nSG){
    fixed <- pbyg[,i] == 1 | pbyg[,i] == 0
    

    exclus <- colnames(Z)[c(which(pbyg[,i]>0 & apply(as.matrix(pbyg[,-i]==0),MARGIN =  1,FUN = all)),
                            which(pbyg[,i]<1 & apply(as.matrix(pbyg[,-i]==1),MARGIN =  1,FUN = all)))]
    
    fixed <- colnames(Z)[which(fixed)]
    
    if(length(exclus) == 0){
      excl <- "There are no exclusive alleles for this group"
    }else{
      excl <- exclus
    }
    
    bygroup[[labelSG[i]]]$exclusive <- excl
    
    if(length(fixed) == 0){
      fix.g <- "There are no fixed alleles for this group"
    }else{
      fix.g <- fixed
    }
    
    bygroup[[labelSG[i]]]$fixed <- fix.g
    
    }
  }
  
    return<-list("whole" = general, "bygroup" = bygroup)
}


