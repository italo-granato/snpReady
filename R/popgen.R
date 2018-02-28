popgen <- function(M, subgroups=NULL)
  {
  if(is.null(colnames(M)))
  stop("Colnames is missing")
  
  hasAllMiss <- colSums(is.na(M)) == nrow(M)
  
  if(any(hasAllMiss))
    warning("There are some markers with all data missing. These markers were removed from dataset")
  
  Z <- as.matrix(M[, !hasAllMiss])
  
  if(is.null(subgroups))
    subgroups <- 1
  
  labelSG <- unique(subgroups)
  nSG <- length(labelSG)

  general <- g.of.p(Z)
  
  bygroup <- c("There are no subgroups")
  
  if(nSG > 1){
    bygroup <- lapply(labelSG, function(i) g.of.p(Z[subgroups == i, ]) )
    names(bygroup) <- labelSG
    
    pbyg <- sapply(X = as.vector(labelSG), FUN = function(x) bygroup[[x]]$Markers$p)
    ## Exclusive alleles
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
    
    # ----- F statistics -----
    ngroups <- as.vector(table(subgroups))
    Hig <- matrix(sapply(bygroup, function(x) x$Population["Ho","mean"]), ncol = nSG)
    His <- sapply(bygroup, function(x) x$Markers[,"Ho"])
    Hss <- sapply(bygroup, function(x) x$Markers[,"He"])
    Hsg <- matrix(colMeans(Hss), ncol = nSG)
    Ht <- matrix(general$Markers[,"He"], ncol = 1, dimnames = list(rownames(general$Markers), NULL))
    
    Fstatsg <- F.stats(Hi = Hig, Hs = Hsg, Ht = mean(Ht), ngroups = ngroups)
    Fstatsm <- F.stats(Hi = His, Hs = Hss, Ht = Ht, ngroups = ngroups)
    
    # ------ Fstats pairwise ------
    pw <- combn(x = labelSG, m = 2)
	matFST <- matrix(0, nrow = nSG, ncol = nSG, dimnames = list(labelSG, labelSG))
    
    Fstspw <- round(data.frame("Fis" = numeric(ncol(pw)+1),
                               "Fst" = numeric(ncol(pw)+1),
                               "Fit" = numeric(ncol(pw)+1),
                               row.names = c("All_Pop", paste(pw[1,], pw[2,], sep = "-") )), 3)
    Fstspw[1,] <- Fstatsg
    
    for(i in 1:ncol(pw)){
      sel <- labelSG %in% pw[,i]
      nsbg <- ngroups[sel] 
      Hisg <- Hig[,sel, drop=FALSE]
      Hssg <- Hsg[,sel, drop=FALSE]
      Fstspw[i+1,] <- F.stats(Hi = Hisg, Hs = Hssg, Ht = mean(Ht), ngroups = nsbg)
	  matFST[pw[1,i], pw[2,i]] <- matFST[pw[2,i], pw[1,i]] <- Fstspw[i+1, 2]						   
    }
    
    Fstats <- list("Genotypes" = Fstspw, "Markers" = Fstatsm)
    
    bygroup <- c(bygroup, list("F.stats" = Fstats))
  }
  
    out <- list("whole" = general, "bygroup" = bygroup)
    return(out)
}

 g.of.p <- function(M){
    m<-ncol(M)
    g<-nrow(M)
    
    p <- colMeans(M, na.rm = T)/2
    fs <- cbind(p, 1-p)
    MAF <- apply(fs, 1, min)
    q <- 1-p
    Hesp <- 2*p*q
    Hobs <- colMeans(M == 1, na.rm = T)
    Dg <- 1- p^2 - q^2
    PIC <- 1-(p^2 + q^2) - (2*p^2*q^2)
    propMiss <- colSums(is.na(M))/g
	counts <- matrix(NA, nrow = ncol(M), ncol = 3, dimnames = list(colnames(M), c(0,1,2)))
  for(i in 1:3){
    counts[,i] <- colSums(M == colnames(counts)[i], na.rm = T)
  }
  hwetest <- chiS(counts = counts)																					
    markers <- cbind(round(cbind(p, q, MAF, "He" = Hesp, "Ho" = Hobs, "GD" = Dg, PIC, "Miss" = propMiss), 2),
                   hwetest)
    markers[is.nan(markers)] <- NA
    markers <- as.data.frame(markers)
    
    
    Hg.obs <- rowMeans(M == 1, na.rm = TRUE)
    Fi <- inbreeding.fun(mat = M, p = p)
    
    
    genotypes <- round(cbind("Ho" = Hg.obs, "Fi" = Fi),2)
    
    meanMrk <- colMeans(markers, na.rm = TRUE)
    rangeMrk <- t(apply(X = markers, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
    
    meanGen <- colMeans(genotypes, na.rm = TRUE)
    rangeGen <- t(apply(X = genotypes, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
    
    population <- round(rbind(cbind(meanMrk, rangeMrk)[c(6,7,3),], cbind(meanGen, rangeGen)), 2)
    rownames(population) <- c(rownames(population)[1:4], "F")
    colnames(population) <- c("mean", "lower", "upper")
    
    Ne <- (1/(2*mean(Fi)))*g
    Va <- sum(2*p*q)
    Vd <- sum((2*p*q)^2)
    variance <- t(round(data.frame(Ne, Va, Vd, "number of genotypes" = g, "number of markers" = m),2))
    colnames(variance) <- ("estimate")
    
    average <- list("Markers" = markers, "Genotypes" = genotypes, "Population" = population, "Variability" = variance)
    return(average)
  }
F.stats <- function(Hi, Hs, Ht, ngroups){
  n.harm <- matrix(ngroups/sum(ngroups), nrow = 1)
  
  if(nrow(Hi) > 1){
    n.harm <- n.harm[rep(1, nrow(Hi)),]
  }
  
  Hs.pop <- rowSums(Hs * n.harm)
  Hi.pop <- rowSums(Hi * n.harm)
  
  Fis.pop <- (Hs.pop - Hi.pop)/Hs.pop
  Fst.pop <- (Ht - Hs.pop)/Ht
  Fit.pop <- (Ht - Hi.pop)/Ht
  FST.pop <- round(data.frame("Fis" = Fis.pop, "Fst" = Fst.pop, "Fit" = Fit.pop), 3) 
  rownames(FST.pop) <- rownames(Ht)
  return(FST.pop)
}
inbreeding.fun <- function(mat, p){
  nOHom <- rowSums(mat != 1, na.rm = T)
  nEHom <- 1 - (2* p *(1-p))
  EH <- as.vector(round(nEHom %*%  t(!is.na(mat))))
  Fi <- round((nOHom - EH)/(rowSums(!is.na(mat)) - EH), 3)
  return(Fi)
}

chiS <- function(counts){
  p <- ((2 * counts[,"2"]) + counts[,"1"])/(2*rowSums(counts))
  
  Expfr <- Vectorize(FUN = function(p){
    return(c("0" = (1-p)**2, "1" = 2*p*(1-p), "2" = p**2))
  })
  
  E <- t(Expfr(p)) * rowSums(counts)
  
  sumChi <- rowSums((E - counts)^2/E)
  
  pvalue <- pchisq(sumChi, 1, lower.tail = FALSE)
  return(cbind("pval" = pvalue, "chiSq" = round(sumChi, 3)))
}