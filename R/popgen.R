popgen <- function(M, subgroups=NULL, plot = FALSE)
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
  # ---- plots by subgroup ----
  if(plot){
    main <- c("MAF", "GD", "PIC", "He") 
    
    # barplot
    for(i in 1:4){
      pdf(paste("whole_", main[i],".pdf", sep=""), width = 14, height = 7)
      param <- general$Markers[, main[i]]
      breaks = seq(0, 0.5, 0.1)
      
      switch(main[i],
             "MAF" = {main.tmp = "Minor allele frequency"},
             "GD" = {main.tmp = "Nei's genetic diversity"},
             "PIC" = {main.tmp = "PIC"
             breaks = seq(0, 0.4, 0.1)},
             "He" = {main.tmp = "Expected heterozygosity"})
      
      barplot.pg(x = param, breaks = breaks, main = main.tmp)
      dev.off()
      
      if(nSG > 1){
        pdf(paste("bygroup_", main[i],".pdf", sep=""), width = 14, height = 7)
        ceiling(nSG/2)
        par(mfrow = c(ceiling(nSG/2), 2))
        for(j in 1:nSG){
          param <- bygroup[[j]]$Markers[, main[i] ]
          barplot.pg(x = param, breaks = breaks, main = paste(main.tmp, names(bygroup)[j], sep = " - "))
        }
        dev.off()
      }
    }
    
    # heatmap Fst
    pdf("heatmap_Fst.pdf", width = 14, height = 7)
    heatmap(matFST, scale = "none", Rowv = NA, Colv = NA, cexRow = 0.9, 
            cexCol = 0.9,  main="Fst pairwise", verbose = F)
    dev.off()
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

barplot.pg <- function(x, space = 0.75, width = 4, breaks = NULL, names.arg=NULL, rot_angle = 60,
                       cex = 0.8, xlim = c(0, 50), main = NULL, plotName='general', ext='jpeg')
{
  tmp <- hist(x = x, breaks = breaks, right = T, plot = FALSE)
  nbars <- length(tmp$counts)
  spacet <- c(2, rep(space,nbars-1))
  if(is.null(names.arg))
    names.arg <- paste(breaks[seq(nbars)], "-", breaks[seq(2,nbars+1)])
  
  ptm <- barplot(tmp$counts, beside = T, xlim = xlim, width = width, 
                 ylim = c(0,max(tmp$counts)+10), ylab = "counts",
                 space = spacet, col = "black", axes = F, border = F)
  
  axis(side = 1, at = ptm, labels = F, col.ticks='white')
  axis(side = 2, pos = 5, las = 2, cex.axis = cex)
  text(x = ptm,y = par("usr")[3], labels = names.arg, 
       srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex = cex)
  title(main = main)
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