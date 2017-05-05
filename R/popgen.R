#' Population genetics from genomic data
#'
#' @description This function allows for estimating parameters of population genetics from genomic data. In addition,
#' it also allows estimations considering subpopulations.
#' 
#' @usage popgen(M, subgroups)
#' 
#' @param M object of class \code{matrix}. A (non-empty) matrix of molecular markers, considering the number favorable alleles per loci (0, 1 or 2). Markers must be in columns and individuals in rows.
#' @param subgroups a \code{vector} with information for subgroups or subpopulations.
#' 
#' @details 
#' The matrix of makers is of dimension \eqn{n} x \eqn{p}, in which individuals are in rows and markers in columns.
#' The number of subgroups is user defined  and accepts any data type (\code{character}, \code{integer}, \code{numeric}...) to distinguish subpopulations.
#' These two dataset must have the same sort for rows (genotypes).

#' @return Two lists are returned (\code{general} and \code{bygroup}), one with general information for markers and individuals and another by group (if applicable).
#' 
#' \code{general}: A four-level list 
#' 
#' \itemize{marker}: For each marker it presents the allelic frequency (\eqn{p} and \eqn{q}),
#' Minor Allele Frequency (\eqn{MAF}), expected heterozygosity (\eqn{H_e}), observed
#' heterozygosity (\eqn{H_o}), Nei's Genetic Diversity (\eqn{DG}) and Polymorphism Informative Content(\eqn{PIC}) 
#' 
#' \itemize{genotypes}: it presents observed heterozygosity (\eqn{H_o}), coefficient of inbreeding (\eqn{F_i}) and selfing index (\eqn{S_i})
#'
#' \itemize{population}: The same parameters produced for markers are returned for general population with its mean, lower and upper limits
#'
#' \itemize{Variability}: shows estimates of effective population size (\eqn{Ne}), additive (\eqn{Va}) and dominance (\eqn{Vd}) variances components, and a
#' summary of number of groups, genotypes and markers
#' 
#' \code{bygroups}
#' 
#' Same outputs produced for general it is created for subpopulations or subgroups. Moreover, two more list are presented each with number of exclusive and fixed alleles per group
#'
#' @examples
#' # hybrid maize data
#' data(maize.hyb)
#' x <- popgen(maize.hyb) 
#'
#' # using subpopulations
#' PS<-c(rep(1,25), rep(2,25))
#' x <- popgen(maize.hyb, subgroups=PS)

#' @export
popgen <- function(M, subgroups)
  {
  if(is.null(colnames(M)))
  stop("Colnames is missing")
  
  hasAllMiss <- colSums(is.na(M)) == nrow(M)
  
  if(any(hasAllMiss))
    warning("There are some markers with all data points missing. They were removed from dataset")
  
  Z<-as.matrix(M[, !hasAllMiss])
  
  if(missing(subgroups))
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
      Hobs <- colMeans(M==1, na.rm = T)/2
      Dg <- 1-p^2-q^2
      PIC <- 1-(p^2+q^2)-(2*p^2*q^2)
      propMiss <- colSums(is.na(M))/g
      markers <- round(cbind(p, q, MAF, "He"=Hesp, "Ho"=Hobs, "DG"=Dg, PIC, "Miss" = propMiss), 2)
      markers[is.nan(markers)] <- NA
      markers <- as.data.frame(markers)
      
      
      Hg.obs <- rowMeans(M==1, na.rm = T)/2
      Fi <- 1- Hg.obs/mean(Hesp, na.rm = TRUE)
      Si <- (2*Fi)/(1+Fi)
      
      genotypes <- round(cbind("Ho"=Hg.obs, Fi, Si),2)
      
      meanMrk <- colMeans(markers, na.rm = TRUE)
      rangeMrk <- t(apply(X = markers, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
      
      meanGen <- colMeans(genotypes, na.rm = TRUE)
      rangeGen <- t(apply(X = genotypes, MARGIN = 2, FUN = function(x) range(x, na.rm = TRUE)))
      
      population <- round(rbind(cbind(meanMrk, rangeMrk)[c(6,7,3),], cbind(meanGen, rangeGen)), 2)
      rownames(population) <- c(rownames(population)[1:4], "F", "S")
      colnames(population) <- c("mean", "lower", "upper")
      
      Ne <- 1/(2*mean(Fi))*g
      Va <- sum(2*p*q)
      Vd <- sum((2*p*q)^2)
      variance <- t(round(data.frame(Ne, Va, Vd, "number of genotypes" = g, "number of markers" = m),2))
      colnames(variance) <- ("estimate")
      
      average <- list("markers" = markers, "genotypes" = genotypes, "population" = population, "variability" = variance)
      return(average)
  }
  
  general <- g.of.p(Z)
  
  bygroup <- c("There are no subgroups")
  
  if(nSG > 1){
    bygroup <- lapply(labelSG, function(i) g.of.p(Z[subgroups == labelSG, ]) )
    names(bygroup) <- labelSG
    
    pbyg <- sapply(X = labelSG, FUN = function(x) bygroup[[x]]$markers$p)
    
    for(i in 1:nSG){
    fixed <- pbyg[,i] == 1 | pbyg[,i] == 0
    

    exclus <- colnames(Z)[c(which(pbyg[,i]>0 & apply(as.matrix(pbyg[,-i]==0),MARGIN =  1,FUN = all)),
                            which(pbyg[,i]<1 & apply(as.matrix(pbyg[,-i]==1),MARGIN =  1,FUN = all)))]
    
    fixed <- colnames(Z)[which(fixed)]
    
    bygroup[[labelSG[i]]]$exclusive <- ifelse(length(exclus) == 0,
                                              "There are no exclusive alleles for this group",
                                              exclus)
    bygroup[[labelSG[i]]]$fixed <- ifelse(length(fixed) == 0,
                                          "There are no exclusive alleles for this group",
                                          fixed)
    }
  }
  
    return<-list("general" = general, "bygroup"=bygroup)
}


