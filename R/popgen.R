popgen <-
function(Z, subgroups){
  Z<-as.matrix(Z) # matrix of markers incidence by genotype
  if(missing(subgroups)) {subgroups <- rep(1, nrow(Z))}
  subgroups<-as.factor(subgroups)
  X<-as.data.frame(cbind(subgroups, Z)) # data frame with subgroups and Z
  
g.of.p<-function(X){
  
  M<-X[,-1]
  G<-X[,1]
  m<-ncol(M) # number of markers
  g<-length(G) # number of genotypes
  
  #markers
  f<-apply(M, 2, sum)/(2*g) # allele frequency
  fs<-cbind(f, 1-f) 
  MAF<-apply(fs, 1, min)
  p<-f
  q<-1-f
  Hesp<-2*p*q # expected heterosigosity
  H<-function(M){sum(M[M==1])} #proportion of heterozygotes
  Hobs<-apply(M, 2, H)/(2*g) # observed heterosigosity
  Dg<-1-p^2-q^2 # genetic diversity
  PIC<-1-(p^2+q^2)-(2*p^2*q^2) # Polimorphism Information Contend
  
  markers<-round(data.frame(p, q, MAF, "He"=Hesp, "Ho"=Hobs, "DG"=Dg, PIC),2)
  
  #genotypes
  Hg.obs<-apply(M, 1, H)/(2*m) # observed heterosigosity
  Fi<-1-Hg.obs/mean(Hesp) # endogamy
  Si<-(2*Fi)/(1+Fi) # self index
  
  genotypes<-round(data.frame("Ho"=Hg.obs, Fi, Si),2)
  
  #population
  DG.pop<-c("mean"=mean(Dg), "lower"=range(Dg)[1], "upper"=range(Dg)[2])
  PIC.pop<-c("mean"=mean(PIC), "lower"=range(PIC)[1], "upper"=range(PIC)[2])
  MAF.pop<-c("mean"=mean(MAF), "lower"=range(MAF)[1], "upper"=range(MAF)[2])
 
  Hg.obs.pop<-c("mean"=mean(Hg.obs), "lower"=range(Hg.obs)[1], "upper"=range(Hg.obs)[2])
  Fi.pop<-c("mean"=mean(Fi), "lower"=range(Fi)[1], "upper"=range(Fi)[2])
  Si.pop<-c("mean"=mean(Si), "lower"=range(Si)[1], "upper"=range(Si)[2])
  
  population<-t(round(data.frame("DG"=DG.pop, "PIC"=PIC.pop, "MAF"=MAF.pop, "Ho"=Hg.obs.pop, "F"=Fi.pop, "S"=Si.pop),2))
  
  #variance
  Ne<-1/(2*mean(Fi))*g # effective size
  Va<-sum(2*p*q) # additive variance component
  Vd<-sum((2*p*q)^2) # dominance variance component
  NG<-length(unique(G))
  variance<-t(round(data.frame(Ne, Va, Vd, "number of groups" = NG, "number of genotypes" = g, "number of markers" = m),2))
  colnames(variance)<-("estimate")
  
  #average
  average<-list("markers" = markers, "genotypes" = genotypes, "population" = population, "variability" = variance)  

  return(average)

 }

general<-g.of.p(X)

#by groups
if(length(unique(subgroups)) != 1){
  groups<-split(X, X$subgroups)
  bygroup<-lapply(groups, g.of.p)
  
  pbyg<-matrix(0, ncol(Z), length(unique(subgroups)), dimnames= list(colnames(Z),c(1:length(unique(subgroups)))))

  for(i in 1:length(unique(X$subgroups))){
    pbyg[,i]<-bygroup[[i]]$markers$p
  }

  fallelle <- list()
  exclusive <- list()
  
  for(i in 1:length(unique(X$subgroups))){
    fixed<-pbyg[,i]==1|pbyg[,i]==0
    present<-pbyg[,i]>0

    exclusive[[i]]<-colnames(Z)[c(which(pbyg[,i]>0 & apply(as.matrix(pbyg[,-i]==0), 1, function(x) all(x))),
                                     which(pbyg[,i]<1 & apply(as.matrix(pbyg[,-i]==1), 1, function(x) all(x))))]
    
    
  fallelle[[i]]<-colnames(Z)[fixed==TRUE]
if (length(exclusive[[i]])==0){exclusive[[i]]<-c("there are no exclusive alleles for this group")}
if (length(fallelle[[i]])==0){exclusive[[i]]<-c("there are no fixed alleles for this group")}
  }


  for(i in 1:length(unique(X$subgroups))){
  fixed<-pbyg[,i]==1

  }

return(list("general" = general, "bygroup"=bygroup, "exclusive_alleles"= exclusive, "fixed_alleles" = fallelle))  

  }

else{
  bygroup<-c("there are no subgroups")
  return<-list("general" = general, "bygroup"=bygroup)
 }
}
