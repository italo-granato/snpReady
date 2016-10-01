#' @title Preparation of genomic data to perform genomic predictions
#' 
#' @description This function gets genomic data ready to be used in packages or softwares that
#' perform genomic predictions
#' 
#' @usage raw.data(data, frame=c("long","wide"),
#'        hapmap, sweep.sample= 0, call.rate=0.95, maf=0.05,
#'        input=TRUE, outfile=c("012","-101","structure"))
#' 
#' @param data object of class \code{matrix}
#' @param frame \code{character}. Format of genomic data to be inputed. Two formats are currently supported. \code{"long"} is used for objects
#' with sample ID (1st column), marker ID (2nd column), fist allele (3rd column) and second allele (4th column). \code{"wide"} inputes
#' a \eqn{n} x \eqn{m} matrix where markers must be in columns and individuals in rows
#' @param hapmap \code{matrix}. Object with information on SNPs, chromosome and position
#' @param sweep.sample \code{numeric}. Threshold for removing samples from data by missing rate. Samples with missing rate above the defined threshold are
#' removed from dataset.
#' @param call.rate \code{numeric}. Threshold for removing marker by missing genotype rate. SNP with \code{"call rate"} below threshold are removed from dataset.
#' @param maf \code{numeric}. Threshold for removing SNP by minor allele frequency. 
#' @param input \code{logical}. If \code{"TRUE"}, imputation of missing data is performed. See details.
#' @param outfile \code{character}. Type of output to be produced. \code{"012"} outputs matrix coded as 0 to \code{AA}, 1 to \code{Aa} and 2 to \code{aa}. \code{"-101"}
#' presents marker matrix coded as -1, 0 and 1 to \code{aa}, \code{Aa} and \code{AA}, respectively. \code{"structure"} returns a matrix suitable for STRUCTURE Software.
#' For this, each remaining marker is splited in two columns, one for each allele. Nitrogenous bases are then recoded to a number, so A is 1, C is 2, G is 3 and T is 4. 
#' 
#' @details The function allows flexible imputation of genomic data. Data might be in long format with 4 columns or in wide 
#' format where markers are in columns and individuals in rows. Samples and markers can be eliminated based on missing data rate. Markers can also be eliminated based on
#' the frequency of the minor allele. Imputation is carried out through combination of allelic frequency and individual
#' inbreeding coefficient. Hence, for missing values, genotypes are imputed based on their probability of occurrence. This probability
#' depends both on genotype frequency and inbreeding of the individual a specific locus.
#' 
#' @return Returns a properly coded marker matrix output and a report specifying which individuals are removed by \code{sweep.sample} and which markers are removed by \code{"call.rate"}
#' and \code{"maf"}.
#' @seealso # # missing
#' @references # missing
#' @examples
#' data <- data(maize.line)
#' hapmap <- data(hapmap)
#' raw.data(data, frame="table", hapmap, sweep.sample= 0, 
#'          call.rate=0.95, maf=0.05, input=TRUE, outfile="-101")
#'
#'


raw.data <- function(data, frame = c("long","wide"), hapmap, sweep.sample= 0, call.rate=0.95, maf=0.05, input=TRUE, outfile=c("012","-101","structure")) {

  if (call.rate < 0 | call.rate > 1 | maf < 0 | maf > 1 | sweep.sample < 0 | sweep.sample > 1)
    stop("Treshold for call rate, maf and sweep.clean must be between 0 and 1")
  
  if(missing(outfile)) {outfile = "012"}
  
  if(!is.matrix(data))
    stop("Data must be matrix class")
  
  match.arg(frame)
  match.arg(outfile)
  
  if (frame=="long"){
    if(ncol(data)>4)
      stop("For format long, the object must have four columns")
    
    bs <- unique(na.omit(as.vector(data[, 3:4])))
    if(!any(all(bs %in% c("A","C", "G", "T")) | all(bs %in% c("A", "B"))))
      stop("SNPs must be coded as nitrogenous bases (ACGT) or as A and B")
    
    sample.id <- sort(unique(data[,1]))
    snp.name <- sort(unique(data[,2]))
    
    col2row <- function(x, data){
      curId <- data[,1] %in% x
      curSnp <- ifelse(is.na(data[curId, 3]) | is.na(data[curId, 4]), NA, 
                       paste(data[curId, 3], data[curId, 4], sep = ""))
      curPos <- match(snp.name, data[curId, 2])
      if(any(is.na(curPos))){
        vec <- rep(NA,length(snp.name))
        vec[which(!is.na(curPos))] <- curSnp[na.omit(curPos)]
      }else{
        vec <- curSnp[curPos]
        return(vec)}
    }
    
    mbase <- sapply(sample.id, function(x) col2row(x, data))
    colnames(mbase) <- sample.id
    rownames(mbase) <- snp.name
    data <- t(mbase)
  } else{
    bs <- unique(na.omit(as.vector(data)))
    if(!any(all(bs %in% c("A","C", "G", "T")) | all(bs %in% c("A", "B"))))
      stop("SNPs must be coded as nitrogenous bases (ACGT) or as A and B")
  }
  

  count_alleles <- function(x){
    if (all(is.na(x))){
      return(x)
    }else{
      snp.col <- do.call(rbind, strsplit(as.character(x), split = ""))
      ref.allel <- sort(unique(as.vector(snp.col)))
      count <- rowSums(snp.col == ref.allel[1])
      return(count)}}
  
  m <- sapply(as.data.frame(data), function(x) count_alleles(x))
  rownames(m) <- rownames(data)

  miss.freq <- rowSums(is.na(m))/ncol(m)
  
  if(sweep.sample==0){
    m1 <- m
  }else{
    m1 <- m[miss.freq <= sweep.sample,]
    data <- data[miss.freq <= sweep.sample,]
  }
  
  CR <- (colSums(!is.na(m1)) - colSums(is.na(m1)))/colSums(!is.na(m1))
  CR[!is.finite(CR)] <- 0
  
  
  p <- colSums(m1, na.rm = TRUE)/(2*colSums(!is.na(m1)))
  minor <- apply(cbind(p, 1-p), 1, min)
  minor[is.nan(minor)] <- 0
  
  if (call.rate == 0 & maf == 0)
  {
    m2 <- m1
  }else
  {
    position <- which(CR >= call.rate & minor >= maf)
    if (length(position)==0L){
      return(message("All markers were removed. Try again with another treshold for CR and MAF"))}
    m2 <- m1[,position]
    data <- data[, position]
  }
  
  
  if(input==TRUE & call.rate==0 & any(CR==0))
    stop("There's markers with all missing data. Try again using call rate
         different from zero")
  
  if(all(CR==1) | !isTRUE(input))
  {
    m3 <- m2}
  else{
    if (any(miss.freq==1) & sweep.sample==0)
      stop("There's individuals with all missing data. there's no way to do
           imputation. Try again using sweep.sample different from zero")

	CR <- (colSums(!is.na(m1)) - colSums(is.na(m1)))/colSums(!is.na(m1))
    CR[!is.finite(CR)] <- 0
    
    
    p <- colSums(m1, na.rm=TRUE)/(2*colSums(!is.na(m1)))
    minor <- apply(cbind(p,1-p), 1, min)
    minor[is.nan(minor)] <- 0
    
    if (call.rate == 0 & maf == 0)
    {
      m2 <- m1
      }else
        {
        position <- which(CR >= call.rate & minor >= maf)
        if (length(position)==0L){
          return(message("No marker selected. Try again with another treshold for call rate and MAF"))}
        m2 <- m1[,position]
        data <- data[,position]
      }
    
    f <- rowSums(m2!=1, na.rm = TRUE)/rowSums(!is.na(m2))
    f[is.nan(f)] <- 1
    
    samplefp <- function(p, f){
      samp <- sample(c(0,1,2), 1,
                     prob=c(((1-p)^2+((1-p)*p*f)), 
                            (2*p*(1-p)-(2*p*(1-p)*f)), 
                            (p^2+((1-p)*p*f))))
      return(as.integer(samp))}
    
    input.fun <- function(m, p, f){
      icol <- unlist(apply(m, 1, function(x) which(is.na(x))))
      posrow <- apply(m, 1, function(x) sum(is.na(x)))
      irow <- rep(seq(length(posrow)), times=posrow)
      m[cbind(irow, icol)] <- mapply(samplefp, p[icol], f[irow])
      return(m)}

	  if(input==TRUE & call.rate==0 & any(CR==0))
      stop("There are markers with all missing data. Try again using call rate
           different from zero")
    
    if(all(CR==1) | !isTRUE(input))
      {
      m3 <- m2}
    else{
      if (any(miss.freq==1) & sweep.sample==0)
        stop("There are individuals with all missing data. there's no way to do
           imputation. Try again using sweep.sample different from zero")
      
      f <- 1 - (rowSums(m2==1, na.rm = TRUE)/rowSums(!is.na(m2)))
      f[is.nan(f)] <- 1
      
      samplefp <- function(p, f){
        samp <- sample(c(0,1,2), 1,
                       prob=c(((1-p)^2+((1-p)*p*f)), 
                              (2*p*(1-p)-(2*p*(1-p)*f)), 
                              (p^2+((1-p)*p*f))))
        return(as.integer(samp))}
      
      input.fun <- function(m, p, f){
        icol <- unlist(apply(m, 1, function(x) which(is.na(x))))
        posrow <- apply(m, 1, function(x) sum(is.na(x)))
        irow <- rep(seq(length(posrow)), times=posrow)
        m[cbind(irow, icol)] <- mapply(samplefp, p[icol], f[irow])
        return(m)}
      
      m3 <- input.fun(m=m2, p=p[position], f=f)
        }
      
    if (outfile=="012")
      {m4 <- m3}
    
    m3 <- input.fun(m=m2, p=p[position], f=f)
  }
  
  if (outfile=="012")
  {m4 <- m3}
  
  if (outfile=="-101"){
    m4 <- m3 - 1
  }
  
  if(outfile=="structure"){
    m4 <- lapply(as.data.frame(data), function(x){
      curCol <- do.call(rbind, strsplit(as.character(x), split = ""))
      if(all(is.na(curCol))) {curCol <- cbind(curCol, curCol)}
      return(curCol)})
    m4 <- as.matrix(do.call(cbind, m4))
    colnames(m4) <- rep(colnames(data), each=2)
    
    m4 <- chartr("ACGT", "1234", m4)
    m4[is.na(m4)] <- -9
  }
  
  report <- list(paste(sum(minor < maf, na.rm = TRUE), "Markers removed by MAF =", maf, sep = " "),
                 colnames(m)[minor < maf],
                 paste(sum(CR < call.rate), "Markers removed by Call Rate =", call.rate, sep=" "),
                 colnames(m)[CR < call.rate],
                 paste(sum(miss.freq > sweep.sample), "Samples removed =", sweep.sample, sep = " "),
                 rownames(m)[miss.freq > sweep.sample],
                 paste(sum(is.na(m2)), "markers were inputed = ", round(sum(is.na(m2)/(dim(m2)[1]*dim(m2)[2])*100),2), "%")
                 
  )
  
  if(missing(hapmap)){
    storage.mode(m4) <- "numeric"
    return(list(Z.cleaned=m4, report=report))
  } else{
    storage.mode(m4)  <- "numeric"
    hapmap <- hapmap[hapmap[,1] %in% as.matrix(snp.name)[position,],]
    hapmap <- as.matrix(hapmap)
    hapmap <- hapmap[order(as.numeric(hapmap[,2]), as.numeric(hapmap[,3]), na.last = TRUE, decreasing = F),]
    colnames(hapmap) <- c("SNP","Chromosome","Position")
    return(list(Z.cleaned=m4, Hapmap=hapmap, report=report))
  }
}

