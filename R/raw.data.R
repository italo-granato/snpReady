#' @title Preparation of genomic data
#' 
#' @description This function gets genomic data ready to be used in packages or softwares that
#' perform genomic predictions
#' 
#' @usage raw.data(data, frame=c("long","wide"), hapmap,
#'        bases=TRUE, sweep.sample= 1, call.rate=0.95, maf=0.05,
#'        input=TRUE, outfile=c("012","-101","structure"))
#' 
#' @param data object of class \code{matrix}
#' @param frame \code{character}. Format of genomic data to be inputed. Two formats are currently supported. \code{"long"} is used for objects
#' with sample ID (1st column), marker ID (2nd column), fist allele (3rd column) and second allele (4th column). \code{"wide"} inputes
#' a \eqn{n} x \eqn{m} matrix where markers must be in columns and individuals in rows
#' @param hapmap \code{matrix}. Object with information on SNPs, chromosome and position
#' @param base \code{logical}. Are gentype coded as nitrogenous bases? if \code{TRUE}, data are converted to numeric. If \code{FALSE}, it follows to clean up.
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
#' format where markers are in columns and individuals in rows. Both numeric and nitrogenous bases are accepted. Samples and markers can be eliminated based on missing data rate. Markers can also be eliminated based on
#' the frequency of the minor allele. Imputation is carried out through combination of allelic frequency and individual
#' inbreeding coefficient. Hence, for missing values, genotypes are imputed based on their probability of occurrence. This probability
#' depends both on genotype frequency and inbreeding of the individual a specific locus.
#' 
#' @return Returns a properly coded marker matrix output and a report specifying which individuals are removed by \code{sweep.sample} and which markers are removed by \code{"call.rate"}
#' and \code{"maf"}.
#' @seealso # # missing
#' @references # missing
#' @examples
#' data(maize.line)
#' raw.data(maize.line, frame="table", base=TRUE, hapmap, sweep.sample= 1, 
#'          call.rate=0.95, maf=0.05, input=TRUE, outfile="-101")
#' 
#'


#' @export
raw.data <- function(data, frame = c("long","wide"), hapmap, base=TRUE, sweep.sample= 1, call.rate=0.95, maf=0.05, input=TRUE, outfile=c("012","-101","structure")) {

  if (call.rate < 0 | call.rate > 1 | maf < 0 | maf > 1)
    stop("Treshold for call rate and maf must be between 0 and 1")
  
  if(missing(outfile))
    outfile = "012"

	if(!is.matrix(data))
    stop("Data must be matrix class")
  
  match.arg(frame)
  match.arg(outfile)
  
  if(isTRUE(base)){
    if (frame=="long"){
      if(ncol(data)>4)
        stop("For format long, the object must have four columns")
    
    bs <- unique(na.omit(as.vector(data[, 3:4])))
    if(!any(all(bs %in% c("A","C", "G", "T")) | all(bs %in% c("A", "B"))))
      stop("SNPs must be coded as nitrogenous bases (ACGT) or as A and B")
    
    sample.id <- sort(unique(data[,1L]))
    snp.name <- sort(unique(data[,2L]))
    
    col2row <- function(x, data){
      curId <- data[,1L] %in% x
      curSnp <- ifelse(is.na(data[curId, 3L]) | is.na(data[curId, 4L]), NA, 
                       paste(data[curId, 3L], data[curId, 4L], sep = ""))
      curPos <- match(snp.name, data[curId, 2L])
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
    bs <- unique(unlist(strsplit(unique(data[!is.na(data)]), "")))
    if(!any(all(bs %in% c("A","C", "G", "T")) | all(bs %in% c("A", "B"))))
      stop("SNPs must be coded as nitrogenous bases (ACGT) or as AB")
  }
  
   count_allele <- function(m){
    #' @importFrom stringr str_count
    A <- matrix(str_count(m, "A"), ncol = ncol(m), byrow = FALSE)
    C <- matrix(str_count(m, "C"), ncol = ncol(m), byrow = FALSE)
    G <- matrix(str_count(m, "G"), ncol = ncol(m), byrow = FALSE)
    C[, colSums(A, na.rm = TRUE)!=0] <- 0
    G[, colSums(A, na.rm = TRUE)!=0 | colSums(C, na.rm = TRUE)!=0] <- 0
    res <- A + C + G
    if (any(colSums(res, na.rm=TRUE) == 0))
      res[,colSums(res, na.rm=TRUE) == 0] <- 2
    res[is.na(m)] <- NA
    rownames(res) <- rownames(m)
    colnames(res) <- colnames(m)
    return(res)
  }
    
  m <- count_allele(data)
  } else{
    if(frame=="long")
      stop("format long only works with nitrogenous bases. Check base argument")
    m <- data
  }

  miss.freq <- rowSums(is.na(m))/ncol(m)
  
  if (sweep.sample < 0 | sweep.sample > 1)
      stop("Treshold for sweep.clean must be between 0 and 1")
	id.rmv <- rownames(m)[miss.freq > sweep.sample]
    m <- m[miss.freq <= sweep.sample,]
    data <- data[miss.freq <= sweep.sample,]
    
  CR <- (colSums(!is.na(m)) - colSums(is.na(m)))/colSums(!is.na(m))
  
  p <- colSums(m, na.rm = TRUE)/(2*colSums(!is.na(m)))
  minor <- apply(cbind(p, 1-p), 1, min)
  minor[is.nan(minor)] <- 0
  
  snp.rmv <- vector("list", 2)
  snp.rmv[[1]] <- colnames(m)[CR < call.rate]
  snp.rmv[[2]] <- colnames(m)[minor < maf]
  
  position <- (CR >= call.rate) & (minor >= maf)
  if (sum(position)==0L)
     stop("All markers were removed. Try again with another treshold for CR and MAF")
    m <- m[, position]
    data <- data[, position]
    
  if (input==TRUE && any(!is.finite(CR[position])))
      stop("There are markers with all missing data. There is no way to do
           imputation. Try again using call rate other than zero")
  
  all.equal_ <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
  
  if (any(CR!=1L) && isTRUE(input))
  {
    if (any(all.equal_(miss.freq[miss.freq <= sweep.sample], 1L)))
       stop("There are samples with all missing data. There is no way to do
           imputation. Try again using sweep.sample other than zero")

    f <- rowSums(m!=1, na.rm = TRUE)/rowSums(!is.na(m))
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
    
    m <- input.fun(m=m, p=p[position], f=f)
  }
  
  if (outfile=="-101")
    m <- m - 1

  if(outfile=="structure"){
    m <- lapply(as.data.frame(data), function(x){
      curCol <- do.call(rbind, strsplit(as.character(x), split = ""))
      if(all(is.na(curCol))) {curCol <- cbind(curCol, curCol)}
      return(curCol)})
    m <- as.matrix(do.call(cbind, m))
    colnames(m) <- rep(colnames(data), each=2)
    
    m <- chartr("ACGT", "1234", m)
    m[is.na(m)] <- -9
  }
  
  report <- list(paste(length(snp.rmv[[2]]), "Markers removed by MAF =", maf, sep = " "),
                 snp.rmv[[2]],
                 paste(length(snp.rmv[[1]]), "Markers removed by Call Rate =", call.rate, sep=" "),
                 snp.rmv[[1]],
                 paste(length(id.rmv), "Samples removed by sweep.sample =", sweep.sample, sep = " "),
                 id.rmv,
                 paste(sum(is.na(data)), "markers were inputed = ", round((sum(is.na(data))/length(data))*100, 2), "%"))
  
  if(missing(hapmap)){
    storage.mode(m) <- "numeric"
    return(list(M.clean=m, report=report))
  } else{
    storage.mode(m)  <- "numeric"
    hap <- hapmap[hapmap[,1L] %in% colnames(m),]
    hap <- hap[order(hap[,2L], hap[,3L], na.last = TRUE, decreasing = F),]
    colnames(hap) <- c("SNP","Chromosome","Position")
    return(list(M.clean=m, Hapmap=hap, report=report))
  }
}

