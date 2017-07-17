raw.data <- function(data, frame = c("long","wide"), hapmap = NULL, base = TRUE, sweep.sample= 1,
                     call.rate=0.95, maf=0.05, imput=TRUE, type = c("wright", "mean"),
                     outfile=c("012","-101","structure"))
{

  if (call.rate < 0 | call.rate > 1 | maf < 0 | maf > 1)
    stop("Treshold for call rate and maf must be between 0 and 1")
  
  if(missing(outfile))
    outfile = "012"

	if(!is.matrix(data))
    stop("Data must be in matrix class")
  
  match.arg(frame)
  match.arg(outfile)
  
  if(isTRUE(base)){
    if (frame == "long"){
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
        }
      else{
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
  }else{
    if(frame=="long")
      stop("format long only accepts nitrogenous bases. Check base argument")
    
    m <- data
  }
  
  if(is.null(colnames(m)))
    stop("Colnames is missing")
  
  miss.freq <- rowSums(is.na(m))/ncol(m)
  
  if (sweep.sample < 0 | sweep.sample > 1)
      stop("Treshold for sweep.clean must be between 0 and 1")
  
	id.rmv <- rownames(m)[miss.freq > sweep.sample]
    m <- m[miss.freq <= sweep.sample,]
    data <- data[miss.freq <= sweep.sample,]
    
  CR <- 1 - colSums(is.na(m))/nrow(m)
  
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
    
  if (isTRUE(imput) & any(!is.finite(CR[position])))
      stop("There are markers with all missing data. There is no way to do
           imputation. Try again using another call rate treshold")
  
  all.equal_ <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})
  
  samplefp <- function(p, f){
    samp <- sample(c(0,1,2), 1,
                   prob=c(((1-p)^2+((1-p)*p*f)),
                          (2*p*(1-p)-(2*p*(1-p)*f)),
                          (p^2+((1-p)*p*f))))
    return(as.integer(samp))
  }
  
  input.fun <- function(m, p, f){
    indicesM <- which(x = is.na(m), arr.ind = TRUE)
    m[indicesM] <- mapply(samplefp, p[indicesM[,2]], f[indicesM[,1]])
    return(m)
  }
  
  if (any(CR!=1L) & isTRUE(imput))
  {
    if (any(all.equal_(miss.freq[miss.freq <= sweep.sample], 1L)))
       stop("There are samples with all missing data. There is no way to do
           imputation. Try again using another sweep.sample treshold")
    
    if(type == "wright"){
      f <- rowSums(m!=1, na.rm = TRUE)/rowSums(!is.na(m))
      f[is.nan(f)] <- 1
      m <- input.fun(m=m, p=p[position], f=f)
    }
    if(type == "mean"){
      tmp <- which(is.na(m), arr.ind = TRUE)
      m[tmp] <- colMeans(m, na.rm = T)[tmp[,2]]
    }
  }
  
  if (outfile=="-101") m <- m - 1

  if(outfile=="structure"){
    tmp <- lapply(as.data.frame(data), function(x){
      curCol <- strsplit(as.character(x), split = "")
      tmp <- lapply(curCol, function(x) if(any(is.na(x))){rep(NA, 2)}else{x})
      curCol <- unlist(tmp)
      return(curCol)})
    m <- as.matrix(do.call(cbind, tmp))
    colnames(m) <- colnames(data)
    rownames(m) <- rep(rownames(data), each=2)
    m <- chartr("ACGT", "1234", m)
    m[is.na(m)] <- -9
  }
  
  report <- list(maf = list(r = paste(length(snp.rmv[[2]]), "Markers removed by MAF =", maf, sep = " "),
                            whichID = snp.rmv[[2]]),
                 cr = list(r = paste(length(snp.rmv[[1]]), "Markers removed by Call Rate =", call.rate, sep=" "),
                           whichID = snp.rmv[[1]]),
                 sweep = list(r = paste(length(id.rmv), "Samples removed by sweep.sample =", sweep.sample, sep = " "),
                              whichID = id.rmv),
                 imput = ifelse(isTRUE(imput), paste(sum(is.na(data)), "markers were inputed = ", round((sum(is.na(data))/length(data))*100, 2), "%"),
                                "No marker was imputed"))
  
  for(i in 1:3){
    if(length(report[[i]]$whichID) == 0)
      report[[i]]$whichID <- NULL
  }
  
  if(is.null(hapmap)){
    storage.mode(m) <- "numeric"
    return(list(M.clean = m, report = report))
  } else{
    storage.mode(m)  <- "numeric"
    hap <- hapmap[hapmap[,1L] %in% colnames(m),]
    hap <- hap[order(hap[,2L], hap[,3L], na.last = TRUE, decreasing = F),]
    colnames(hap) <- c("SNP","Chromosome","Position")
    m <- m[, match(hap[,1L], colnames(m))]
    return(list(M.clean = m, Hapmap = hap, report = report))
  }
}

