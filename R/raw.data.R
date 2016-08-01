raw.data <-
function(data, frame=c("table","matrix"), hapmap, sweep.sample= 0, call.rate=0.95, maf=0.05, input=TRUE, outfile=c("012","-101","structure")) {
  
  if (call.rate < 0 | call.rate > 1 | maf < 0 | maf > 1 | sweep.sample < 0 | sweep.sample > 1)
    stop("Treshold for call rate, maf and sweep.clean must be between 0 and 1")
  
  if(missing(outfile)) {outfile = "012"}
  
  if (frame=="table"){
    if(ncol(data)>4)
      stop("For format table, the object must have four columns")
    
    sample.id <- unique(data[,1])
    snp.name <- unique(data[,2])
    mbase <- sapply(snp.name,
               function(x){curRows <- data[,2] %in% x
               if(all(is.na(data[curRows, c(3,4)])))
               {curSnp <- matrix(NA, nrow(data[curRows, c(3,4)]), 1)
               }else{
                 curSnp <- ifelse(is.na(data[curRows, 3]) | is.na(data[curRows, 4]), NA, 
                                  paste(data[curRows, 3], data[curRows, 4], sep = ""))}
               return(curSnp)})
    colnames(mbase) <- snp.name
    rownames(mbase) <- sample.id
    data <- mbase}
  
  if (frame=="matrix"){
   data <- data }
  
  m <- sapply(as.data.frame(data), function(x){
    if (all(is.na(x))){
      return(x)
      }else{
        snp.col <- do.call(rbind, strsplit(as.character(x), split = ""))
        ref.allel <- names(which.max(table(unlist(snp.col))))
        count <- rowSums(snp.col == ref.allel)
        return(count)}
    })
    
  miss.freq <- rowSums(is.na(m))/ncol(m)
      
    if(sweep.sample==0){
      m1 <- m
    }else{
      m1 <- m[miss.freq <= sweep.sample,]
      data <- data[miss.freq <= sweep.sample,]
      }
  
    CR <- (colSums(!is.na(m1)) - colSums(is.na(m1)))/colSums(!is.na(m1))
    CR[!is.finite(CR)] <- 0
    
    
    p <- (2*colSums(m1==2, na.rm=TRUE) +
            colSums(m1==1, na.rm=TRUE))/(2*colSums(!is.na(m1)))
    minor <- apply(cbind(p,1-p), 1, min)
    minor[is.nan(minor)] <- 0
    
    if (call.rate == 0 & maf == 0)
    {
      m2 <- m1
      }else
        {
        position <- which(CR >= call.rate & minor >= maf)
        if (length(position)==0L){
          return(message("No marker selected. Try again with another treshold for CR and MAF"))}
        m2 <- m1[,position]
        data <- data[,position]
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
    
    if (outfile=="-101"){
      m4 <- apply(m3, 2, function(x){ 
        chan <- ifelse(is.na(x),NA, x-1)
        return(chan)})
      }
  if(outfile=="structure"){
    m4 <- lapply(as.data.frame(data), function(x){
      curCol <- do.call(rbind, strsplit(as.character(x), split = ""))
      if(all(is.na(curCol))) {curCol <- cbind(curCol, curCol)}
      return(curCol)})
    m4 <- as.matrix(do.call(cbind, m4))
    colnames(m4) <- rep(colnames(data), each=2)
    m4 <- gsub("A", 1, m4)
    m4 <- gsub("C", 2, m4)
    m4 <- gsub("G", 3, m4)
    m4 <- gsub("T", 4, m4)
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
