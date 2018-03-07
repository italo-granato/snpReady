raw.data <- function(data, frame = c("long","wide"), hapmap = NULL, base = TRUE, sweep.sample= 1,
                     call.rate=0.95, maf=0.05, imput=TRUE, imput.type = c("wright", "mean","knni"),
                     outfile=c("012","-101","structure"), plot = FALSE)
{

  if (call.rate < 0 | call.rate > 1 | maf < 0 | maf > 1)
    stop("Treshold for call rate and maf must be between 0 and 1")
  
 if (sweep.sample < 0 | sweep.sample > 1)
      stop("Treshold for sweep.sample must be between 0 and 1")										  
  if(missing(outfile))
    outfile = "012"

	if(!is.matrix(data))
    stop("Data must be in matrix class")
  
  match.arg(frame)
  match.arg(imput.type)
  
  if(base){
    if (frame == "long"){
      if(ncol(data)>4)
        stop("For format long, the object must have four columns")
    
    bs <- unique(na.omit(as.vector(data[, 3:4])))
    if(!any(all(bs %in% c("A","C", "G", "T")) | all(bs %in% c("A", "B"))))
      stop("SNPs must be coded as nitrogenous bases (ACGT) or as A and B")
    
    if(!all(bs %in% c("A","C", "G", "T")) & outfile == "structure")
      stop("For outfile 'structure', SNPs must be coded as nitrogenous bases (ACGT)")
    
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
    
    if(!all(bs %in% c("A","C", "G", "T")) & outfile == "structure")
      stop("For outfile 'structure', SNPs must be coded as nitrogenous bases (ACGT)")
  }
  
	whichBase <- all(bs %in% c("A","C", "G", "T"))											  
   count_allele <- function(m, nitrBase = T){
    #' @importFrom stringr str_count
	if(nitrBase){		   
    A <- matrix(str_count(m, "A"), ncol = ncol(m), byrow = FALSE)
    C <- matrix(str_count(m, "C"), ncol = ncol(m), byrow = FALSE)
    G <- matrix(str_count(m, "G"), ncol = ncol(m), byrow = FALSE)
    C[, colSums(A, na.rm = TRUE)!=0] <- 0
    G[, colSums(A, na.rm = TRUE)!=0 | colSums(C, na.rm = TRUE)!=0] <- 0
    res <- A + C + G
	}else{
        A <- matrix(str_count(m, "A"), ncol = ncol(m), byrow = FALSE)
        res <- A
    }	
    if (any(colSums(res, na.rm=TRUE) == 0))
      res[,colSums(res, na.rm=TRUE) == 0] <- 2
    res[is.na(m)] <- NA
    rownames(res) <- rownames(m)
    colnames(res) <- colnames(m)
    return(res)
  }
    
  m <- count_allele(data, nitrBase = whichBase)
  }else{
    if(frame == "long")
      stop("format long only accepts nitrogenous bases. Check base argument")
    
    if(outfile == "structure")
      stop("output for 'structure' only accepts nitrogenous bases. Check base argument")
    
    m <- data
  }
  
  if(is.null(colnames(m)))
    stop("Marker names are missing")
	
  if(!is.null(hapmap)){
    if(ncol(m) != nrow(hapmap))
      stop("The number of markers differ in the hapmap and the data")
    
    hapmap <- hapmap[order(hapmap[,2L], hapmap[,3L], na.last = TRUE, decreasing = F), ]
    
    tmp <- match(hapmap[,1L], colnames(m)) #match between marker names
    
    if(any(is.na(tmp)))
      stop("Marker names in data does not match the marker names in hapmap")
											
    
    m <- m[,tmp]
    
    tmpd <- match(hapmap[,1L], colnames(data)) #match between marker names
    data <- data[,tmpd]
  }
  
  idName <- rownames(m)
  markerName <- colnames(m)
	
  missID <- rowMeans(is.na(m))
  posSS <- missID <= sweep.sample #selected by ID call rate
  
  #id.rmv <- rownames(m)[missID > sweep.sample]
  m <- m[posSS,]
  data <- data[posSS,]
  
  CR <- 1 - colMeans(is.na(m))
  poscr <- CR >= call.rate #selected by CR							  
  p <- colMeans(m, na.rm = T)/2
  minor <- apply(cbind(p, 1-p), 1, min)
  minor[is.nan(minor)] <- 0
  posmaf <- minor >= maf #selected by maf	
  #snp.rmv <- vector("list", 2)
  #snp.rmv[[1]] <- colnames(m)[CR < call.rate]
  #snp.rmv[[2]] <- colnames(m)[minor < maf]
  #position <- (CR >= call.rate) & (minor >= maf)
  if (sum(poscr & posmaf)==0L)
     stop("All markers were removed. Try again with another treshold for CR and MAF")
   
  m <- m[, poscr & posmaf]
  data <- data[, poscr & posmaf]
  
  if (imput){
    
  if (any(CR[poscr & posmaf] == 0))
      stop("There are markers with all missing data. There is no way to
           impute. Try again using another call rate treshold")
  
   if(imput.type == "wright"){
      
      if (any(missID[posSS] == 1L))
        stop("There are samples with all missing data. There is no way to do
             imputation. Try again using another sweep.sample treshold")
      
							   
      f <- rowMeans(m != 1, na.rm = TRUE)
      f[is.nan(f)] <- 1
      m <- input.fun(m=m, p=p[poscr & posmaf], f=f)
    } 
    if(imput.type == "mean"){
      tmp <- which(is.na(m), arr.ind = TRUE)
      m[tmp] <- colMeans(m, na.rm = T)[tmp[,2]]
    }
    if(imput.type == "knni"){
      m <- t(impute.knn(t(m), 2)$data)
    }
  }
  
  switch(outfile,
         "-101" = {
           m[m == 0] <- -1
           m[m == 1] <- 0
           m[m == 2] <- 1
           },
         "structure" = {
           tmp <- lapply(as.data.frame(data), function(x){
             curCol <- strsplit(as.character(x), split = "")
             curCol[ which(is.na(curCol))] <- list(rep(NA, 2))
             res <- unlist(curCol)
             return(res)} )
           m <- as.matrix(do.call(cbind, tmp))
           colnames(m) <- colnames(data)
           rownames(m) <- rep(rownames(data), each=2)
           m <- chartr("ACGT", "1234", m)
           m[is.na(m)] <- -9
         },
         "012" = {
           m <- m
         },
         {
         stop("Output selected is not available")
           }
         )
  
  report <- list(maf = list(r = paste(sum(!posmaf), "Markers removed by MAF =", maf, sep = " "),
                            whichID = markerName[!posmaf]),
                 cr = list(r = paste(sum(!poscr), "Markers removed by Call Rate =", call.rate, sep=" "),
                           whichID = markerName[!poscr]),
                 sweep = list(r = paste(sum(!posSS), "Samples removed by sweep.sample =", sweep.sample, sep = " "),
                              whichID = idName[!posSS]),
                 imput = ifelse(imput, paste(sum(is.na(data)), "data points were inputed = ", round((sum(is.na(data))/length(data))*100, 2), "%"),
                                "No data point was imputed"))
  
  for(i in 1:3){
    if(length(report[[i]]$whichID) == 0)
      report[[i]]$whichID <- NULL
  }
  
  storage.mode(m) <- "numeric"
  if(is.null(hapmap)){
    return(list(M.clean = m, report = report))
  } else{
    
    if(plot){
      nimp <- aggregate(colSums(is.na(data)), by = list(hapmap[poscr & posmaf, 2]), FUN = sum)$x
      
      bychrom <- t(cbind("cr" = table(hapmap[!poscr | (!poscr & !posmaf) , 2])/table(hapmap[,2]), 
                         "maf" = table(hapmap[!posmaf & poscr, 2])/table(hapmap[,2]),
                         "n.imput" = nimp/(table(hapmap[, 2])*nrow(data))) )
      barplot.rd(bychrom)
    }
    
    hap <- hapmap[poscr & posmaf, ]															  
    colnames(hap) <- c("rs","chrom","pos")	  
    return(list(M.clean = m, Hapmap = hap, report = report))
  }
}

barplot.rd <- function(data, col = gray.colors(3), plotName = "QCreport"){
  pdf = pdf(paste(plotName,".pdf", sep = ""), width = 10, height = 4)
  
  barplot(data, beside = TRUE, horiz = FALSE, ylab = "frequency", axes = F,
          xlab = "chromosome", ylim = c(0, max(data) ), 
          main = "Quality control", col = col)
  axis(side = 2, las = 2)
  par(new = T, mar = c(2,2,2,0))
  plot(1, type = "n", xlab = "", ylab = "", bty="n", axes = F)
  legend(x = "topright", fill = col, legend = c("Maf", "CR", "Imput"), 
         bty = "n")
  dev.off()
}

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
