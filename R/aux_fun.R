install.dep <- function(package='impute'){
  ## if impute isn't installed the function will try to install it
  if(!package %in% installed.packages()[,'Package']){
    
    ## if the bioconductor manager package is not installed try to install it
    if (!requireNamespace("BiocManager", quietly = TRUE))
    {
      ## messaging that this is been done 
      cat("The package 'impute', used as dependency, is being installed from Bioconductor", '\n')
      ## try to install the BiocManager
      install.packages("BiocManager", quiet = T)
      ## in case the bioc manager was successfuly installed
      if(requireNamespace("BiocManager", quietly = TRUE)){
        BiocManager::install(package, quietly=TRUE, update = FALSE) ## installing inpute
        ## remove.packages('BiocManager') ## removing the biocmanager
        
      }else stop("The package 'impute' could not be installed. The package need to be installed manually. A possible solution can be found in https://github.com/italo-granato/snpReady#installation")
      
    }else{
      BiocManager::install(package, quietly=TRUE, update = FALSE)  
    }
  }
}