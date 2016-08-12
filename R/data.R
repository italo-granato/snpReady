#' maize hybrids
#' 
#' 50 hybrids of maize genotyped with 500 SNP markers
#'
#' @rdname maize.hyb 
#' @name maize.hyb
#' @docType data
#' @format A matrix with 50 rows and 500 columns, where hybrids are in rows and SNP markers in columns
#' @keywords datasets; hybrids
#' @examples
#' #' data(maize.hyb)
#' 
#' 
NULL

#' maize lines
#' 
#' A raw dataset of maize lines genotyped with 768 markers
#' 
#' @rdname maize.line
#' @name maize.line
#' @docType data
#' @format A data frame with 70656 observations on the following 4 variables.
#' \itemize{ 
#' \item{sample: }{identification of samples (name of individuals)}
#' \item{marker: }{identification of SNP markers}
#' \item{allele.1: }{Allele 1}
#' \item{allele.2: }{Allele 2}}
#' @source Lines genotyped from allogamous breeding laboratory - ESALQ/USP \url{http://www.genetica.esalq.usp.br/alogamas/index2.html}
#' @keywords datasets; lines
#' @examples
#' 
#' data(maize.line)
#' ## str(maize.line)
#' 
NULL
