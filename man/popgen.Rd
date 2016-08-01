\name{popgen}
\alias{popgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
popgen(Z, subgroups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{subgroups}{
%%     ~~Describe \code{subgroups} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z, subgroups) 
{
    Z <- as.matrix(Z)
    if (missing(subgroups)) {
        subgroups <- rep(1, nrow(Z))
    }
    subgroups <- as.factor(subgroups)
    X <- as.data.frame(cbind(subgroups, Z))
    g.of.p <- function(X) {
        M <- X[, -1]
        G <- X[, 1]
        m <- ncol(M)
        g <- length(G)
        f <- apply(M, 2, sum)/(2 * g)
        fs <- cbind(f, 1 - f)
        MAF <- apply(fs, 1, min)
        p <- f
        q <- 1 - f
        Hesp <- 2 * p * q
        H <- function(M) {
            sum(M[M == 1])
        }
        Hobs <- apply(M, 2, H)/(2 * g)
        Dg <- 1 - p^2 - q^2
        PIC <- 1 - (p^2 + q^2) - (2 * p^2 * q^2)
        markers <- round(data.frame(p, q, MAF, He = Hesp, Ho = Hobs, 
            DG = Dg, PIC), 2)
        Hg.obs <- apply(M, 1, H)/(2 * m)
        Fi <- 1 - Hg.obs/mean(Hesp)
        Si <- (2 * Fi)/(1 + Fi)
        genotypes <- round(data.frame(Ho = Hg.obs, Fi, Si), 2)
        DG.pop <- c(mean = mean(Dg), lower = range(Dg)[1], upper = range(Dg)[2])
        PIC.pop <- c(mean = mean(PIC), lower = range(PIC)[1], 
            upper = range(PIC)[2])
        MAF.pop <- c(mean = mean(MAF), lower = range(MAF)[1], 
            upper = range(MAF)[2])
        Hg.obs.pop <- c(mean = mean(Hg.obs), lower = range(Hg.obs)[1], 
            upper = range(Hg.obs)[2])
        Fi.pop <- c(mean = mean(Fi), lower = range(Fi)[1], upper = range(Fi)[2])
        Si.pop <- c(mean = mean(Si), lower = range(Si)[1], upper = range(Si)[2])
        population <- t(round(data.frame(DG = DG.pop, PIC = PIC.pop, 
            MAF = MAF.pop, Ho = Hg.obs.pop, F = Fi.pop, S = Si.pop), 
            2))
        Ne <- 1/(2 * mean(Fi)) * g
        Va <- sum(2 * p * q)
        Vd <- sum((2 * p * q)^2)
        NG <- length(unique(G))
        variance <- t(round(data.frame(Ne, Va, Vd, `number of groups` = NG, 
            `number of genotypes` = g, `number of markers` = m), 
            2))
        colnames(variance) <- ("estimate")
        average <- list(markers = markers, genotypes = genotypes, 
            population = population, variability = variance)
        return(average)
    }
    general <- g.of.p(X)
    if (length(unique(subgroups)) != 1) {
        groups <- split(X, X$subgroups)
        bygroup <- lapply(groups, g.of.p)
        pbyg <- matrix(0, ncol(Z), length(unique(subgroups)), 
            dimnames = list(colnames(Z), c(1:length(unique(subgroups)))))
        for (i in 1:length(unique(X$subgroups))) {
            pbyg[, i] <- bygroup[[i]]$markers$p
        }
        fallelle <- list()
        exclusive <- list()
        for (i in 1:length(unique(X$subgroups))) {
            fixed <- pbyg[, i] == 1 | pbyg[, i] == 0
            present <- pbyg[, i] > 0
            exclusive[[i]] <- colnames(Z)[c(which(pbyg[, i] > 
                0 & apply(as.matrix(pbyg[, -i] == 0), 1, function(x) all(x))), 
                which(pbyg[, i] < 1 & apply(as.matrix(pbyg[, 
                  -i] == 1), 1, function(x) all(x))))]
            fallelle[[i]] <- colnames(Z)[fixed == TRUE]
            if (length(exclusive[[i]]) == 0) {
                exclusive[[i]] <- c("there are no exclusive alleles for this group")
            }
            if (length(fallelle[[i]]) == 0) {
                exclusive[[i]] <- c("there are no fixed alleles for this group")
            }
        }
        for (i in 1:length(unique(X$subgroups))) {
            fixed <- pbyg[, i] == 1
        }
        return(list(general = general, bygroup = bygroup, exclusive_alleles = exclusive, 
            fixed_alleles = fallelle))
    }
    else {
        bygroup <- c("there are no subgroups")
        return <- list(general = general, bygroup = bygroup)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
