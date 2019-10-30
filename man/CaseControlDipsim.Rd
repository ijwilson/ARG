\name{CaseControlDipsim}
\alias{CaseControlDipsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Simulate Diploids from a Case Control model.
}
\usage{
CaseControlDipsim(n = 20, penetrance = c(0, 0.7, 1), sites = 2010, var = 200, growthmodel = "constant", freq = 0.1, rho = 0.1, panel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{The number of cases and controls }
  \item{penetrance}{ a vector which gives the penetrance of the disease
    for 0, 1 and 2 copies of the mutant allele.}
  \item{sites}{ The number of sites to simulate.  Recombinations can
    only occur between sites.}
  \item{var}{ The number of variable sites.  Should be less than 10\% of
    the sites}
  \item{growthmodel}{ ~~Describe \code{growthmodel} here~~ }
  \item{freq}{ ~~Describe \code{freq} here~~ }
  \item{rho}{ ~~Describe \code{rho} here~~ }
  \item{panel}{ ~~Describe \code{panel} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
\note{ ~~further notes~~ }

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n = 20, penetrance = c(0, 0.7, 1), sites = 2010, var = 200, 
    growthmodel = "constant", freq = 0.1, rho = 0.1, panel) 
{
    selectdip <- function(u, penetrance, n) {
        genotypes <- matrix(u, ncol = 2, byrow = T)
        n.geno <- nrow(genotypes)
        m <- rowSums(genotypes) + 1
        pu <- penetrance[m]
        sel <- matrix(runif(n.geno) < pu)
        if (sum(sel) < n) 
            stop("need larger penetrance or more sporadic cases to get ", 
                n)
        seli <- as.vector(t(matrix(sel, nrow = n.geno, ncol = 2)))
        list(case = which(seli)[1:n], control = which(!seli)[1:n])
    }
    fmn <- freq - 0.02
    f <- c((1 - fmn)^2, 2 * fmn * (1 - fmn), fmn^2) \%*\% penetrance
    ss <- 2 * round((n + 4 * sqrt(n))/f)
    cat("Choose a total sample size of ", ss/2, "individuals,\n", 
        "we expect to get", round(ss * f/2), "affected individuals\nsampling tree\n")
    a <- simARG(ss, sites, rho, growthmodel = growthmodel)
    if (missing(panel)) 
        ascert <- paste("panel(", as.integer(n/2), ")", sep = "")
    else ascert = paste("panel(", as.integer(panel), ")", sep = "")
    cat("mutating tree with ascertainment", ascert, "\n")
    b <- mutate(a, (var + 1), ascert = ascert)
    f1 <- apply(b$d, 2, mean)
    w1 <- which(abs(f1 - freq) < 0.02)
    if (length(w1) == 1) 
        stop("no sites with the correct frequency\n")
    wh <- w1[which(abs(w1 - var/2) == min(abs(w1 - var/2)))]
    if (length(wh) > 1) 
        wh <- sample(wh, 1)
    pos <- b$position[wh]
    u1 <- b$d[, wh]
    cc <- selectdip(u1, penetrance, n)
    b <- prune(b, c(cc$case, cc$control))
    b$SNP <- b$d[, wh]
    b$ss <- c(n, n)
    b$d <- b$d[, -wh]
    b$DiseasePos <- pos
    b$cc <- cc
    b$locations <- rep(1:2, c(n, n))
    b$positions <- b$positions[-wh]
    class(b) <- c("CCDip", "mutatedstructuredARG", class(b))
    b
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
