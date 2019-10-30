\name{PopCaseControlHapsim}
\alias{PopCaseControlDipsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Simulate a Case and Control sample }
\description{
 Simulate a Case and Control sample with an additional random sampled population.
}
\usage{
PopCaseControlHapsim(n = c(100, 100), popsize = 100, sites = 1000, penetrance = c(0, 1), growthmodel = "constant", thetapars = c(1, 1), diseasefreq = 0.1, fdiff = 0.02, rho = 0, loud = FALSE)
PopCaseControlDipsim(n=c(100,100),popsize=100,sites=1000,penetrance=c(0.0,0.1,1.0),growthmodel="constant"
           ,thetapars=c(1,1),diseasefreq=0.1,fdiff=0.02,rho=0.0,loud=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ A vector of length two containing the number of cases and controls }
  \item{popsize}{ The size of the population sample }
  \item{sites}{The number of sites to simulate }
  \item{penetrance}{ A vector of length two (for Hapsim) or 3 (for
    Dipsim) giving the probability of
    the disease given the wild type allele (0) and the disease allele (1) }
  \item{growthmodel}{ The growth model for the population }
  \item{thetapars}{The parameters for the gamma distribution that gives
    the mutation rate --separate mutation rate for each site }
  \item{diseasefreq}{ The frequency of the disease allele}
  \item{fdiff}{ maximum allowable difference between the frequency of
    the disease allele ( a random quantity) and \code{diseasefreq}}
  \item{rho}{ The recombination rate between sites (either a scalar or a
    vector of length sites-1)}
  \item{loud}{ make the function talkative?}
}
\details{
  Simulate a Case-Control haploid or diploid sample with an addition
  sample from
  the population(that may include cases and control).  This function is
  used to try to mirror the SNP discovery/Sample Case-Controls for the
  discovered SNPs.  This particular function uses a different mutation
  model than the usual UEP SNP mutation model, rather there is a prior
  on the mutation rate theta, that varies across sites.
} \value{
 A List with elements:
  \item{haplotype}{A matrix of SNPs}
  \item{position }{The positions of the variable SNPs}
  \item{diseaseSNP}{The disease marker}
  \item{location}{The status of each sample, Case, Control or
    Population}
  \item{diseasemuts}{The number of mutations at the disease site}
  \item{muts}{The number of mutations at each of the sites}
  
}
\references{ \url{http://www.staff.ncl.ac.uk/i.j.wilson/}}
\author{ Ian Wilson \email{I.J.Wilson@ncl.ac.uk}}
\section{warning}{
  Each call of this function creates a data structure within the
  underlying C++ library.  This is not currently destroyed when it goes
  out of scope -- it shall be rewritte using the \code{.Call} interface
  in the future.  Hence one should explicity delete the underlying
  structure using \code{\link{ARG.remove}}.  Note also that other
  functions such as \code{\link{prune}} can affect the underlying
  C++ structure so that it may no longer match the returned
  object.  
}
\seealso{  \code{\link{ARG.remove}}  }
\examples{
a <- PopCaseControlHapsim(n=c(100, 100),popsize=12,sites=200, penetrance
= c(0.1,0.6),thetapars = c(1, 1),diseasefreq = 0.1, fdiff = 0.02, rho =
0.15, loud =TRUE)
image(t(a$haplotype))
## Now use the population tree to get markers
pop  <- a$haplotype[201:212,]
u <- colSums(pop)>0&colSums(pop)<12
## test markers against case control status
dd <- a$haplotype[1:200,u]
pos <- a$position[u]
ub <- colSums(dd)>0 & colSums(dd)<200
dd <- dd[,ub]
cc <- a$location[1:200]

pv <- apply(dd,2,function(x) chisq.test(x,cc,simulate=TRUE)$p.value)
plot(pos[ub],pv,xlab="position",ylab="p-value",axes=FALSE,log="y")
axis(1)
axis(2)

arrows(a$diseasepos,0.01,a$diseasepos,0.0005,col="blue")
text(a$diseasepos,0.005,"Position of Disease SNP",col="blue",pos=4)

a <- PopCaseControlDipsim(rho=0.1,loud=TRUE,penetrance=c(0.2,0.4,0.3))
 pop  <- a$haplotype[401:412,]
 u <- colSums(pop)>0&colSums(pop)<12
dd <- a$haplotype[1:400,u]
pos <- a$position[u]
 ub <- colSums(dd)>0 & colSums(dd)<400
 dd <- dd[,ub]
cc <- a$location[1:400]
 pv <- apply(dd,2,function(x) chisq.test(x,cc,simulate=TRUE)$p.value)
plot(pos[ub],pv,xlab="position",ylab="p-value",axes=FALSE,log="y")
axis(1)
axis(2)

arrows(a$diseasepos,0.01,a$diseasepos,0.0005,col="blue")
text(a$diseasepos,0.005,"Position of Disease SNP",col="blue",pos=4)
}
\keyword{ datagen }

