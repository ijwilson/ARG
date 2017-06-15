\name{mutateSTR}
\alias{mutateSTR}
\alias{mutatetheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mutate an ARG }
\description{
 Mutate an already existing ARG.
}
\usage{
mutateSTR(graph, theta, positions)
mutatetheta(graph, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{The Ancestral Recombination Graph to mutate }
  \item{theta}{ The scaled mutation rate }
  \item{positions}{ The positions to mutate (for mutateSTR) }
}
\details{
  Mutates an ARG.  For mutateSTR just a few positions are mutated, for
  mutatetheta
  all sites are mutated.  \code{mutatetheta} assumes a biallelic model
  with mutation rate of theta between alleles.  \code{mutateSTR} assumes
  a stepping stone mutation model suitable for modelling an SR.
}
\value{
  A list with components
  \item{ss }{The sample size}
  \item{sites }{The number of sites in the ARG}
  \item{rec}{The recombination rate between sites}
  \item{d}{A matrix giving the haplotypes (by row)}
  \item{positions}{The positions of the variable sites}
  \item{which}{The label of the underlying ARG}
  \item{var}{The number of variable sites}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Ian Wilson \email{Ian.Wilson@newcastle.ac.uk}}

\seealso{ \code{\link{mutate}}}
\examples{
aa <- simARG(ss=100,sites=200,rec=0.5,growthmodel="exponential(4)")
bb  <- mutatetheta(aa,rgamma(200,1,10))
image(t(bb$haplotype),axes=FALSE)
ARG.remove(aa)
#
aa <- simARG(ss=100,sites=200,rec=0.5,growthmodel="constant")
dd <- mutateSTR(aa,theta=rep(10,200),position=seq(1,200))

lengthvar <- apply(dd$haplotype,2,var)
het <- 1-unlist(apply(dd$haplotype,2,function(x) {t <- table(x);sum(t^2)/100^2}))
par(mfrow=c(2,2))
 plot(het~treelength(aa)[dd$position],ylab="Expected Heterozygosity",xlab="Length of Tree under Position")
 plot(lengthvar~treelength(aa)[dd$position],ylab="Variance of Allele Length",xlab="Length of Tree under Position")
 plot(het~treeheight(aa)[dd$position],ylab="Expected Heterozygosity",xlab="Height of Tree under Position")
 plot(lengthvar~treeheight(aa)[dd$position],ylab="Variance of Allele Length",xlab="Height of Tree under Position")
}

\keyword{ datagen }

