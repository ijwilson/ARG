\name{mutate}
\alias{mutate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mutate an ARG}
\description{
  Mutates an ancestral recombination graph.
}
\usage{
mutate(graph, var, ascert = "panel(3)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{ The ARG to be mutated }
  \item{var}{ The number of variable sites (must be less than a tenth
  of the number of sites in the ARG)}
\item{ascert}{ The type of ascertainment.  A string that describes
  the ascertainment model.  Can be panel, doublepanel, random, older or
  bylength.  See details for more information about panel specification
}
}
\details{
 The \code{ascert} types are quoted strings that give the ascertainment method
 to use:
\describe{
 \item{panel(a)}{ Uses a panel of \code{a} diploids (the first 2*a)
 chromosomes in the sample.  For a site to be accepted it must be
 variable in the panel. }
 \item{bylength}{ Chooses sites for mutations proportional to their
 length.}
 \item{older(t)}{ Chooses mutations by length, but only keeps those that
 are older than time \code{t}.}
 \item{random}{Chooses sites at random.}
 \item{doublepanel(a,b,c)}{Uses two panels of size \code{a} and
 \code{c} diploids where the first panel is from 1 to 2\code{a} and the
 second panel is from 2\code{b} to 2\code{b}+2\code{c}.  SNPs must be
 variable in both panels to be accepted.}
}
}
\value{
  A list that has class \code{mutatedARG}.  This inherits other
  information from the ARG generated in \code{\link{simARG}}.
  \item{haplotype }{A matrix of haplotypes}
  \item{position }{The positions of the variable sites}
}
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
\examples{
a <- simARG(200,1000,.01)
am <- mutate(a,50,"panel(4)")
##check that the panel is correct
## there should be no monomorphic sites in the first 8
cs <- colSums(am$haplotype[1:8,])
sum(cs==0|cs==8)
image(t(am$haplotype))
## now with a different ascertainment scheme
## must be polymorphic in the first 6 samples and in samples 41 to 46
bb <- mutate(a,50,"doublepanel(3,20,3)")
cs1 <- colSums(bb$haplotype[1:6,])
cs2 <- colSums(bb$haplotype[41:46,])
sum(cs1==0|cs1==6)
sum(cs2==0|cs2==6)
image(t(bb$haplotype))
}
\keyword{datagen}

