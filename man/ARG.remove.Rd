\name{ARG.remove}
\alias{ARG.remove}
\alias{prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to manipulate, extract and remove ARG data.}
\description{
  Functions to remove and manipulate ARG data.  
}
\usage{
ARG.remove(graph)
prune(graph,keep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{An object that inherits from the \code{ARG} class}
  \item{keep}{A vector that gives the indices of the end points to keep.}
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ \url{http://www.mas.ncl.ac.uk/~nijw}}
\author{Ian Wilson }
\seealso{ \code{\link{treelength}},\code{\link{simARG}} }
\examples{
## generate 9 trees
a <- replicate(9,simARG(100,1000,0.01),simplify=FALSE)
b <- lapply(a,function(x) prune(x,sample(100,10)))
par(mfrow=c(3,3),mar=rep(1,4))
lapply(b,function(x) plot(tree(x,1)))
}
\keyword{datagen }% at least one, from doc/KEYWORDS

