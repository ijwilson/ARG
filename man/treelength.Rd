\name{treelength}
\alias{treeheight}
\alias{treelength}
\title{Calculate the tree lengths on an ARG}
\description{
  Calculate the lengths of trees under each site of an ancestral
  recombination graph
}
\usage{
treelength(graph)
treeheight(graph)
}
\arguments{
  \item{graph}{An ancestral recombination graph}
}
\value{
  Returns an array of times (one for each site).
}
\seealso{
  \code{\link{simARG}}
}
\examples{
a <- simARG(ss=100,sites=1000,rec=0.01)
par(mfrow=c(2,1))
plot(treeheight(a),type="l");
plot(treelength(a),type="l");
}
\author{Ian Wilson <I.J.Wilson@ncl.ac.uk>}
\keyword{utilities}
