\name{tree}
\alias{tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract a tree from an Ancestral Recombination Graph }
\description{
  This function extracts the tree at a particular position from
  an ancestral recombination graph (ARG).   
}
\usage{
tree(graph, position, correct=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{ The ancestral recombination graph considered }
  \item{position}{ The position of the tree to extract }
  \item{correct}{Correct the node labels to read 1 to n?}
}
\details{
  Returns an object of class \code{phylo} which can be plotted using the
  \code{ape} class.
}
\value{
  An object of class \code{phylo}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Ian Wilson <I.J.Wilson@ncl.ac.uk>}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
a <- simARG(20,9,2)
h <- max(treeheight(a))
par(mfrow=c(3,3),mar=c(1,1,1,1))
lapply(1:9,function(x) plot(tree(a,x)))
}
\keyword{ datagen }
