\name{TMRCA}
\alias{TMRCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get TMRCAs from an ARG }
\description{
  Get the time to the most recent common ancestor for a pair of individuals.
}
\usage{
TMRCA(graph, position, samps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{ The ARG }
  \item{position}{ The \code{position} to look at. }
  \item{samps}{ The samples number of the pair that you want to find the TMRCA. }
}

\value{
The TMRCA at a single position.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Ian Wilson \email{Ian.Wilson@ncl.ac.uk} }

\seealso{ \code{\link{treeheight}}, \code{\link{treelength}} }
\examples{
a = simARG(100,100,0.1)
tmrca12 <- sapply(1:100,TMRCA,graph=a,samps=c(1,2))
#TMRCA(a,1,c(1,2))
th <- treeheight(a)
plot(th,type="l",xlab="Position",ylab="height",ylim=c(0,max(th)))
lines(tmrca12,col="red")
}

\keyword{ datagen }

