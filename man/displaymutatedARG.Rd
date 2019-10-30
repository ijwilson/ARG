\name{displaymutatedARG}
\alias{displaymutatedARGsmall}

\title{ Function to produce an plot of an ARG}
\description{
  Produces an interactive plot of a mutated Ancestral
  Recombination Graph.  This plot can be interactive.
}
\usage{
displaymutatedARG(b, interactive = TRUE)
}
\arguments{
  \item{b}{ an object of class mutatedARG }
  \item{interactive}{ Use an interactive plot?}
}
%'\details{
% # ~~ If necessary, more details than the description above ~~
%}
\value{
  Does not return a value.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ simARG,mutate}
\examples{
 bb <- mutate(simARG(40,1000,0.1),100)
 displaymutatedARG(bb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }

