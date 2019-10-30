\name{simARG}
\alias{simARG}
\title{Simulate an ancestral recombination graph}
\description{
  Simulate an ancestral recombination graph.  The 
}
\usage{
simARG(ss,sites,rec,growthmodel="constant",migmatrix,poptree)
}
\arguments{
  \item{ss}{Sample size - either a vector or number}
  \item{sites}{Number of sites between which recombination can occur}
  \item{rec}{Recombination rate between sites, either a single number or
    an array of length sites-1}
  \item{growthmodel}{growth model, e.g. "constant" or "exponential(2)"}
  \item{migmatrix}{The migration matrix.  Either an actual matrix or
    a string such as "Island(4,1)" that indicates an island model for 4
    populations with a migration rate of 1 individual per generation}
  \item{poptree}{A population tree - a string in Newick format.  Assumes
  a single population if missing.}
}
\value{
  The data is not stored within R, rather it is allocated within a
  static data structure within the C code, and information about it is
  extracted using R functions.  The list contains information about
  where it is stored in R, so that multiple structures may be stored.
  A list is produced with
  the some or all of the components:
  \item{ss}{The sample size}
  \item{sites}{The number of sites sampled}
  \item{rec}{The recombination rate between sites}
  \item{wh}{A label that gives which ARG this is (for internal package
    use)}
  \item{growthmodel}{The growthmodel used for simulation.  Can be
    "constant", "exponential(alpha)",
    "piecewise(a1,t1,a2,t2,...,ak,tk,a)" where ai is the scaled
    population size from time t(i-1) to t(i) and a is the population
    size for the rest of time.}
  \item{location}{A vector of the locations of samples}
  \item{migmatrix}{The migration matrix used}
}

\section{Warning}{ The Ancestral Recombination Graph structures may take
  up a large amount  of memory with large recombination rates.  If you want to use many
  large structures it is best to use \code{\link{ARG.remove}} to
  remove unwanted structures.
  }
\references{
  The Newick tree format:  
  \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
}
\seealso{
  \code{\link{treelength}}
  \code{\link{treeheight}}
  \code{\link{ARG.remove}}
}
\examples{
a <- simARG(ss=100,sites=1000,rec=0.1)
b <- simARG(ss=100,sites=1000,rec=0.1,growthmodel="exponential(1)")
plot(treeheight(a),type="l",ylim=c(0,10));
lines(treeheight(b),col="blue");
}
\author{Ian Wilson <I.J.Wilson@ncl.ac.uk>}
\keyword{datagen}
