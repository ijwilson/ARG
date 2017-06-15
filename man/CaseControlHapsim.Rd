\name{CaseControlHapsim}
\alias{CaseControlHapsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CaseControlHapsim(n = 60, penetrance = c(0, 1), sites = 2010, var = 200, growthmodel = "constant", freq = 0.1, rho = 0, panel, fdiff = 0.02, loud = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ The number of case and control haplotypes to simulate}
  \item{penetrance}{ The probability of the disease given the wild type
    and disease allele}
  \item{sites}{ The number of sites to simulate }
  \item{var}{ The number of variable sites to simulate}
  \item{growthmodel}{What is the growth model }
  \item{freq}{The frequency of the disease allele }
  \item{rho}{Recombination rate between sites }
  \item{panel}{ The size of the panel which determines which SNPs are seen }
  \item{fdiff}{ The absolute value of the difference between the
    required disease frequency and that accepted.  }
  \item{loud}{ Do you want plenty of messages? }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{locations }{1s give cases and 2 controls}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Ian Wilson <I.J.Wilson@ncl.ac.uk>}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{simARG}},\code{\link{CaseControlDipsim}}}
\examples{
m <- CaseControlHapsim(penetrance=c(0.2,0.7))
library(genomic)
plot(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

