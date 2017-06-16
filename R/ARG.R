

"simARG" <- function(ss,sites,rec,growthmodel="constant"
                     ,migmatrix,poptree)
{
  if (length(rec)==1) recomb <- rep(rec,sites-1)
  
  if (missing(migmatrix)&&missing(poptree)) {
    graph <- simplesim(ss, sites, rec, growthmodel) 
    r <- list(ss=ss,sites=sites,rec=rec,graph=graph)
    class(r) <- "ARG"
    return(r)
  } else if (missing(poptree)) {
    if (length(ss)==1) {
      warning("using same sample size for each population")
      ss <- rep(ss,npops(migmatrix))
    }
    if (length(ss)!=npops(migmatrix)) 
      stop("Length of ss not the same as the number of populations")
    
    loc <- rep(1:length(ss), rep(ss, length(ss)))
    graph <- mmsim(loc, sites, recomb, growthmodel, migmatrix)

    r <- list(ss=ss, sites=sites, rec=recomb
              ,graph=graph, location=loc
              ,growthmodel=growthmodel
              ,migmatrix=migmatrix)
    class(r) <- c("structuredARG","ARG")
    return(r)
  } else if (missing(migmatrix)) {
    stop("poptree not defined yet")
  } else {
    stop("poptree not defined yet")
  }
}