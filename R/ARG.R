

"npops" <- function(migmattext) {
   b <- unlist(strsplit(migmattext,"\\("))
   if (b[1]=="Island")
     return(as.integer(unlist(strsplit(b[2],","))[1]))
   if (b[1]=="General") {
     n <- length(unlist(strsplit(b[2],",")))
     rn <- floor(sqrt(n))
     if (n!=rn^2)
       stop("A general model should give a square matrix")
     return(rn);
   }
   stop(paste("'",b[1],"'"," not recognised as a migration type",sep=""))
}

makemigmatrix <- function(type="Island",mig=1,npops=3) {
  paste(type,"(",npops,",",mig,")")
}


"simARG" <- function(ss,sites,rec,growthmodel="constant"
                     ,migmatrix,poptree)
{
  if (length(rec)==1) recomb <- rep(rec,sites-1)
 
  if (missing(migmatrix)&&missing(poptree)) {
    wh <- .C("simplesim"
            ,ret=as.integer(ss)
            ,as.integer(sites)
            ,as.double(recomb)
            ,as.character(growthmodel)
            ,as.integer(1)
             ,PACKAGE="ARG"
            )[[5]]
    if (wh==-1) {
      warning("error, no ARG produced")
      return(NA)
    }
    r <- list(ss=ss,sites=sites,rec=recomb,which=wh)
    class(r) <- "ARG"
    return(r)
  } else if (missing(poptree)) {
    if (length(ss)==1) {
      warning("using same sample size for each population")
      ss <- rep(ss,npops(migmatrix))
    }
    if (length(ss)!=npops(migmatrix)) 
      stop("Length of ss not the same as the number of populations")
    
    loc <- rep(1:length(ss),ss)-1
    wh <- .C("mmsim"
       ,ret=as.integer(loc)
       ,as.integer(length(loc))
       ,as.integer(sites)
       ,as.double(recomb)
       ,as.character(growthmodel)
       ,as.character(migmatrix)
       ,as.integer(1)
             ,PACKAGE="ARG"
       )[[7]]
  if (wh==-1) {
      warning("error, no ARG produced")
      return(NA)
    }
    r <- list(ss=ss
              ,sites=sites
              ,rec=recomb
              ,which=wh
              ,location=loc+1
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
### function to extract the tree heights 
"treeheight" <- function(graph)
  {
    .C("treeheight",as.integer(graph$which),as.double(1:graph$sites)                                                     , PACKAGE="ARG")[[2]]
  }
##
##
##
### function to extract the tree lengths
"treelength" <- function(graph)
  {
    .C("treelength",as.integer(graph$which),as.double(1:graph$sites)
       ,PACKAGE="ARG")[[2]]
  }
##
##
##
"tree" <-  function(graph,position,correct=TRUE) {
  if (position<1||position>graph$sites) {
    stop("position should be between 1 and sites")
  }
  nchars=sum(graph$ss*64);
  txt <- .C("extracttree"
            ,as.integer(graph$which)
            ,as.integer(position-1)
            ,raw(nchars)
         #   ,character(nchars)
            ,as.integer(nchars)
            ,as.integer(correct)
            ,PACKAGE="ARG")[[3]]


   if (!correct)
    return(read.tree(text=rawToChar(txt)))
  t <- read.tree(text=rawToChar(txt))
  if (correct) {
    tips <- as.numeric(t$tip.label)+1
    t$tip.label <- tips
  }
  t
}
##
##
##
"TMRCA" <- function(graph,position,samps)
  {
    if (length(samps)!=2)
      stop("Expect two samples")
    if (max(samps)>sum(graph$ss)||min(samps)<1)
       stop("Samples should be between 1 and number of samples")
    if (samps[1]==samps[2]) {
      warning("samples should be different")
      return(0.0)
    }
    if (position<1||position>graph$sites) {
      stop("position should be between 1 and sites")
    }
    
    .C("TMRCA"
            ,as.integer(graph$which)
            ,as.integer(position-1)
            ,as.integer(samps-1)
            ,as.double(numeric(1))
            ,PACKAGE="ARG"
       )[[4]];
  }
##
##
##
"mutate" <- function(graph,var,ascert="panel(3)")
  {
    if (var>0.1*graph$sites||var<1) {
       stop("variable sites must be between 1 and 10% of available sites (but 1% better)")
    }
    ##  void mutate(int *which, int *var, int *data, int *u, char **ascertainment)
    res <- .C("mutate"
            ,as.integer(graph$which)
            ,as.integer(var)
            ,as.integer(numeric(var*sum(graph$ss)))
            ,as.integer(numeric(var))
            ,as.character(ascert)
            ,PACKAGE="ARG");
    
    r <- list(ss=graph$ss
              ,sites=graph$sites
              ,rec=graph$rec
              ,haplotype=matrix(res[[3]],ncol=var,byrow=TRUE)
              ,position=res[[4]]+1
              ,which=graph$which
              ,var=var);
    if (inherits(graph,"structuredARG")) {
      r$location <- graph$location
      class(r) <- c("mutatedstructuredARG","mutatedARG","structuredARG","ARG","GenomicHaplotype")
    } else {
      class(r) <- c("mutatedARG","ARG","GenomicHaplotype")
    }
    r
  }
## mutate 
##
##
##
"mutatetheta" <- function(graph,theta)
  {
    if (length(theta)!=graph$sites) {
       stop("length of theta must be the same as the number of variable sites")
    }
    samplesize <- sum(graph$ss)
    maxvar <- graph$sites

    r <- .C("mutatetheta"
            ,as.integer(graph$which)
            ,as.double(theta)
            ,data=as.integer(numeric(maxvar*samplesize))
            ,position=as.integer(numeric(maxvar))
            ,varsites=as.integer(maxvar)
            ,muts=as.integer(numeric(maxvar))
            ,PACKAGE="ARG")

   # cat("mutated tree")
    
    var <- r$varsites
    res <- list(ss=graph$ss
                ,sites=graph$sites
                ,rec=graph$rec
                ,haplotype= matrix(r[[3]][1:(var*samplesize)],ncol=var,byrow=TRUE)
                ,position=r[[4]][1:var]+1
                ,which=graph$which
                ,var=var
                ,muts=r$muts[1:var]
                )
    if (inherits(graph,"structuredARG")) {
      res$location <- graph$location
      class(res) <- c("mutatedstructuredARG","mutatedARG","structuredARG","ARG","GenomicHaplotype")
    } else {
      class(res) <- c("mutatedARG","ARG","GenomicHaplotype")
    }
    res
  }
## a much quicker way to do this - do not try to
## allocate lotes of memory, instead do each node separately
##
"mutatethetapos" <- function(graph,pos,theta)
  {
    samplesize <- sum(graph$ss)
    maxvar <- length(theta)
    if (pos<1||pos>graph$sites)
      stop("pos must be between 1 and",graph$sites)

    r <- .C("mutatethetapos"
            ,as.integer(graph$which)
            ,as.integer(pos-1)
            ,as.double(theta)
            ,data=as.integer(numeric(maxvar*samplesize))
            ,position=as.integer(numeric(maxvar))
            ,varsites=as.integer(maxvar)
            ,muts=as.integer(numeric(maxvar))
            ,PACKAGE="ARG")
    var <- r$varsites
    list(haplotype= matrix(r$data[1:(var*samplesize)],ncol=var,byrow=FALSE)
         ,position=r$position[1:var]+1
         ,muts=r$muts[1:var])
  }
##
##
##
"mutateSTR" <- function(graph,theta,position)
  {
    if (length(theta)>graph$sites) {
       stop("length of theta must be less than or equal to the number of variable sites")
    }
    if (length(theta)!=length(position))       {
       stop("length of theta must be the same as positions given")
    }
    samplesize <- sum(graph$ss)
    var <- length(theta)

    r <- .C("STRmutate"
            ,as.integer(graph$which)
            ,as.double(theta)
            ,as.integer(position-1)
            ,as.integer(var)
            ,as.integer(numeric(var*samplesize))
            ,PACKAGE="ARG")
    
    res <- list(ss=graph$ss,sites=graph$sites
              ,rec=graph$rec
              ,haplotype=matrix(r[[5]],ncol=var,byrow=FALSE)
              ,position=position
              ,which=graph$which
                )
    if (inherits(graph,"structuredARG")) {
      res$location <- graph$location
      class(res) <- c("mutatedstructuredARG","mutatedARG","structuredARG","ARG","GenomicHaplotype")
    } else {
      class(res) <- c("mutatedARG","ARG","GenomicHaplotype")
    }
    res
  }
##
##
##
"mutateAboveNodes" <- function(graph,pos)
  {
    n=sum(graph$ss)
    r <- .C("mutateAboveNodes"
            ,as.integer(graph$which)
            ,as.integer(pos-1)
            ,as.integer(numeric(n*(n-3)))
            ,PACKAGE="ARG");
    matrix(r[[3]],ncol=n-3,byrow=TRUE)
  }
##
##
##
 "prune" <- function(graph,keep)
  {
    if (length(keep)>graph$ss||length(keep)<2)
      stop("array of samples to keep must have length > 2 and less than (or equal to) sample size")
    r <- range(keep)
    if (r[1]<1||r[2]>graph$ss)
      stop("value to keep must be between 1 and sample size");
    
    ukeep <- unique(keep)
    if (length(ukeep)<length(keep))
      warning("repeated values in keep array, removed and continuing")

    cat("starting\n")

    
    .C("pruneARG"
       ,as.integer(graph$which)
       ,as.integer(ukeep-1)
       ,as.integer(length(ukeep))
       ,PACKAGE="ARG")

    cat("finished\n")
    
    graph$keep <- ukeep
    graph$ss <- length(ukeep)
    if (inherits(graph,"mutatedARG"))
      graph$haplotype <- graph$haplotype[ukeep,]
    class(graph) <- c("prunedARG",class(graph))
    graph
  }
##
##
##
### what are the lengths of haplotypes around position
"haplengths" <- function(graph,position)
  {

    l <- .C("haplengths"
            ,as.integer(graph$which)
            ,as.integer(position-1)
            ,as.double(numeric(2*graph$ss))
            ,PACKAGE="ARG")[[3]]
    matrix(l,nrow=graph$ss)
  }

###
"sharedsectionARG" <- function(graph,position,depth=2)
  {
  #  if (is.null(graph$keep))
      wh <- 0:(sum(graph$ss)-1)
  #  else wh <- graph$keep-1;

    if (depth==2) {
      l <- .C("sharedsection"
              ,as.integer(graph$which)
              ,as.integer(wh)
              ,as.integer(length(wh))
              ,as.integer(position-1)
              ,as.double(numeric(2*length(wh)))
              ,PACKAGE="ARG")[[5]]+1
      return(matrix(l,nrow=length(wh)))
    } else {
      l <- .C("sharedsectionK"
              ,as.integer(graph$which)
                 ,as.integer(wh)
                 ,as.integer(length(wh))
                 ,as.integer(position-1)
                 ,as.double(numeric(2*length(wh)))
                 ,as.integer(depth)
              ,PACKAGE="ARG")[[5]]+1
      return(matrix(l,nrow=length(wh)))
    }
  }
##
##
##
##
"ARG.remove" <- function(graph)
  {
    .C("removeARG"
       ,as.integer(graph$which)
       ,PACKAGE="ARG")
  }
    
