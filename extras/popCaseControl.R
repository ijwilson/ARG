.selecthap <- function(u,penetrance,n) {
  n.geno <- length(u)
  pu <- penetrance[u+1]
  repeat {
    sel <-matrix(runif(n.geno)<pu)
  #  cat(sum(sel),n,"\n")
    if (sum(sel)>=n[1]) break
    else warning("did not get enought cases - trying again but perhaps need larger penetrance or more sporadic cases to get ",n)
  }
  list(case=which(sel)[1:n[1]],control=which(!sel)[1:n[2]])
}

.selectdip <- function(u,penetrance,n,offset=0) {
  
    genotypes <- matrix(u,ncol=2,byrow=TRUE)
    indices = matrix(1:length(u),ncol=2,byrow=TRUE)
    isCase <- runif(nrow(genotypes)) < penetrance[rowSums(genotypes)+1] 
    if (sum(isCase)<n[1]) {
      warning("need larger penetrance or more sporadic cases to get ",n)
      return(0)
    }

    case <- as.vector(t(indices[isCase,]))
    

    ret <- list(case=case[1:(2*n[1])]+offset)
    if (n[2]>0) {
      control <- as.vector(t(indices[!isCase,]))
      ret$control <- control[1:(2*n[2])]+offset
    }
   return(ret)
    
  }

## and now a function taken for haploid mitochondrial diseases - hence generally rho=0.0
## note that we shall use theta rather than variable sites for this simulation
## for now use a gamma distribution for the mutation rate with parameters
## thetapars[1] and thetapars[2]
##
PopCaseControlHapsim <-
  function(n=c(100,100),popsize=100,sites=1000,penetrance=c(0.0,1.0),growthmodel="constant"
           ,thetapars=c(1,1),diseasefreq=0.1,fdiff=0.02,rho=0.0,loud=FALSE)
  {
        
    sampsize <- .CaseSampleSize(n[1],diseasefreq-fdiff,penetrance,0.00001)+popsize
    if (n[1]+n[2]+popsize>sampsize) sampsize <- n[1]+n[2]+popsize
    if (loud) {
      pcase <- penetrance*c(1-diseasefreq,diseasefreq)
      cat("Choose a total sample size of ", sampsize,"individuals,\n",
          "we expect to get"
          ,round((sampsize-popsize)*sum(pcase)),"affected individuals\nsampling tree\n")
    }

    a <- simARG(sampsize,sites,rho,growthmodel=growthmodel)
    ## note that the last opopsize samples are the population sample
    theta <- rgamma(sites,thetapars[1],thetapars[2])
       
    repeat {
      b <- mutatetheta(a,theta)
      f1 <- apply(b$haplotype,2,mean)
      w1 <- which(abs(f1-diseasefreq)<fdiff)
      if (length(w1)==0) 
        warning("no sites with the correct frequency")
      if (length(w1)>0) break
    }
    ## find one as close to the centre as possible
    wh <- w1[which.min(abs(w1-sites/2))]
    if (length(wh)>1)
        wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]
    
    u1 <- b$haplotype[,wh][1:(sampsize-popsize)]          ## the disease locus for the CC set
    cc <- .selecthap(u1,penetrance,n)
    pop <- (sampsize-popsize+1):sampsize
    bb <- prune(b,c(cc$case,cc$control,pop))
    bb$muts <- b$muts[wh]
    bb$diseaseSNP <- bb$haplotype[,wh]
    bb$ss <- c(n,popsize)
    bb$haplotype <- bb$haplotype[,-wh]                            ## remove the disease mutation
    bb$diseasepos <- pos
    bb$cc <- cc
    bb$location <- factor(rep(1:3,c(n,popsize)),labels=c("Case","Control","Pop"))
    bb$muts <- b$muts
    bb$position <- bb$position[-wh]
    class(bb) <-c("CCPopHap","CCHap","mutatedstructuredARG",class(b))
    bb
  }
##
## everything is diploids!
PopCaseControlDipsim <-
  function(n=c(100,100),popsize=100,sites=1000,penetrance=c(0.0,0.1,1.0),growthmodel="constant"
           ,thetapars=c(1,1),diseasefreq=0.1,fdiff=0.02,rho=0.0,loud=FALSE)
  {
    sampsize <- 2*.DipCaseSampleSize(n[1],diseasefreq-fdiff,penetrance,0.00001)+2*popsize
    if (n[1]+n[2]+popsize>sampsize) sampsize <- n[1]+n[2]+popsize
    if (loud) {
      pcase <- penetrance*c((1-diseasefreq)^2,2*(1-diseasefreq)*diseasefreq,diseasefreq^2)
      cat("Choose a total sample size of ", sampsize,"individuals,\n",
          "we expect to get"
          ,round((sampsize-popsize)*sum(pcase)),"affected individuals\nsampling tree\n")
    }

    a <- simARG(sampsize,sites,rho,growthmodel=growthmodel)
    ## note that the last opopsize samples are the population sample
    theta <- rgamma(sites,thetapars[1],thetapars[2])

    ## repeat the mutations until we can get at least one with the correct
    ## disease allele frequency
    repeat {
      b <- mutatetheta(a,theta)
      f1 <- apply(b$haplotype,2,mean)
      w1 <- which(abs(f1-diseasefreq)<fdiff)
      if (length(w1)==0) 
        warning("no sites with the correct frequency")
      if (length(w1)>0) break
    }
    ## find one as close to the centre as possible
    wh <- w1[which.min(abs(w1-sites/2))]
    if (length(wh)>1)
        wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]
    
    u1 <- b$haplotype[,wh][1:(sampsize-popsize)]          ## the disease locus for the CC set
    cc <- .selectdip(u1,penetrance,n)
    pop <- (sampsize-(2*popsize)+1):sampsize
    bb <- prune(b,c(cc$case,cc$control,pop))

    
    bb$diseasemuts <- b$muts[wh]
    bb$diseaseSNP <- bb$haplotype[,wh]
    bb$diseasepos <- pos

    bb$ss <- c(n,popsize)
    bb$haplotype <- bb$haplotype[,-wh]            ## remove the disease mutation
    bb$cc <- cc
    bb$diseaseTree <- tree(bb,pos)
    bb$location <- factor(rep(1:3,c(2*n,2*popsize)),labels=c("Case","Control","Pop"))
    bb$muts <- b$muts[-wh]
    bb$position <- bb$position[-wh]
    class(bb) <-c("CCPopDip","CCDip","mutatedstructuredARG",class(b))
    bb
  }
