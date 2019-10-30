## How may samples do I need to take to ensure that
## I get enough cases for a haploid model?
".CaseSampleSize" <-
function(n,minfreq,penetrance,pfail=0.00001)
  {  
    if (minfreq>0.5)
      stop("Expect the minimum population frequency to be less than 0.5")
    if (length(penetrance) !=2)
      stop("expect two penetrance parameters in .CaseSampleSize")
    if (penetrance[2]<penetrance[1])
      stop("penetrance must be greater than proportion of sporadic cases")
    pcase <- penetrance[2]*minfreq+penetrance[1]*(1-minfreq)
    q <- qnorm(pfail)
  
    a <- polyroot(c(n*n,-2*n*pcase-q*q*pcase*(1-pcase),pcase*pcase))
    ceiling(Re(a[2]))
  }
## returns the *Diploid" sample size required for a simulation
".DipCaseSampleSize" <-
function(n,minfreq,penetrance,pfail=0.00001)
  {
    if (minfreq>0.5)
      stop("Expect the minimum population frequency to be less than 0.5")
    if (length(penetrance) !=3)
      stop("expect three penetrance parameters in .DipCaseSampleSize")
    
    pcase <-  penetrance[3]*minfreq*minfreq+2*penetrance[2]*(1-minfreq)*minfreq
    +penetrance[1]*(1-minfreq)*(1-minfreq)

  #  cat(pcase,"\n")
    q <- qnorm(pfail)
    
   a <- polyroot(c(n*n,-2*n*pcase-q*q*pcase*(1-pcase),pcase*pcase))
    ceiling(Re(a[2]))
    
  }
##
## n and popsize are the number of diploids
##
"CaseControlDipTree" <-
  function(n=c(20,20),penetrance=c(0.00,0.7,1.0),sites=2010,var=200,growthmodel="constant"
           ,freq=0.1,fdiff=0.02,rho=0.1,popsize=0,loud=TRUE,age=0)
  {
    
    SampleSize <- .DipCaseSampleSize(n[1],freq-fdiff,penetrance)+2*popsize
    ss <- 2*SampleSize
    
    if (loud) cat("Choose a total sample size of ", SampleSize,"diploids,\nsimulating tree ...\n")
    a <- simARG(ss,sites,rho,growthmodel=growthmodel)
  
    if (age>0) ascert="older(age)"
    else ascert="bylength"

    repeat {
      if (loud) cat("mutating ...\n")
      b <- mutate(a,var,ascert=ascert)
      f1 <- colMeans(b$haplotype)
      
      w1 <- which(abs(f1-freq)<fdiff)
      if (length(w1)==1)
         warning("no sites with the correct frequency")
       if (length(w1)>0) break
    }
    ## find one as close to the centre as possible
    wh <- w1[which(abs(w1-var/2) == min(abs(w1-var/2)))]
    if (length(wh)>1)
      wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]

    maxCC <- ss-2*popsize
    
    u1 <- b$haplotype[,wh][1:maxCC]   ## the disease locus for the CC set 
    if (loud) cat("selecting cases ...\n")
    cc <- .selectdip(u1,penetrance,n)
    if (popsize>0) pop <- maxCC:ss
    else pop <- NULL
    if (loud) cat("pruning ARG ...\n")
    bb <- prune(b,c(cc$case,cc$control,pop))
   
    bb$diseaseSNP <- bb$haplotype[,wh]
    bb$ss <- c(n,popsize)
    bb$diseasepos <- pos
    bb$cc <- cc
    if (popsize>0) 
      bb$location <- factor(rep(1:3,c(2*n,2*popsize)),labels=c("Case","Control","Pop"))
    else bb$location <- factor(rep(1:2,2*n),labels=c("Case","Control"))
    ## b$position <- b$position[-wh]
    class(bb) <-c("CCPopDip","CCDip","mutatedstructuredARG","prunedARG",class(b))
    bb

  }
##
## n and popsize are the number of diploids
##
"CaseControlDipSamples" <-
  function(n=c(20,20),penetrance=c(0.00,0.7,1.0),sites=2010,var=200,growthmodel="constant"
           ,freq=0.1,fdiff=0.02,rho=0.1,popsize=0,loud=TRUE,age=0)
  {
    
    SampleSize <- .DipCaseSampleSize(n[1],freq-fdiff,penetrance)+2*popsize
    ss <- 2*SampleSize
    
    if (loud) cat("Choose a total sample size of ", SampleSize,"diploids,\nsimulating tree ...\n")
    a <- simARG(ss,sites,rho,growthmodel=growthmodel)
  
    if (age>0) ascert="older(age)"
    else ascert="bylength"

    repeat {
      if (loud) cat("mutating potential disease sites ...\n")
      b <- mutate(a,var,ascert=ascert)
      f1 <- colMeans(b$haplotype)
      
      w1 <- which(abs(f1-freq)<fdiff)
      if (length(w1)==1)
         warning("no sites with the correct frequency")
       if (length(w1)>0) break
    }
    ## find one as close to the centre as possible
    wh <- w1[which(abs(w1-var/2) == min(abs(w1-var/2)))]
    if (length(wh)>1)
      wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]

    maxCC <- ss-2*popsize
    
    u1 <- b$haplotype[,wh][1:maxCC]   ## the disease locus for the CC set 
    if (loud) cat("selecting cases ...\n")
    cc <- .selectdip(u1,penetrance,n)
    if (popsize>0)
      cc$pop <- (maxCC+1):ss
    cc$diseasepos <- pos
    cc$diseaseSNP <- b$haplotype[c(cc$case,cc$control,cc$pop),wh]
    cc
  }
##
##
##
"CaseControlDipsim" <-
  function(n=20,penetrance=c(0.00,0.7,1.0),sites=2010,var=200,growthmodel="constant"
           ,freq=0.1,rho=0.1,panel)
  {
    selectdip <- function(u,penetrance,n) {
      genotypes <- matrix(u,ncol=2,byrow=T)
      n.geno <- nrow(genotypes)
      m <- rowSums(genotypes)+1 
      pu <- penetrance[m]
      sel <-matrix(runif(n.geno)<pu) 
      
      if (sum(sel)<n)
        stop("need larger penetrance or more sporadic cases to get ",n)
      seli <- as.vector(t(matrix(sel,nrow=n.geno,ncol=2))) ## double up
      
      list(case=which(seli)[1:n],control=which(!seli)[1:n])
    }
    
    SampleSize <- .DipCaseSampleSize(n,freq,penetrance)
    ss <- 2*SampleSize
   
    
    cat("Choose a total sample size of ", ss/2,"individuals,\n")
    a <- simARG(ss,sites,rho,growthmodel=growthmodel)
    if (missing(panel))
      ascert <- paste("panel(",as.integer(n/2),")",sep="")
    else
      ascert= paste("panel(",as.integer(panel),")",sep="")
    
    cat("mutating tree with ascertainment",ascert,"\n")
    b <- mutate(a,(var+1),ascert=ascert)
    
    f1 <- apply(b$haplotype,2,mean)
     
    w1 <- which(abs(f1-freq)<0.02)
    if (length(w1)==1)
      stop("no sites with the correct frequency\n")
    wh <- w1[which(abs(w1-var/2) == min(abs(w1-var/2)))]
    if (length(wh)>1)
      wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]

    b$diseaseSNP <- b$haplotype[,wh]
    cc <- selectdip(b$diseaseSNP,penetrance,n)
    b <- prune(a,c(cc$case,cc$control))
  
    b$ss <- c(n,n)
    b$haplotype <- b$haplotype[,-wh]
    b$diseasepos <- pos
    b$cc <- cc
    b$location <- factor(rep(1:2,c(n,n)),labels=c("Case","Control"))
   
    b$position <- b$position[-wh]
    class(b) <-c("CCDip","mutatedstructuredARG",class(b))
    b
  }
##
##
##




## and now a function taken for haploid mitochondrial diseases - hence generally rho=0.0
CaseControlHapsim <-
  function(n=60,penetrance=c(0.0,1.0),sites=2010,var=200,growthmodel="constant"
           ,freq=0.1,rho=0.0,panel,fdiff=0.02,loud=FALSE)
  {
    require(ARG)
    selecthap <- function(u,penetrance,n) {
      n.geno <- length(u)
      pu <- penetrance[u+1]

      repeat {
        sel <-matrix(runif(n.geno)<pu) 
        if (sum(sel)>=n) break
        else warning("did not get enought cases - trying again but perhaps need larger penetrance or more sporadic cases to get ",n)
      }
      list(case=which(sel)[1:n],control=which(!sel)[1:n])
    }

    
    pcase <-     penetrance*c(1-freq,freq)
    pcontrol <- (1-penetrance)*c(1-freq,freq)
    ss <- max(c(n/sum(pcase), n/sum(pcontrol)))
    ss <- round(ss+8*sqrt(ss))
    if (loud)
      cat("Choose a total sample size of ", ss,"individuals,\n",
          "we expect to get"
          ,round(ss*sum(pcase)),"affected individuals\nsampling tree\n")

    a <- simARG(ss,sites,rho,growthmodel=growthmodel)
    if (missing(panel))
      ascert <- paste("panel(",as.integer(n/2),")",sep="")
    else
      ascert= paste("panel(",as.integer(panel),")",sep="")
    
    if (loud) cat("mutating tree with ascertainment",ascert,"\n")
  
    repeat {
      b <- mutate(a,(var+1),ascert=ascert)
      f1 <- apply(b$haplotype,2,mean)
      w1 <- which(abs(f1-freq)<fdiff)
      if (length(w1)==0) 
        warning("no sites with the correct frequency")
      if (length(w1)>0) break
    }
    
    wh <- w1[which(abs(w1-var/2) == min(abs(w1-var/2)))]
    if (length(wh)>1)
        wh <- sample(wh,1) ## possible to get two equally spaced
    pos <- b$position[wh]
    
    u1 <- b$haplotype[,wh]
    cc <- selecthap(u1,penetrance,n)
  
    bb <- prune(b,c(cc$case,cc$control))
   # print(class(b))
    bb$diseaseSNP <- bb$haplotype[,wh]
    bb$ss <- c(n,n)
    bb$haplotype <- bb$haplotype[,-wh]                            ## remove the disease mutation
    bb$diseasepos <- pos
    bb$cc <- cc
    bb$location <- factor(rep(1:2,c(n,n)),labels=c("Case","Control"))
   
    bb$position <- bb$position[-wh]
    class(bb) <-c("CCHap","mutatedstructuredARG",class(b))
    bb
  }
##
##
##
fixphylo <- function(phy,cc,casechar="*")
  {
    lab <- as.numeric(phy$tip.label)+1
    case <- lab %in% cc$case
    newlab <- lab
    newlab[case] <- casechar
    newlab[!case] <- " "
    phy$tip.label <- newlab
    phy
  }
## produce a vector of p-values for case control data from sima
CaseControl <- function(sar,test="Chisq")
{
  cc <- function(cl) {
    chisq.test(table(cl,sar$location))$p.value
  }
  ccf <- function(cl) {
    t <- table(cl,sar$location)
    fisher.test(t)$p.value
  }
  if (test=="Chisq")
    return(apply(sar$haplotype,2,cc))
  else if (test=="Fisher"||test=="fisher")
    return(apply(sar$haplotype,2,ccf))
  
}
