
displaymutatedARG <- function(b,interactive=FALSE)
  {
    if(length(b$position)!=ncol(b$haplotype))
      stop("problem with data - number of columns (",ncol(b$haplotype)
           ,") does not match positions (",length(b$position),").\n")
    rows <- nrow(b$haplotype)
    columns <- ncol(b$haplotype)
    getlabels <- function(ttp) {
      labels <- as.numeric(ttp$tip.label)
      }
    
    plotshared <- function(m,currorder,pos) {
      par(mar=c(0,0,0,0))
      y <- 1:nrow(b$haplotype)
      mo <- m[currorder,]
      plot(-100,-100,ylim=c(1,rows),xlim=c(1,b$sites),axes=F,ylab="",xlab="")
      segments(mo[,1],y,mo[,2],y,col="blue",lwd=2)
      segments(1,1,1,rows,col="lightblue");
      segments(b$sites,1,b$sites,rows,col="lightblue");
      segments(pos,1,pos,rows,col="salmon");
    }
    
    plotlen <- function(pos) {
      screen(5,TRUE)
      par(mar=c(0,2,0,0))
      plot(1:b$sites,h,type="l",lwd=2,axes=F)
      segments(pos,0,pos,10,col="salmon");
      axis(2)
      screen(6,TRUE)
      par(mar=c(0,2,0,0))
      plot(1:b$sites,l,type="l",lwd=2,axes=F)
      segments(pos,0,pos,30,col="salmon");
      axis(2)
      rug(b$position,lwd=2);
    }
    
    treeplot <- function(up) {  
            
      if (up<1||up>b$sites)
        cat("this position outside region, try again")
      
      screen(3,new=T); 
      par(mar=c(0,0,0,0))
      tp <- tree(b,up)
      
      currentorder <- getlabels(tp)
      tp$tip.label <- currentorder
      plot(tp,edge.width=2,bg="white",x.lim=th)
      plotlen(up);
      
      m <- sharedsectionARG(b,up);
      
      screen(2,TRUE)
     # print(sort(currentorder))
      genoplot(currentorder,m)
      
    }

    genoplot <- function(currorder,m) {
      posmat <- matrix(b$position,nrow=rows,ncol=columns,byrow=TRUE)
      mod <- posmat>m[currorder,1]&posmat<m[currorder,2]
      dat <- b$haplotype[currorder,] + mod*0.1
      if (inherits(b,"mutatedstructuredARG")) {
          
        loc <- matrix(as.numeric(b$location[currorder]),nrow=rows,ncol=columns)
        
        if (length(table(b$location)==2)) {
          cols <- c("white","wheat1","cyan","black")
          breaks <- c(-1,0.015,0.025,0.15,1.15)
        } else if (length(table(b$location)==3)) {
          cols <- c("white","wheat1","skyblue","cyan","black")
          breaks <- c(-1,0.015,0.025,0.035,0.15,1.15)
        } else  stop("can only have 3 locations")
        dat <- dat+loc*0.01
        d <- t(dat)          
        image(x=1:columns,y=1:rows,d,axes=F,xlab="",ylab="",col=cols,breaks=breaks)
        } else {
          d <- t(dat)
          image(x=1:columns,y=1:rows,d,axes=F,xlab="",ylab=""
                ,col=c("white","cyan","black"),breaks=c(-1,0.05,0.15,1.15))
        }
       
      axis(2,at=1:rows,labels=currorder,tick=F,lwd=0,las=1,cex.lab=1.2)
    }
    
    l <- treelength(b);
    h <- treeheight(b)    
    
    if (interactive) {
      split.screen(c(2,1));
      split.screen(c(1,2),screen=1)
      split.screen(c(2,1),screen=4)
      screen(2);
    }
    th <- max(treeheight(b))*1.1;
    opar <- par(mar=c(0,2,0,0),bg="white")
    if (inherits(b,"CCDip")) {
      start <- b$diseasepos
    } else start=b$sites/2
    treeplot(start)
     
    if (interactive) {
      while(TRUE) {
        v <- locator(1)
        if (v$y>sum(b$ss)) break
        usepos <- round(v$x,0)
        treeplot(b$position[usepos])
      }
    }
    
    par(opar)
  }


displaymutatedARGsmall <- function(b,interactive=TRUE,highlight)
  {
    getlabels <- function(tp) {
        labels <- as.numeric(tp$tip.label)+1
      }
    
    plotshared <- function(m,currorder,pos) {
      par(mar=c(0,0,0,0))
      y <- 1:sum(b$ss)
      mo <- m[currorder,]
      plot(-100,-100,ylim=c(1,sum(b$ss)),xlim=c(1,b$sites),axes=F,ylab="",xlab="")
      segments(mo[,1],y,mo[,2],y,col="blue",lwd=2)
      segments(1,1,1,sum(b$ss),col="lightblue");
      segments(b$sites,1,b$sites,sum(b$ss),col="lightblue");
      segments(pos,1,pos,sum(b$ss),col="salmon");
    }
    
    plotlen <- function(pos) {
      screen(5,TRUE)
      par(mar=c(0,2,0,0))
      plot(1:b$sites,h,type="l",lwd=2,axes=F)
      segments(pos,0,pos,10,col="salmon");
      axis(2)
      screen(6,TRUE)
      par(mar=c(0,2,0,0))
      plot(1:b$sites,l,type="l",lwd=2,axes=F)
      segments(pos,0,pos,30,col="salmon");
      axis(2)
      rug(b$position,lwd=2);
    }
    
    treeplot <- function(m) {  
      usepos <- round(v$x,0)
      usepos <- b$position[usepos]
     
      if (usepos<1||usepos>b$sites)
        cat("this position outside region, try again")
      else {
        screen(3,new=T); 
        par(mar=c(0,0,0,0))
        tp <- tree(b,usepos)
        currentorder <- getlabels(tp)
        tp$tip.label <- currentorder
        plot(tp,edge.width=2,bg="white",x.lim=th)
        plotlen(usepos);
        m <- sharedsectionARG(b,usepos)
        screen(2,TRUE)
        genoplot(currentorder,m)
      }
    }

    
    genoplot <- function(currorder,m) {
      
        plot(-100,-1,ylim=c(1,sum(b$ss)),xlim=c(1,b$sites),axes=F,xlab="",ylab="",xaxs="i")
        abline(h=1:sum(b$ss),col="lightcyan",lwd=4)
        abline(v=b$position,col="lightcyan",lwd=1)
        y <- matrix((1:sum(b$ss)),nrow=sum(b$ss),ncol=length(b$position))[b$haplotype==1]
        x <- matrix(b$position,nrow=sum(b$ss),ncol=length(b$position),byrow=T)[b$haplotype==1]
        if (!missing(highlight)) {
          col1 <-  matrix(highlight,nrow=sum(b$ss),ncol=length(b$position),byrow=T)
          col2 <- matrix("black",nrow=sum(b$ss),ncol=length(b$position))
          col2[col1] <- "red"
          col2 <- col2[b$haplotype==1]
        } else col2 <- "black"
        points(x,y,pch=16,cex=0.3,col=col2)
     
      axis(2,at=1:sum(b$ss),labels=currorder,tick=F,lwd=0,las=1,cex.lab=1.2)
    }
    
    l <- treelength(b);
    h <- treeheight(b)    
    
    if (interactive) {
      split.screen(c(2,1));
      split.screen(c(1,2),screen=1)
      split.screen(c(2,1),screen=4)
      screen(2);
    }
    th <- max(treeheight(b))*1.1;
    opar <- par(mar=c(0,2,0,0),bg="white")
    if (inherits(b,"CCDip")) {
      start <- b$DiseasePos
    } else start=b$sites/2
    treeplot(1:sum(b$ss))
    m <- sharedsectionARG(b,start)
    genoplot(1:sum(b$ss),m);
    
    if (interactive) {
      while(TRUE) {
        v <- locator(1);
        
        if (v$y>sum(b$ss)) break
        treeplot(v);
        
      }
    }
    
    par(opar)
  }
