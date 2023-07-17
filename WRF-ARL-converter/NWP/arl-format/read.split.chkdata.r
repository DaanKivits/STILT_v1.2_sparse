audit.all.split.chkdata <- function(dir,verbose=FALSE) {
  cat('audit results for dir=',dir,'\n')
  test <- read.all.split.chkdata(dir,audit.only=TRUE,verbose=verbose)
  if (length(test) > 0) {
    badtimes <- NULL
    for (x in names(test)) badtimes <- unique(c(badtimes,names(test[[x]])))
    for (badtime in sort(badtimes)) {
      cat(badtime,':')
      for (x in names(test)) if (badtime %in% names(test[[x]])) cat(' ',x,'=',test[[x]][badtime],sep='')
      cat('\n')
    }
  } else {
    cat('Found no errors\n')
  }
}
  
plot.all.split.chkdata <- function(tmp.list,tmp.list2=NULL,ps.or.pdf='pdf',type='l',...) {
  #!# $Id: read.split.chkdata.r,v 1.11 2016/03/15 18:21:37 trn Exp $
  # Plot all output from split_chk_data.csh
  # Input: tmp.list - output from read.all.split.chkdata
  for (tmp.var in names(tmp.list)) {
    levs.all <- NULL
    if (tmp.var == "NULL" || is.null(tmp.list[[tmp.var]])) {
      cat ("Warning: tmp.list contains empty component named",tmp.var,"indicates\n",
           "         problem in generation of ARL file\n")
    } else {
      if (dim(tmp.list[[tmp.var]])[1] > 1) {
        levs <- 1:dim(tmp.list[[tmp.var]])[1]
        if (is.null(levs.all)) levs.all <- levs
# ignore top level of convective up/downdraft fluxes (should always be zero)
        if ((tmp.var == 'CFU1' || tmp.var == 'CFD1') && length(levs) == length(levs.all))
          levs <- levs[-length(levs)]
        plot.6panel.split.chkdata(tmp.list,tmp.var,tmp.list2=tmp.list2,levs=levs,ps.or.pdf=ps.or.pdf,type=type,...)
      } else {
        plot.3panel.split.chkdata(tmp.list,tmp.var,tmp.list2=tmp.list2,ps.or.pdf=ps.or.pdf,type=type,...)
      }
    }
  }
}

read.all.split.chkdata <- function(dir,template='chk_data_????.txt',audit.only=FALSE,...) {
  #!# $Id: read.split.chkdata.r,v 1.11 2016/03/15 18:21:37 trn Exp $
  # Return all output from split_chk_data.csh in one list object
  # Input:
  # dir - directory name of input split_chk_data.csh output files
  if (substring(dir,nchar(dir),nchar(dir)) != '/') dir <- paste(dir,'/',sep='')
  ifirst <- 0
  ilast <- 0
  for (i in 1:nchar(template)) {
    if (substring(template,i,i) == '?') {
      if (ifirst == 0) ifirst <- i + nchar(dir)
      ilast <- i + nchar(dir)
    }
  }
  fnames <- system(paste('/bin/ls -1 ',dir,template,sep=''),intern=T)
  vars <- substring(fnames,ifirst,ilast)
  out.list <- vector('list',length(vars))
  names(out.list) <- vars
  for (i in 1:length(vars)) {
    tmp.var <- vars[i]
    if (audit.only) {
      out.list[[tmp.var]] <- audit.one.split.chkdata(tmp.arr=NULL,fname=fnames[i],...)
    } else {
      out.list[[tmp.var]] <- read.one.split.chkdata(fnames[i],...)
    }
  }
  out.list
}

read.one.split.chkdata <- function(fname,
                                   col.names=c('yr','mo','da','hr','fh','lev','fmn','min','max','mean'),
                                   col.widths=c(rep(2,6),3,rep(15,3)),remove.repeated.times=TRUE,
                                   verbose=FALSE)
{
  #!# $Id: read.split.chkdata.r,v 1.11 2016/03/15 18:21:37 trn Exp $
  tmp.what <- vector('list',length(col.names))
  for (i in 1:length(col.names)) tmp.what[[i]] <- numeric()
  tmp.dat <- scan(fname,what="",sep='\n',quiet=!verbose)
  nchar.tmp <- nchar(tmp.dat[1])
  if (substring(tmp.dat[1],11,11) == 'X') {
    tmp.dat <- paste(substring(tmp.dat,1,10),substring(tmp.dat,12,nchar.tmp),sep='')
    col.widths=c(rep(2,5),3,3,rep(15,3))
  }
  
  tmp.arr <- array(NA,dim=c(length(tmp.dat),length(col.widths)),dimnames=list(NULL,
                                                                    col.names))
  cum.widths <- c(0,cumsum(col.widths))
  for (i in 1:dim(tmp.arr)[2]) tmp.arr[,i] <- as.numeric(substring(tmp.dat,cum.widths[i]+1,cum.widths[i+1]))
  out.arr <- NULL
  keep.trying <- T
  kount <- 0
  while (keep.trying) {
    kount <- kount+1
    tmp.levs <- unique(tmp.arr[,'lev'])
    tmp.times <- (paste(tmp.arr[,'yr'],tmp.arr[,'mo'],tmp.arr[,'da'],tmp.arr[,'hr'],tmp.arr[,'fh'],tmp.arr[,'fmn'],sep='.'))[seq(1,dim(tmp.arr)[1],by=length(tmp.levs))]
    if (length(tmp.times)*length(tmp.levs) != dim(tmp.arr)[1]) {
      cat('read.one.split.chkdata error for file',fname,'\n')
      cat('length(tmp.levs)=',length(tmp.levs),' length(tmp.times)=',length(tmp.times),
          'dim(tmp.arr)=',dim(tmp.arr),'\n')
      if (kount == 1 && length(tmp.levs) > 1) {
    # try it without the top level:
        top.lev <- tmp.levs[length(tmp.levs)]
        cat ('Trying again without the top level=',top.lev,'\n')
        tmp.arr <- tmp.arr[tmp.arr[,'lev'] != top.lev,]
        tmp.levs <- unique(tmp.arr[,'lev'])
        tmp.times <- (paste(tmp.arr[,'yr'],tmp.arr[,'mo'],tmp.arr[,'da'],tmp.arr[,'hr'],tmp.arr[,'fh'],tmp.arr[,'fmn'],sep='.'))[seq(1,dim(tmp.arr)[1],by=length(tmp.levs))]
        keep.trying <- T
      } else {
        cat ('Could not be fixed, not stored\n')
        audit.one.split.chkdata(tmp.arr)
        keep.trying <- F
      }
    } else {
      out.arr <- array(as.vector(tmp.arr[,c('min','max','mean')]),dim=c(length(tmp.levs),length(tmp.times),3),
                       dimnames=list(as.character(tmp.levs),tmp.times,c('min','max','mean')))
      keep.trying <- F
    }
  }
  if (remove.repeated.times && !is.null(out.arr)) {
    repeated.times <- c(FALSE,dimnames(out.arr)[[2]][-1]==dimnames(out.arr)[[2]][-dim(out.arr)[2]])
    out.arr <- out.arr[,!repeated.times,,drop=FALSE]
  }
  out.arr
}

audit.one.split.chkdata <- function(tmp.arr=NULL,fname=NULL,
                                   col.names=c('yr','mo','da','hr','fh','lev','fmn','min','max','mean'),
                                   col.widths=c(rep(2,6),3,rep(15,3)),verbose=TRUE)
{
  #!# $Id: read.split.chkdata.r,v 1.11 2016/03/15 18:21:37 trn Exp $
  if (is.null(tmp.arr)) {
    if (is.null(fname)) stop('audit.one.split.chkdata: Need to supply non-NULL tmp.arr or fname')
    tmp.what <- vector('list',length(col.names))
    for (i in 1:length(col.names)) tmp.what[[i]] <- numeric()
    tmp.dat <- scan(fname,what="",sep='\n',quiet=!verbose)
    nchar.tmp <- nchar(tmp.dat[1])
    if (substring(tmp.dat[1],11,11) == 'X') {
      tmp.dat <- paste(substring(tmp.dat,1,10),substring(tmp.dat,12,nchar.tmp),sep='')
      col.widths=c(rep(2,5),3,3,rep(15,3))
    }
  
    tmp.arr <- array(NA,dim=c(length(tmp.dat),length(col.widths)),dimnames=list(NULL,
                                                                    col.names))
    cum.widths <- c(0,cumsum(col.widths))
    for (i in 1:dim(tmp.arr)[2]) tmp.arr[,i] <- as.numeric(substring(tmp.dat,cum.widths[i]+1,cum.widths[i+1]))
  }
  all.times <- paste(sprintf('%2.2i',tmp.arr[,'yr']),
                     sprintf('%2.2i',tmp.arr[,'mo']),
                     sprintf('%2.2i',tmp.arr[,'da']),
                     sprintf('%2.2i',tmp.arr[,'hr']),tmp.arr[,'fh'],tmp.arr[,'fmn'],sep='.')
  tmp.levs <- unique(tmp.arr[,'lev'])
  lentimes <- tapply(tmp.arr[,'lev'],all.times,FUN=function(x) length(x),simplify=TRUE)
  lenrange <- range(lentimes)
  if (verbose) cat('Audit results')
  if (!is.null(fname) && verbose) cat(' for fname=',fname,sep='')
  if (verbose) cat(':\n rows=',length(all.times),' levs=',length(tmp.levs),' times=',length(lentimes),
                   ' from ',all.times[1],' to ',all.times[length(all.times)],'\n',sep='')
  lenbreaks <- seq(lenrange[1]-0.5,lenrange[2]+0.5,by=1)
  lenhist <- hist(lentimes,breaks=lenbreaks,plot=FALSE)
  lensumm <- cbind(mids=lenhist$mids,counts=lenhist$counts)
  lensumm <- lensumm[lensumm[,'counts'] > 0,,drop=FALSE]
  if (verbose) print(lensumm)
  badlevs <- lensumm[,'mids'][!(lensumm[,'mids'] %in% c(length(tmp.levs),2*length(tmp.levs)))]
  badtimes <- NULL
  if (length(badlevs) == 0) {
    if (verbose) cat('ok: all times have either levs or 2*levs (06Z) rows\n')
  } else {
    status <- 'bad'
    if (verbose) cat('not ok: following times have unexpected numbers of levs\n')
    badtimes <- lentimes[lentimes %in% badlevs]
    if (verbose) print(badtimes)
  }
  badtimes
}

plot.3panel.split.chkdata <- function(tmp.list=NULL,tmp.var=NULL,tmp.arr=tmp.list[[tmp.var]],
                                      tmp.list2=NULL,tmp.arr2=tmp.list2[[tmp.var]],
                                      diff.flag=FALSE,ps.or.pdf='ps',
                                      psname=paste(tmp.var,'3panel',ps.or.pdf,sep='.'),
                                      xlab='Time Period',panels=c('min','max','mean'),
                                      ylabs=paste(tmp.var,panels),...) {
  if (is.null(tmp.arr)) stop ('Missing args: either tmp.list and tmp.var, or tmp.arr are needed\n')
  if (!is.null(psname)) {
    valid.ps.or.pdf <- FALSE
    if (ps.or.pdf == 'ps') {
      valid.ps.or.pdf <- TRUE
      postscript(psname,horizontal=F,paper='letter')
    }
    if (ps.or.pdf == 'pdf') {
      valid.ps.or.pdf <- TRUE
      pdf(psname,paper='letter',width=0,height=0)
    }
    if (!valid.ps.or.pdf) stop('Invalid ps.or.pdf=',ps.or.pdf)
    par(mfrow=c(length(panels),1))
  }
  for (i in 1:length(panels)) {
    x <- panels[i]
    yvals <- tmp.arr['0',,x]
    ylim <- range(yvals,na.rm=TRUE)
    yvals2 <- NULL
    if (!is.null(tmp.arr2)) {
      yvals2 <- tmp.arr2['0',,x]
      if (diff.flag) {
        yvals <- yvals-yvals2
        ylim <- range(0,yvals,na.rm=TRUE)
      } else {
        ylim <- range(ylim,yvals2,na.rm=TRUE)
      }
    }
    plot(1:dim(tmp.arr)[2],yvals,ylim=ylim,,xlab='Time period',ylab=ylabs[i],...)
    if (!is.null(tmp.arr2) && !diff.flag) points(1:dim(tmp.arr2)[2],yvals2,pch=2,col=2)
  }
  if (!is.null(psname)) dev.off()
}

plot.6panel.split.chkdata <- function(tmp.list=NULL,tmp.var=NULL,tmp.arr=tmp.list[[tmp.var]],
                                      levs=1:dim(tmp.arr)[1],
                                      tmp.list2=NULL,tmp.arr2=tmp.list2[[tmp.var]],
                                      diff.flag=FALSE,ps.or.pdf='ps',
                                      psname=paste(tmp.var,'6panel',ps.or.pdf,sep='.'),
                                      xlab='Time Period',panels=c('min','max','mean'),
                                      ylabs=paste(tmp.var,panels),zlab='Level',...) {
  if (is.null(tmp.arr)) stop ('Missing args: either tmp.list and tmp.var, or tmp.arr are needed\n')
  if (!is.null(psname)) {
    valid.ps.or.pdf <- FALSE
    if (ps.or.pdf == 'ps') {
      valid.ps.or.pdf <- TRUE
      postscript(psname,horizontal=F,paper='letter')
    }
    if (ps.or.pdf == 'pdf') {
      valid.ps.or.pdf <- TRUE
      pdf(psname,paper='letter',width=0,height=0)
    }
    if (!valid.ps.or.pdf) stop('Invalid ps.or.pdf=',ps.or.pdf)
    par(mfrow=c(length(panels),2))
  }
  for (i in 1:length(panels)) {
    x <- panels[i]
    yvals <- apply(tmp.arr[levs,,x,drop=F],MAR=2,FUN=x)
    ylim <- range(yvals,na.rm=TRUE)
    yvals2 <- NULL
    if (!is.null(tmp.arr2)) {
      yvals2 <- apply(tmp.arr2[levs,,x,drop=F],MAR=2,FUN=x)
      ylim <- range(ylim,yvals2,na.rm=TRUE)
    }
    plot(1:dim(tmp.arr)[2],yvals,ylim=ylim,xlab='Time period',ylab=ylabs[i],...)
    if (!is.null(tmp.arr2)) points(1:dim(tmp.arr2)[2],yvals2,pch=2,col=2)
    yvals <- apply(tmp.arr[levs,,x,drop=F],MAR=1,FUN=x)
    ylim <- range(yvals,na.rm=TRUE)
    yvals2 <- NULL
    if (!is.null(tmp.arr2)) {
      yvals2 <- apply(tmp.arr2[levs,,x,drop=F],MAR=1,FUN=x)
      if (diff.flag) {
        yvals <- yvals-yvals2
        ylim <- range(0,yvals,na.rm=TRUE)
      } else {
        ylim <- range(ylim,yvals2,na.rm=TRUE)
      }
    }
    plot(yvals,as.numeric(dimnames(tmp.arr)[[1]])[levs],xlim=ylim,
         ylab=zlab,xlab=ylabs[i],...)
    if (!is.null(tmp.arr2) && !diff.flag)
      points(yvals2,as.numeric(dimnames(tmp.arr2)[[1]])[levs],pch=2,col=2)
  }
  if (!is.null(psname)) dev.off()
}

