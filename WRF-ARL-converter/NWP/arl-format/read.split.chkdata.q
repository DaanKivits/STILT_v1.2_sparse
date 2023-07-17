plot.all.split.chkdata <- function(tmp.list,...) {
  #!# $Id: read.split.chkdata.q,v 1.2 2006/01/16 13:04:08 trn Exp $
  # Plot all output from split_chk_data.csh
  # Input: tmp.list - output from read.all.split.chkdata
  for (tmp.var in names(tmp.list)) {
    if (dim(tmp.list[[tmp.var]])[1] > 1) {
      levs <- 1:dim(tmp.list[[tmp.var]])[1]
# ignore top level of convective up/downdraft fluxes (should always be zero)
      if (tmp.var == 'CFU1' || tmp.var == 'CFD1') levs <- levs[-length(levs)]
      plot.6panel.split.chkdata(tmp.list,tmp.var,levs=levs,...)
    } else {
      plot.3panel.split.chkdata(tmp.list,tmp.var,...)
    }
  }
}


read.all.split.chkdata <- function(dir,template='chk_data_????.txt',...) {
  #!# $Id: read.split.chkdata.q,v 1.2 2006/01/16 13:04:08 trn Exp $
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
  fnames <- unix(paste('/bin/ls -1 ',dir,template,sep=''))
  vars <- substring(fnames,ifirst,ilast)
  out.list <- vector('list',length(vars))
  names(out.list) <- vars
  for (i in 1:length(vars)) {
    tmp.var <- vars[i]
    out.list[[tmp.var]] <- read.one.split.chkdata(fnames[i],...)
  }
  out.list
}

read.one.split.chkdata <- function(fname,
                                   col.names=c('yr','mo','da','hr','fh','lev','min','max','mean'),
                                   col.widths=c(rep(2,6),rep(15,3))
                                   )
{
  #!# $Id: read.split.chkdata.q,v 1.2 2006/01/16 13:04:08 trn Exp $
  tmp.what <- vector('list',length(col.names))
  for (i in 1:length(col.names)) tmp.what[[i]] <- numeric()
  tmp.dat <- scan(fname,what=tmp.what,widths=col.widths)
  tmp.arr <- array(NA,dim=c(length(tmp.dat[[1]]),length(tmp.dat)),dimnames=list(NULL,
                                                                    col.names))
  for (i in 1:dim(tmp.arr)[2]) tmp.arr[,i] <- tmp.dat[[i]]
  tmp.levs <- unique(tmp.arr[,'lev'])
  tmp.times <- unique(paste(tmp.arr[,'yr'],tmp.arr[,'mo'],tmp.arr[,'da'],tmp.arr[,'hr'],tmp.arr[,'fh'],sep='.'))
  out.arr <- array(as.vector(tmp.arr[,c('min','max','mean')]),dim=c(length(tmp.levs),length(tmp.times),3),
                 dimnames=list(as.character(tmp.levs),tmp.times,c('min','max','mean')))
  out.arr
}

plot.3panel.split.chkdata <- function(tmp.list=NULL,tmp.var=NULL,tmp.arr=tmp.list[[tmp.var]],
                                      psname=paste(tmp.var,'3panel','ps',sep='.'),
                                      xlab='Time Period',panels=c('min','max','mean'),
                                      ylabs=paste(tmp.var,panels),...) {
  if (is.null(tmp.arr)) stop ('Missing args: either tmp.list and tmp.var, or tmp.arr are needed\n')
  if (!is.null(psname)) {
    postscript(psname,horizontal=F)
    par(mfrow=c(length(panels),1))
  }
  for (i in 1:length(panels)) {
    x <- panels[i]
    plot(1:dim(tmp.arr)[2],tmp.arr['0',,x],xlab='Time period',ylab=ylabs[i],...)
  }
  if (!is.null(psname)) dev.off()
}

plot.6panel.split.chkdata <- function(tmp.list=NULL,tmp.var=NULL,tmp.arr=tmp.list[[tmp.var]],
                                      levs=1:dim(tmp.arr)[1],
                                      psname=paste(tmp.var,'6panel','ps',sep='.'),
                                      xlab='Time Period',panels=c('min','max','mean'),
                                      ylabs=paste(tmp.var,panels),zlab='Level',...) {
  if (is.null(tmp.arr)) stop ('Missing args: either tmp.list and tmp.var, or tmp.arr are needed\n')
  if (!is.null(psname)) {
    postscript(psname,horizontal=F)
    par(mfrow=c(length(panels),2))
  }
  for (i in 1:length(panels)) {
    x <- panels[i]
    plot(1:dim(tmp.arr)[2],apply(tmp.arr[levs,,x],MAR=2,FUN=x),xlab='Time period',ylab=ylabs[i],...)
    plot(apply(tmp.arr[levs,,x],MAR=1,FUN=x),as.numeric(dimnames(tmp.arr)[[1]])[levs],
         ylab=zlab,xlab=ylabs[i],...)
  }
  if (!is.null(psname)) dev.off()
}

