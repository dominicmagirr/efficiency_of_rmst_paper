#############################
# rmst1 (one-arm) -- hidden
#############################
rmst1=function(time, status, tau, alpha=0.05){
  #-- time
  #-- statuts
  #-- tau -- truncation time
  #-- alpha -- gives (1-alpha) confidence interval
  
  ft= survfit(Surv(time, status)~1)
  idx=ft$time<=tau
  
  wk.time=sort(c(ft$time[idx],tau))
  wk.surv=ft$surv[idx]
  wk.n.risk =ft$n.risk[idx]
  wk.n.event=ft$n.event[idx]
  
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  
  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var =c(wk.var,0)
  rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)
  
  #--- check ---
  # print(ft, rmean=tau)
  
  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se, rmst+qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL")
  colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
  
  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  Z$fit=ft
  class(Z)="rmst1"
  
  return(Z)
  
}

#########################################
# rmst2 (2-arm) contrast (main function)
#########################################
#rmst2=function(time, status, arm, tau=NULL, covariates=NULL, adjust.method="reg", alpha=0.05){
rmst2=function(time, status, arm, tau=NULL, covariates=NULL,                      alpha=0.05){
  #-- time
  #-- statuts
  #-- arm (1 or 0)
  #-- covariates (matrix)
  #-- adjust = "reg"-- regression ("aug" -- augumentation)
  #-- alpha=0.05
  
  #==================================
  #  initial check
  #==================================
  
  #===== tau =====
  idx=arm==0; tt=time[idx]; tt0max=max(tt); ss=status[idx]; ss0max=min(ss[tt==tt0max]);
  idx=arm==1; tt=time[idx]; tt1max=max(tt); ss=status[idx]; ss1max=min(ss[tt==tt1max]);
  
  ttmax = max(tt0max, tt1max)
  ttmin = min(tt0max, tt1max)
  
  #--case 1: the last obs (smaller one)=event, the last obs (longer one)=event
  if(ss0max==1 & ss1max==1){
    if(!is.null(tau)){
      if(tau>ttmax){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmax, digits=2)))}
      if(tau<=ttmax){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmax
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmax, digits=2)," is used."))
    }
  }
  
  #--case 2: the last obs (smaller one)=event, the last obs (longer one)=censor
  if((ss0max==0 & ss1max==1 & tt0max>=tt1max) | (ss0max==1 & ss1max==0 & tt1max>tt0max)){
    if(!is.null(tau)){
      if(tau>ttmax){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmax, digits=2)))}
      if(tau<=ttmax){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmax
      NOTE=paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmax, digits=2)," is used.")
    }
  }
  
  #--case 3: the last obs (smaller one)=censor, the last obs (longer one)=event
  if((ss0max==1 & ss1max==0 & tt0max>=tt1max) | (ss0max==0 & ss1max==1 & tt1max>tt0max)){
    if(!is.null(tau)){
      if(tau>ttmin){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmin, digits=2)))}
      if(tau<=ttmin){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmin
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmin, digits=2)," is used."))
    }
  }
  
  #--case 4: the last obs (smaller one)=censor, the last obs (longer one)=censor
  if(ss0max==0 & ss1max==0){
    if(!is.null(tau)){
      if(tau<=ttmin){
        NOTE=paste("The truncation time: tau =", tau, " was specified.")
      }
      if(tau>ttmin){
        stop(paste("The truncation time, tau, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(ttmin, digits=2)))
      }
    }else{
      tau = ttmin
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmin, digits=2)," is used."))
    }
  }
  
  
  Z=list()
  Z$tau=tau
  Z$note=NOTE
  
  #==================================
  #  unadjusted analysis
  #==================================
  if(is.null(covariates)){
    
    wk1=rmst1(time[arm==1], status[arm==1], tau, alpha)
    wk0=rmst1(time[arm==0], status[arm==0], tau, alpha)
    
    Z$RMST.arm1=wk1
    Z$RMST.arm0=wk0
    
    
    #--- contrast (RMST difference) ---
    rmst.diff.10     = wk1$rmst[1]-wk0$rmst[1]
    rmst.diff.10.se  = sqrt(wk1$rmst.var + wk0$rmst.var)
    rmst.diff.10.low = rmst.diff.10 - qnorm(1-alpha/2)*rmst.diff.10.se
    rmst.diff.10.upp = rmst.diff.10 + qnorm(1-alpha/2)*rmst.diff.10.se
    rmst.diff.pval   = pnorm(-abs(rmst.diff.10)/rmst.diff.10.se)*2
    rmst.diff.result = c(rmst.diff.10, rmst.diff.10.low, rmst.diff.10.upp, rmst.diff.pval)
    
    #--- contrast (RMST ratio) ---
    rmst.log.ratio.10     = log(wk1$rmst[1]) - log(wk0$rmst[1])
    rmst.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmst[1]/wk1$rmst[1] + wk0$rmst.var/wk0$rmst[1]/wk0$rmst[1])
    rmst.log.ratio.10.low = rmst.log.ratio.10 - qnorm(1-alpha/2)*rmst.log.ratio.10.se
    rmst.log.ratio.10.upp = rmst.log.ratio.10 + qnorm(1-alpha/2)*rmst.log.ratio.10.se
    rmst.log.ratio.pval   = pnorm(-abs(rmst.log.ratio.10)/rmst.log.ratio.10.se)*2
    rmst.ratio.result     = c(exp(rmst.log.ratio.10), exp(rmst.log.ratio.10.low), exp(rmst.log.ratio.10.upp),rmst.log.ratio.pval)
    
    #--- contrast (RMTL ratio  0/1) ---
    # rmtl.log.ratio.01     = log(wk0$rmtl[1]) - log(wk1$rmtl[1])
    # rmtl.log.ratio.01.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
    # rmtl.log.ratio.01.low = rmtl.log.ratio.01 - qnorm(1-alpha/2)*rmtl.log.ratio.01.se
    # rmtl.log.ratio.01.upp = rmtl.log.ratio.01 + qnorm(1-alpha/2)*rmtl.log.ratio.01.se
    # rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.01)/rmtl.log.ratio.01.se)*2
    # rmtl.ratio.result     = c(exp(rmtl.log.ratio.01), exp(rmtl.log.ratio.01.low), exp(rmtl.log.ratio.01.upp),rmtl.log.ratio.pval)
    
    #--- contrast (RMTL ratio  1/0) ---
    rmtl.log.ratio.10     = log(wk1$rmtl[1]) - log(wk0$rmtl[1])
    rmtl.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
    rmtl.log.ratio.10.low = rmtl.log.ratio.10 - qnorm(1-alpha/2)*rmtl.log.ratio.10.se
    rmtl.log.ratio.10.upp = rmtl.log.ratio.10 + qnorm(1-alpha/2)*rmtl.log.ratio.10.se
    rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.10)/rmtl.log.ratio.10.se)*2
    rmtl.ratio.result     = c(exp(rmtl.log.ratio.10), exp(rmtl.log.ratio.10.low), exp(rmtl.log.ratio.10.upp),rmtl.log.ratio.pval)
    
    
    
    #--- results ---
    out=rbind(rmst.diff.result, rmst.ratio.result , rmtl.ratio.result )
    rownames(out)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)")
    colnames(out)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")
    
    
    #--- output ---
    Z$unadjusted.result = out
    # Z$RMST.difference=out[1,]
    # Z$RMST.ratio=out[2,]
    # Z$RMTL.ratio=out[3,]
    # Z
  }
  
  #==================================
  #  Adjusted analysis
  #==================================
  if (!is.null(covariates)){
    
    ## if (adjust.method=="reg"){
    
    aa=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="difference", alpha=alpha)
    bb=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="ratio", alpha=alpha)
    cc=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="lossratio", alpha=alpha)
    
    #--- output ---
    out.adj=rbind(aa[2,c(1,5,6,4)], bb[2, c(5,6,7,4)], cc[2, c(5,6,7,4)])
    rownames(out.adj)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)")
    colnames(out.adj)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")
    
    
    #--- output ---
    Z$adjusted.result = out.adj
    
    Z$RMST.difference.adjusted = aa
    Z$RMST.ratio.adjusted      = bb
    Z$RMTL.ratio.adjusted      = cc
    
    ## }else{
    ##  stop "Please sepcify adjust.method"
    ## }
    
  }
  
  
  class(Z)="rmst2"
  
  Z
  
}
