#####################################################
#  consSurv Paper Outputs                           #
#  Angus Jennings					27Jun23	                  #
#####################################################


###############################################
# SETUP       				                        #
###############################################

#### SEED/WD/LIBS/RUNVARS ####

set.seed(34556)

#update to working directory; where to export outputs
wd <- "WORKING DIRECTORY; OUTPUT"
#ConstrainNR.R
cnr <- "WORKING DIRECTORY; ConstrainNR.R"

#where plots/datasets will output
setwd(wd)

#install required packages
list.of.packages <- c("doParallel", "simsurv", "ggplot2", "eha", "rlist", "tidyverse", "ggnewscale", "RColorBrewer", "patchwork", "officer", "ggthemes", "survPen", "pracma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

#import functions from ConstrainNR.R
source(cnr)

#can select here to only run certain sections of code

# (if F will not run and will just skip)
flag_hazplot=F #plot of dgm hazards/survs

# (if F will not run and will just import from drive)
flag_datagen=F #generate data for all sims for all dgms
flag_truermstgen=F #generate large dataset for each dgm (to calc true rmsts later)
flag_hazardgen=T #fit all models and calc hazards/rmsts
flag_rmsttruth=T #calc true rmsts using large datasets

# (if F will not run and will just skip)
flag_panelplots=T #produce all pannel plots
flag_rmstplots=F #produce plots of rmst ests/biases
flag_rmsttable=T #produce table of rmst biases

#GGTHEME

tt_theme <- theme_base() +
  theme(strip.background = element_rect(fill="grey90",color="black"), 
        text = element_text(size=11), 
        plot.title = element_text(size=10),
        plot.caption = element_text(hjust = 0),
        plot.background =  element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey95"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', colour = "grey95"))


###############################################
# HAZARDS     				                        #
###############################################

#HAZARD FUNCTIONs

weib <- function(t,p1,p2){#lambda,gamma
  p1*p2*t^(p2-1)
}
lnorm <- function(t,p1,p2){#depreciated
  eha::hlnorm(t,meanlog=p1,sdlog=p2)
}

hr <- function(t,x,betas,...){
  t0 <- as.numeric(t<tstart)
  t1 <- as.numeric(t>=tstart & t<tend)
  t2 <- as.numeric(t>=tend)
  
  if(x[["trt"]]==0){
    rep(1, length(t))
  }else{
    if(tend-tstart==0){
      diff <- 0.00000001
    }else{
      diff <- tend-tstart
    }
    t0*betas[["trt"]] + t1*( (cos(pi+(t-tstart)/(diff)*pi)+1)/2*(betas[["trt"]]) + betas[["trt"]] ) + t2
  }
  
}

haz <- function(t,x,betas,...){
  params <- lookuplist[paste0(struct,lev)][[1]]
  
  bl <- do.call(ifelse(struct%in%c("exp","weib"),"weib","lnorm"), list(t=t,p1=params[1],p2=params[2]))
  
  bl * hr(t,x,betas,...) * exp(log(betas[["zk"]])*x[["zk"]]) * exp(log(betas[["zu"]])*x[["zu"]])
}

if(flag_hazplot){#not included in manuscript
  t_vec <- seq(0.01,50,length.out=1000)
  b <- c(trt=0.5,zk=1,zu=1.5)
  par(mfrow=c(3,2), mai = c(0.5, 0.5, 0.1, 0.1))
  tstart=3
  tend=5
  lookuplist <- list("explow"=c(0.15,1),
                     "weiblow"=c(0.3,0.5),
                     "lnormlow"=c(1.5,1.1),
                     "exphigh"=c(0.3,1),
                     "weibhigh"=c(0.7,0.5),
                     "lnormhigh"=c(0.8,1))
  #exp haz
  struct="exp"
  lev="low"
  plot(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b), lty=2, type="l",ylim=c(0,0.4), xlab = "time",ylab="haz")
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b))
  lev="high"
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b), lty=2,col=2)
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b),col=2)
  #exp surv
  lev="low"
  plot(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))), type="l", ylim=c(0,1), xlab = "time",ylab="surv")
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))), lty=2)
  lev="high"
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))),col=2)
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))),col=2, lty=2)
  
  #weib haz
  struct="weib"
  lev="low"
  plot(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b), type="l", lty=2,ylim=c(0,0.4), xlab = "time",ylab="haz")
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b))
  lev="high"
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b), lty=2,col=2)
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b),col=2)
  #weib surv
  lev="low"
  plot(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))), type="l", ylim=c(0,1), xlab = "time",ylab="surv")
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))), lty=2)
  lev="high"
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))), lty=2,col=2)
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))),col=2)
  
  #lnorm haz
  struct="lnorm"
  lev="low"
  plot(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b), type="l",ylim=c(0,0.4), xlab = "time",ylab="haz")
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b), lty=2)
  lev="high"
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=0,zk=0,zu=0), betas=b),col=2)
  lines(x=t_vec,y=haz(t=t_vec, x=c(trt=1,zk=0,zu=0), betas=b),col=2, lty=2)
  #lnorm surv
  lev="low"
  plot(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))), type="l", ylim=c(0,1), xlab = "time",ylab="surv")
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))), lty=2)
  lev="high"
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=1,zk=0,zu=0),betas=b)$value))),col=2, lty=2)
  lines(x=t_vec,y=exp(-unlist(lapply(t_vec,function(i)integrate(haz, 0, i, x=c(trt=0,zk=0,zu=0),betas=b)$value))),col=2)
  
  plot.new()
  legend("center",col=c(1,1,2,2), lty=c(1,2,1,2), legend=c("pbo low","trt low","pbo high","trt high"))
  
}


###############################################
# DATA GEN     				                        #
###############################################

#SIM DATA 

start_time0 <- Sys.time()
start_time <- Sys.time()

nsim=300

tstart=3
ssparams <- expand.grid(nobs=c(600, 1000), #take ~IQR from P?ron 2012 ???
                        bx=c(0.5), #marginality means attenuation of HR, hence fully conditional will be stronger than observed in a comparable trial
                        bzk=c(1),
                        bzu=c(1.5), #only unknown het atm (all adjusted for in analysis step anyway)
                        tend=c(3),#,5,10), #UPDATEHERE
                        haz=c("low","high"),
                        hazshape=c("exp","weib"))#removed lnorm

lookuplist <- list("explow"=c(0.15,1),
                   "weiblow"=c(0.3,0.5),
                   "lnormlow"=c(1.5,1.1),
                   "exphigh"=c(0.3,1),
                   "weibhigh"=c(0.7,0.5),
                   "lnormhigh"=c(0.8,1)) #chosen to match 5 yr survival in OS vs PFS for ta903

ssparams$seed <- floor(runif(nrow(ssparams), 0, 1000000))
ssparams$seed2 <- floor(runif(nrow(ssparams), 0, 1000000))

no_cores <- min(detectCores() - 1, nrow(ssparams))

if(flag_datagen){
  
  cl <- makeCluster(no_cores, type="PSOCK")
  registerDoParallel(cl)
  
  datasets <- foreach(i=1:nrow(ssparams), .packages = 'simsurv') %dopar% {
    
    tend <- ssparams[i,"tend"]
    
    set.seed(ssparams[i,"seed"])
    
    b <- c(ssparams[i,"bx"], ssparams[i,"bzk"], ssparams[i,"bzu"])
    names(b) <- c("trt", "zk", "zu")
    
    n <- ssparams[i,"nobs"]*nsim
    df <- data.frame(ptid=1:n, trt=rbinom(n, 1, 0.5), zk=rnorm(n), zu=rnorm(n))
    
    struct <- ssparams[i,"hazshape"]
    lev <- ssparams[i,"haz"]
    
    weib <- function(t,p1,p2){#lambda,gamma
      p1*p2*t^(p2-1)
    }
    lnorm <- function(t,p1,p2){
      eha::hlnorm(t,meanlog=p1,sdlog=p2)
    }
    
    sim <- simsurv(betas = b, x = df,
                   hazard = haz,
                   maxt = 100,
                   interval = c(1E-16,10000))[2:3]
    
    df$eventtime <- ceiling(sim$eventtime*365.25)/365.25
    
    df$ptid <- rep(1:ssparams[i,"nobs"], nsim)
    df$ssid <- rep(1:nsim, each=ssparams[i,"nobs"])
    df$paramid <- i
    
    df$enroltime = runif(nrow(df),0,1)
    
    df[,c("paramid","ssid","ptid","trt","zk","zu","enroltime","eventtime")]
    
  }
  stopCluster(cl)
  
  names(datasets) <- paste0("params_",1:nrow(ssparams))
  
}else{
  load("datasets.RData")
}

save(datasets, file="datasets.RData")

end_time <- Sys.time()
t_datagen <- round(end_time - start_time,3)

#TRUE RMSTS DATA

start_time <- Sys.time()

rmstparams <- distinct(ssparams[,setdiff(colnames(ssparams),c("nobs","seed","seed2"))])

if(flag_truermstgen){
  
  cl <- makeCluster(no_cores, type="PSOCK")
  registerDoParallel(cl)
  
  truermstdatasets <- foreach(i=1:nrow(rmstparams), .packages = 'simsurv') %dopar% {
    
    tend <- rmstparams[i,"tend"]
    
    set.seed(rmstparams[i,"seed2"])
    
    b <- c(rmstparams[i,"bx"], rmstparams[i,"bzk"], rmstparams[i,"bzu"])
    names(b) <- c("trt", "zk", "zu")
    
    n <- 1000000
    df <- data.frame(ptid=1:n, trt=rbinom(n, 1, 0.5), zk=rnorm(n), zu=rnorm(n))
    
    struct <- rmstparams[i,"hazshape"]
    lev <- rmstparams[i,"haz"]
    
    weib <- function(t,p1,p2){#lambda,gamma
      p1*p2*t^(p2-1)
    }
    lnorm <- function(t,p1,p2){
      eha::hlnorm(t,meanlog=p1,sdlog=p2)
    }
    
    sim <- simsurv(betas = b, x = df,
                   hazard = haz,
                   maxt = 100,
                   interval = c(1E-16,10000))[2:3]
    
    df$eventtime <- ceiling(sim$eventtime*365.25)/365.25
    
    df$ptid <- rep(1:n)
    df$ssid <- 1
    df$paramid <- i
    
    df[,c("paramid","ssid","ptid","trt","zk","zu","eventtime")]
    
  }
  stopCluster(cl)
  
  names(truermstdatasets) <- paste0("rmst_",1:nrow(rmstparams))
  
}else{
  load("truermstdatasets.RData")
}

save(truermstdatasets, file="truermstdatasets.RData")

end_time <- Sys.time()
t_truermstgen <- round(end_time - start_time,3)


###############################################
# DATA ANALYSIS				                        #
###############################################

#SIM HAZARD / RMST CALCULATION

start_time <- Sys.time()

tt <- sort(c(3,5, seq(0.5,10, length.out=60), seq(10,40, length.out=89)[-1]))

ttdf <- data.frame(eventtime=tt, trt=rep(c(0,1), each=length(tt)), zk=0, zu=0)

modelparams <- expand.grid(bldf=c(2),
                           trtdf=c(3,4,5,6),
                           cons=c(5,10,20),
                           k95=c(F,T), t1=c(T), l95=c(0.95),#dont vary l95 for k95 F
                           form=c("cond"),
                           cens=c(3.1)) #20,3,3.1)) #UPDATEHERE
modelparams<-modelparams[!(modelparams$cens==3 & modelparams$trtdf==6) & !(modelparams$cens==3.1 & modelparams$trtdf==6) & !(modelparams$cens==20 & modelparams$trtdf==3),]

rownames(modelparams) <- 1:nrow(modelparams)

zero <- 0.00000000000001

if(flag_hazardgen){
  
  cl <- makeCluster(no_cores, type="PSOCK", outfile="")
  registerDoParallel(cl)
  
  analyses <- foreach(i=1:length(datasets), .packages = c("survPen","rlist","pracma","tidyverse")) %dopar% {
    source("~/Angus Jennings/4. U o Leicester/Paper 2_V1.1/Paper 2 - ConstrainNR_V1.1 5aug2024.R")
    
    df <- datasets[[i]]
    
    if(!(ssparams[i,"tend"]==3)){
      modelparams_sub <- modelparams[!(floor(modelparams$cens)==3),]
    }else{
      modelparams_sub <- modelparams
    }
    
    if(nrow(modelparams_sub)>0){#deleteme
    
    hazs <- lapply(1:nsim, FUN=function(j){
      lkuptbl <- c("cond"=as.formula("~ RevResCubicSpline(log(eventtime), k=k_bl) + trt + zk + zu + trt:RevResCubicSpline(log(eventtime), k=k_trt)"), 
                   "part"=as.formula("~ RevResCubicSpline(log(eventtime), k=k_bl) + trt + zk + trt:RevResCubicSpline(log(eventtime), k=k_trt)"), 
                   "marg"=as.formula("~ RevResCubicSpline(log(eventtime), k=k_bl) + trt + trt:RevResCubicSpline(log(eventtime), k=k_trt)"))
      lkuptbl_n <- c("cond"=2, 
                     "part"=1, 
                     "marg"=0)
      
      hs <- lapply(1:nrow(modelparams_sub), FUN=function(k){
        d <- df[df$ssid==j,]
        d$dead <- d$eventtime<modelparams_sub[k,"cens"]
        d$eventtime[!(d$dead)] <- modelparams_sub[k,"cens"]
        
        if(modelparams_sub[k,"cens"]==3.1){
          d$calcens <- (1-d$enroltime)+2 #time from enrolment to start of 2y clock
          d$dead <- d$eventtime<d$calcens
          d$eventtime[!(d$dead)] <- d$calcens[!(d$dead)]
        }
        
        c <- if(modelparams_sub[k,"cons"]=="NULL"){1000000}else{as.numeric(as.character(modelparams_sub[k,"cons"]))}
        k_ext <- if(modelparams_sub[k,"k95"]){modelparams_sub[k,"l95"]}else{NULL}
        
        if(modelparams_sub[k,"form"]=="cond-id"){
          .GlobalEnv$k_bl <- quantile(d$eventtime[d$dead], seq(0,1,length.out=modelparams_sub[k,"bldf"]+1) )
          .GlobalEnv$k_trt <- quantile(d$eventtime[d$dead & d$trt==1], sort(c( seq(0,1,length.out=modelparams_sub[k,"trtdf"]+1) , k_ext)))
        }else{
          .GlobalEnv$k_bl <- log(quantile(d$eventtime[d$dead], seq(0,1,length.out=modelparams_sub[k,"bldf"]+1) ))
          .GlobalEnv$k_trt <- log(quantile(d$eventtime[d$dead & d$trt==1], sort(c( seq(0,1,length.out=modelparams_sub[k,"trtdf"]+1) , k_ext))))
        }
        
        f <- lkuptbl[[modelparams_sub[k,"form"]]]
        
        sP_uncons <- tryCatch({ survPen_cons(f,data=d,t1=eventtime,event=dead, max.it.beta = 5000) }, error=function(e){NA} )
        
        if(all(is.na(sP_uncons))){
          
          h <- data.frame(rep(NA,nrow(ttdf)))
          
          sSP_pbocons <- list(mrmst=list(trt1=data.frame(t= ifelse(modelparams_sub[k,"cens"]=="3"|modelparams_sub[k,"cens"]=="3.1",3,20) , Est=NA), trt0=data.frame(t= ifelse(modelparams_sub[k,"cens"]=="3"|modelparams_sub[k,"cens"]=="3.1",3,20) , Est=NA) ))
          
          h <- cbind(h, data.frame(rep(NA,nrow(ttdf))), c(sSP_pbocons$mrmst$trt0$Est, sSP_pbocons$mrmst$trt1$Est))
          
        }else{
          
          h <- data.frame(predict(sP_uncons, newdata=ttdf, do.surv=F)$haz)
          
          .GlobalEnv$k_trt <- log(quantile(d$eventtime[d$dead & d$eventtime<c & d$trt==1], sort(c( seq(0,1,length.out=modelparams_sub[k,"trtdf"]+1) , k_ext))))
          if(c<1000000){
            .GlobalEnv$k_trt[length(.GlobalEnv$k_trt)] <- log(c)
          }
          
          p <- 1+modelparams_sub[k,"bldf"]+1+lkuptbl_n[[modelparams_sub[k,"form"]]]+modelparams_sub[k,"trtdf"]+as.numeric(modelparams_sub[k,"k95"])
          t1 <- if(modelparams_sub[k,"t1"] & !(c==1000000)){modelparams_sub[k,"bldf"]+2}else{NULL}
          c <- if(c==1000000){NULL}else{p}
          
          uncons.keep <- 1:(modelparams_sub[k,"bldf"]+1)
          cons.zero <- c(t1, c)
          
          b.ini <- rep(0,p)
          b.ini[uncons.keep] <- sP_uncons$coefficients[uncons.keep]
          b.ini[cons.zero] <- zero
          
          sP_pbocons <- tryCatch({ survPen_cons(f,data=d,t1=eventtime,event=dead, beta.ini=b.ini, cons=c(uncons.keep, cons.zero), max.it.beta = 5000) }, error=function(e){NA} )
          
          if(all(is.na(sP_pbocons))){
            
            sSP_pbocons <- list(mrmst=list(trt1=data.frame(t= ifelse(modelparams_sub[k,"cens"]=="3"|modelparams_sub[k,"cens"]=="3.1",3,20) , Est=NA), trt0=data.frame(t= ifelse(modelparams_sub[k,"cens"]=="3"|modelparams_sub[k,"cens"]=="3.1",3,20) , Est=NA) ))
            
            h <- cbind(h, data.frame(rep(NA,nrow(ttdf))), c(sSP_pbocons$mrmst$trt0$Est, sSP_pbocons$mrmst$trt1$Est))
            
          }else{
            
            sSP_pbocons <- standSurvPen(object=sP_pbocons, data=d, at=list(trt1=list(trt=1),trt0=list(trt=0)), times=ifelse(modelparams_sub[k,"cens"]=="3"|modelparams_sub[k,"cens"]=="3.1",3,20), type=c("rmst"), n.legendre=20)
            
            h <- cbind(h,predict(sP_pbocons, newdata=ttdf, do.surv=F)$haz, c(sSP_pbocons$mrmst$trt0$Est, sSP_pbocons$mrmst$trt1$Est))
            
          }
        }
        
        names(h) <- paste0("s_",j,
                           "_",modelparams_sub[k,"cens"],
                           "_",modelparams_sub[k,"form"],
                           "_",modelparams_sub[k,"bldf"],
                           "_",modelparams_sub[k,"trtdf"],
                           "_",modelparams_sub[k,"cons"],
                           "_",modelparams_sub[k,"t1"],
                           "_",modelparams_sub[k,"k95"],
                           "_",modelparams_sub[k,"l95"],
                           "_", c("uncons","cons","rmsts"))
        
        print(paste0("analysis ",i,", modelparams ",k,", ssid ",j," completed."))
        h
        
      })
      
      list.cbind(hs)
      
    })
    
    hazs <- cbind(ttdf, list.cbind(hazs))
    
    }#deleteme
    
  }
  stopCluster(cl)
  
  names(analyses) <- paste0("params_",1:length(analyses))
  
  rmsts <- lapply(analyses,FUN=function(df){
    d <- select(df, ends_with("rmsts"))[1:4,]
    rownames(d) <- c("trt0_uncons", "trt1_uncons", "trt0_cons", "trt1_cons")
    d
  })
  analyses <- lapply(analyses,FUN=function(df){
    select(df, !(ends_with("rmsts")))
  })
  
}else{
  load("analyses.RData")
  load("rmsts.RData")
}

names(analyses) <- names(datasets)

save(analyses, file="analyses.RData")
save(rmsts, file="rmsts.RData")

end_time <- Sys.time()
t_hazardgen <- round(end_time - start_time,3)

#TRUE RMST CALCULATION

start_time <- Sys.time()

if(flag_rmsttruth){
  truermsts <- lapply(truermstdatasets, FUN=function(df){
    rmst40 <- survRM2::rmst2(df$eventtime, !(df$eventtime==100), df$trt, tau=20)
    rmst3 <- survRM2::rmst2(df$eventtime, !(df$eventtime==100), df$trt, tau=3)
    c(r40t0=unname(rmst40$RMST.arm0$rmst[1]),r40t1=unname(rmst40$RMST.arm1$rmst[1]),r3t0=unname(rmst3$RMST.arm0$rmst[1]),r3t1=unname(rmst3$RMST.arm1$rmst[1]))
    
  })
}else{
  load("truermsts.RData")
}

save(truermsts, file="truermsts.RData")

end_time <- Sys.time()
t_rmsttruth <- round(end_time - start_time,3)

end_time0 <- Sys.time()

t_totalgen <- round(end_time0 - start_time0,3)

#TOTAL TIMES

print(t_datagen)
print(t_truermstgen)
print(t_hazardgen)
print(t_rmsttruth)
print(t_totalgen)


###############################################
# HAZARD VIS (panels)	                        #
###############################################

start_time <- Sys.time()

cols <- c(brewer.pal(9,"Blues")[c(5,8)], brewer.pal(9,"Greens")[c(5,8)])

#PANELS

if(flag_panelplots){
  truehaz <- ttdf
  for(i in 1:nrow(ssparams)){
    struct=ssparams[i,"hazshape"]
    lev=ssparams[i,"haz"]
    tend <- ssparams[i,"tend"]
    truehaz[,paste0("ssparam_",i)] <-
      apply(truehaz, MARGIN = 1, FUN = function(x){haz(t=x[1], x=data.frame("trt"=x[2], "zk"=x[3], "zu"=x[4]), betas=data.frame("trt"=ssparams[i,"bx"], "zk"=ssparams[i,"bzk"], "zu"=ssparams[i,"bzu"]))})
    
  }
  rm(struct, lev, tend)
  
  panelvars <- c("trtdf","cons")
  panels <- unique(modelparams[,setdiff(colnames(modelparams),panelvars)])
  
  lapply(seq_along(analyses), FUN=function(i){
    if(!(ssparams[i,"tend"]==3)){
      panels_sub <- panels[!(floor(panels$cens)==3),]
    }else{
      panels_sub <- panels
    }
    
    analysis <- analyses[[i]]
    #corresponding true hazards
    t <- truehaz[,c("eventtime","trt",paste0("ssparam_",i))]
    colnames(t)[3] <- "true"
    #sims to rows
    analysis_long <- pivot_longer(analysis, starts_with("s_"), names_to="modelparams",values_to="Haz")
    #split/add sim name components
    analysis_long_Sim <- str_split(analysis_long$modelparams, "_", simplify=T)[,-1]
    colnames(analysis_long_Sim) <- c("sim", "cens", "form", "bldf", "trtdf", "cons", "t1", "k95", "l95", "conspbo")
    analysis_long <- cbind(analysis_long,analysis_long_Sim) %>% filter(sim %in% 1:100)
    #add truths
    analysis_long <- merge(analysis_long, t, by=c("eventtime","trt"))
    #calc avgs
    analysis_wider <- pivot_wider(select(analysis_long,-modelparams), names_from=sim, names_prefix="s_", values_from=Haz)
    analysis_wider$avg <- rowMeans(select(analysis_wider, starts_with("s_")), na.rm=T)
    #relengthen
    analysis_long <- pivot_longer(analysis_wider, starts_with("s_"),names_to="s_sim",values_to="Haz")
    analysis_long$sim <- str_split(analysis_long$s_sim, "_", simplify=T)[,-1]
    #produce panels
    lapply(1:nrow(panels_sub), FUN=function(j){
      p <- panels_sub[j,]
      #matching scenarios
      a <- analysis_long[rowSums(!(analysis_long[,colnames(p)] == p[rep(1, nrow(analysis_long)),]))==0,]
      #form no constraint col from 10y cons uncons vals
      a_uncons <- filter(a, cons=="10" & conspbo=="uncons")
      a_uncons <- transform(a_uncons, cons="NULL", conspbo="cons")
      a <- rbind(a, a_uncons)
      #caption
      cap <- paste0(nsim,"x",ssparams[i,"nobs"]," samples, ",
                    ssparams[i,"haz"],"/",ssparams[i,"hazshape"]," hazard, ",
                    "waning 3-",ssparams[i,"tend"],"years,",
                    ifelse(p["k95"],"with ","no "),p["l95"]*100,"th %ile knot.")
      if(panels_sub[j,"cens"]==20){
        xlims <- c(20,0)
        aim="Aim 1"
      }else{
        xlims <- c(10,40,0)
        aim="Aim 2"
      }
      #3 plot types
      for(type in xlims){
        if(type==0){
          #hr
          a_0 <- filter(a, trt==0)
          a_1 <- filter(a, trt==1)
          a_hr <- merge(a_0, a_1, suffixes=c(".pbo",".trt"), by=setdiff(colnames(a_0),c("trt","Haz","true","avg")),all=T)
          #calc hrs/avgs/truths
          a_hr$HR <- a_hr$Haz.trt/a_hr$Haz.pbo
          a_hr$avg <- a_hr$avg.trt/a_hr$avg.pbo
          a_hr$true <- a_hr$true.trt/a_hr$true.pbo
          #lims
          xmax <- 20
          ymin <- min(a_hr[a_hr$eventtime<=xmax,]$HR, a_hr[a_hr$eventtime<=xmax,]$true, na.rm=T)
          ymax <- min(max(a_hr[a_hr$eventtime<=xmax,]$HR, a_hr[a_hr$eventtime<=xmax,]$true, na.rm=T) , 1.5)
          #labs
          labs <- expand.grid(cons=c(5,10,20,"NULL"),cens=p$cens,label=c(" Constraint"), line=c(ymin,ymax), y=c(ymax-ymax/20))
          #labs[labs$label==" Constraint","y"] <- labs[labs$label==" Constraint","y"] - ymax/20
          labs$x <- ifelse(labs$label==" Constraint", as.numeric(as.character(labs$cons)), as.numeric(as.character(labs$cens)))
          #caption
          caption <- paste0("Constrained HR, censored from ",p$cens,"y: ",cap)
          #full plot
          g_all <- ggplot(mapping=aes(x=eventtime)) + 
            geom_line(data=labs, aes(x=as.numeric(as.character(cens)), y=line), color="gray70") +
            geom_line(data=labs, aes(x=as.numeric(as.character(cons)), y=line), linetype="dashed", color="gray50") +
            #geom_text(data=filter(labs, line==ymin), aes(x=x,y=y,label=label), color="gray70", hjust=0, size=3) +
            geom_line(data=filter(a_hr,conspbo=="cons"), aes(y=HR, group=sim),alpha=0.3, color=brewer.pal(9, "BuGn")[6]) + 
            geom_line(data=filter(a_hr,conspbo=="cons", sim==1 & eventtime<panels_sub[j,"cens"]), aes(y=true), color="gray30") + 
            geom_line(data=filter(a_hr,conspbo=="cons", sim==1), aes(y=avg), size=0.7, color=brewer.pal(9, "BuGn")[9]) +
            facet_grid(switch="y", paste0("Trt df: ", trtdf) ~ factor(factor(paste0("cons=",cons),levels=paste0("cons=",c(5,10,20,"NULL"))),labels=c(paste0("HR Constraint: ",c(5,10,20),"y"),"No constraint"))) +
            coord_cartesian(ylim=c(ymin,ymax)) + xlim(0,xmax) + labs(y="HR",x="Time (years)",title=caption, caption=paste0("Bold line corresponds to average of ",nsim," repeats (paler lines)."))+
            tt_theme
          ggsave(paste0(aim,"/HR Plots/DGM",i,"_panel",j,".png"), g_all, width = 26, height = 15, units = "cm")
          
        }else{
          #hazs (10/40y)
          #lims
          xmax <- type
          ymin <- min(a[a$eventtime<=xmax,]$Haz, a[a$eventtime<=xmax,]$true, na.rm=T)
          ymax <- min(max(a[a$eventtime<=xmax,]$Haz, a[a$eventtime<=xmax,]$true,na.rm=T) , ifelse(ssparams[i,"haz"]=="low",0.3,0.5), na.rm=T)
          #labs
          labs <- expand.grid(cons=c(5,10,20,"NULL"),cens=p$cens,label=c(" Constraint"), line=c(ymin,ymax), y=c(ymax-ymax/20))
          #labs[labs$label==" Constraint","y"] <- labs[labs$label==" Constraint","y"] - ymax/20
          labs$x <- ifelse(labs$label==" Constraint", as.numeric(as.character(labs$cons)), as.numeric(as.character(labs$cens)))
          #caption
          caption <- paste0("Constrained hazards, censored from ",p$cens,"y: ",cap)
          #full plot
          g_all <- ggplot(mapping=aes(x=eventtime)) + 
            geom_line(data=labs, aes(x=as.numeric(as.character(cens)), y=line), color="gray70") +
            geom_line(data=labs, aes(x=as.numeric(as.character(cons)), y=line), linetype="dashed", color="gray50") +
            #geom_text(data=filter(labs, line==ymin), aes(x=x,y=y,label=label), color="gray70", hjust=0) +
            geom_line(data=filter(a,conspbo=="cons"), aes(y=Haz,group=paste0(sim,trt), color=factor(trt)),alpha=0.2) +
            geom_line(data=filter(a,conspbo=="cons", sim==1 & eventtime<panels_sub[j,"cens"]), aes(y=true,group=trt), color="gray30") + 
            scale_color_manual(values = cols[c(1,3)]) +
            guides(color = "none") +
            new_scale_colour() +
            geom_line(data=filter(a,conspbo=="cons", sim==1), aes(y=avg,color=factor(trt)), size=0.7) +
            facet_grid(switch="y", paste0("Trt df: ", trtdf) ~ factor(factor(paste0("cons=",cons),levels=paste0("cons=",c(5,10,20,"NULL"))),labels=c(paste0("HR Constraint: ",c(5,10,20),"y"),"No constraint"))) +
            coord_cartesian(ylim=c(ymin,ymax)) + xlim(0,xmax) + labs(y="Hazards",x="Time (years)",title=caption, color="Treatment", caption=paste0("Bold line corresponds to average of ",nsim," repeats (paler lines).")) +
            scale_color_manual(values = cols[c(2,4)]) + tt_theme
          ggsave(paste0(aim,"/",type,"y Haz Plots/DGM",i,"_panel",j,".png"), g_all, width = 26, height = 17, units = "cm")
          
        }
        
      }
      
    })
    print(paste0(nrow(panels_sub)," panels for DGM ",i," completed."))
  })
}

end_time <- Sys.time()
t_panelplots <- round(end_time - start_time,3)


###############################################
# RMST VIS (RMST plots/tables)                #
###############################################

start_time <- Sys.time()

# FORMAT RMST DATA
# TABLULATE CONVERGENCE RATES

if(flag_rmstplots | flag_rmsttable){
  #reformat rmsts dataset
  load("rmsts.RData")
  rmsts <- lapply(seq_along(rmsts), FUN=function(i){
    df <- rmsts[[i]]
    d <- rownames_to_column(as.data.frame(t(df)))
    d_rowname <- str_split(d$rowname, "_", simplify = T)[,-c(1,11)]
    colnames(d_rowname) <- c("sim","cens","form","bldf","trtdf","cons","t1","k95","l95")
    d <- cbind(d, d_rowname)
    d$rowname <- apply(d_rowname[,-1], MARGIN = 1, FUN=function(x){paste0(x, collapse = "_")})
    d$true0 <- ifelse(d$cens=="3"|d$cens=="3.1",truermsts[[floor((i+1)/2)]][3],truermsts[[floor((i+1)/2)]][1]) #duplicate i's for the nobs, as long as this alternates 600/1000/600/1000etc
    d$true1 <- ifelse(d$cens=="3"|d$cens=="3.1",truermsts[[floor((i+1)/2)]][4],truermsts[[floor((i+1)/2)]][2])
    d
  })
  #calc avg rmsts
  rmsts_avg <- lapply(rmsts, FUN=function(df){
    df %>% group_by(rowname) %>% summarise(trt0_uncons=mean(trt0_uncons, na.rm=T),
                                           trt1_uncons=mean(trt1_uncons, na.rm=T),
                                           trt0_cons=mean(trt0_cons, na.rm=T),
                                           trt1_cons=mean(trt1_cons, na.rm=T)) })
  
  # CONVERGENCE
  
  #convergence by dgm, by model
  missing <- lapply(seq_along(rmsts), function(i){
    df <- rmsts[[i]]
    df %>% group_by(rowname) %>% summarise(n=n(),missing=sum(is.na(trt0_uncons))) %>% as.data.frame() %>% mutate(dgm=i)
  })
  #limit to only cases with non-convergence
  missingN0 <- lapply(missing, function(df){
    df[df$missing>0,]
  }) %>% list.rbind()
  #always with the 95th percentile knot
  #all are 20y models
  #all constrained at 5yrs
  #most using 6df except for 2 on 5df
  
  #covergence by dgm
  missing_bydgm <- lapply(missing, function(df){
    d <- colSums(df[,-1])
    paste0(d["missing"]," / ",d["n"]," (",round(d["missing"]*100/d["n"],1),")")
  }) %>% list.rbind()
  missing_bydgm <- ssparams %>% mutate(missing=missing_bydgm)
  #vast majority in low exp 3y end
  #more missing in exp
  #more missing with a faster waning
  #only 14 occurred in the non-3y waning (/172800 = 0.008%)
  
  #total convergence rates
  missing_total <- missing %>% list.rbind %>% .[,2:3] %>% colSums
  #0.15% not converged (/172800)
  
}

# RMST PLOTS

if(flag_rmstplots){
  
  cols <- c(brewer.pal(9,"BuPu")[c(6,8)], brewer.pal(9,"OrRd")[c(5,7)])
  
  #calc mins/maxs
  minmax <- list.rbind(rmsts) %>% 
    pivot_longer(starts_with("trt") & ends_with("cons")) %>% 
    summarise(min=min(value, na.rm=T),
              max=max(value, na.rm=T))
  minmax_bias <- list.rbind(rmsts) %>% 
    pivot_longer(starts_with("trt") & ends_with("cons"))
  minmax_bias$bias <- minmax_bias$value - ifelse(str_starts(minmax_bias$name,"trt0"),minmax_bias$true0,minmax_bias$true1)
  minmax_bias <- minmax_bias %>% 
    summarise(min=min(bias, na.rm=T),
              max=max(bias, na.rm=T))
  
  lapply(seq_along(rmsts), FUN=function(i){
    if(!(ssparams[i,"tend"]==3)){
      modelparams_sub <- modelparams[!(modelparams$cens==3|modelparams$cens==3.1),]
    }else{
      modelparams_sub <- modelparams
    }
    
    rmst <- rmsts[[i]]
    rmst <- pivot_longer(rmst, starts_with("trt0")|starts_with("trt1"), values_to = "rmst", names_to = "type")
    rmst_avg <- rmsts_avg[[i]]
    rmst_avg <- pivot_longer(rmst_avg, starts_with("trt"), values_to = "rmst", names_to = "type")
    
    rn_order <- filter(rmst,type=="trt0_uncons"&sim==1)$rowname
    rn_order <- rn_order[order(str_split(rn_order,"_",simplify=T)[,1])]
    
    rmst$rowname <- factor(rmst$rowname, levels=rn_order)
    rmst_avg$rowname <- factor(rmst_avg$rowname, levels=rn_order)
    
    theme_ssparams <- theme(axis.title.y=element_text(angle=0,vjust=0.5), axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.caption = element_text(hjust=0))
    
    l <- length(levels(factor(modelparams_sub$cens)))*length(levels(factor(modelparams_sub$form)))
    
    #l1 = levels(factor(modelparams_sub$cens))[1]; l2 = levels(factor(modelparams_sub$form))[1]
    for(j in 0:1){
      for(l1 in levels(factor(modelparams_sub$cens))){#"3" or "20"
        for(l2 in levels(factor(modelparams_sub$form))){#depreciated
          if(l1=="20"){aim="Aim 1"}else{aim="Aim 2"}
          
          r <- filter(rmst, str_starts(type,paste0("trt",j)) & cens==l1 & form==l2)
          r_avg <- filter(rmst_avg, str_starts(type,paste0("trt",j)) & str_starts(rowname,paste(l1,l2,sep="_")))
          
          ss <- as.data.frame(str_split(levels(factor(r$rowname)),"_",simplify=T))
          ss <- ss[c(1:nrow(ss),nrow(ss)),]
          colnames(ss) <- c("cens","form","bldf","trtdf","cons","t1","k95","l95")
          
          g_trtdf <- ggplot() + geom_step(data=ss, mapping=aes(x=0:(nrow(ss)-1),y=as.numeric(trtdf)), size=0.8, color="gray25") +
            scale_x_continuous(limits=c(0,(nrow(ss)-1)), expand=c(0,0)) + scale_y_continuous(expand=c(0.1,0.1), breaks = as.numeric(levels(factor(ss$trtdf)))) + labs(y="Trt df",title="Simulation Parameters:") +
            theme_ssparams
          g_cons <- ggplot() + geom_step(data=ss, mapping=aes(x=0:(nrow(ss)-1),y=as.numeric(cons)), size=0.8, color="gray25") +
            scale_x_continuous(limits=c(0,(nrow(ss)-1)), expand=c(0,0)) + scale_y_continuous(expand=c(0.1,0.1), breaks = sort(as.numeric(levels(factor(ss$cons)))), labels=c(paste0(sort(as.numeric(levels(factor(ss$cons)))),"y"))) + labs(y="Constrained") +
            theme_ssparams
          g_k95 <- ggplot() + geom_step(data=ss, mapping=aes(x=0:(nrow(ss)-1),y=as.numeric(as.logical(k95))), size=0.8, color="gray25") +
            scale_x_continuous(limits=c(0,(nrow(ss)-1)), expand=c(0,0)) + scale_y_continuous(limits=c(0,1), expand=c(0.1,0.1), breaks = c(0,1), labels=c("No","Yes")) + labs(y="95th Knot?") +
            theme_ssparams
          
          cap <- paste0(ifelse(j==0,"PBO ","TRT "),ifelse(l1=="3"|l1=="3.1",3,20),"y RMSTs (censored from ",l1,"y): ","TRT RMSTs: ",nsim,"x",ssparams[i,"nobs"]," samples, ",
                        ssparams[i,"haz"],"/",ssparams[i,"hazshape"]," hazard, ",
                        "waning 3-",ssparams[i,"tend"],"years.")
          
          g_jitter <- ggplot() + 
            geom_point(data=r,mapping=aes(x=as.numeric(sim),y=rmst,color=type), alpha=0.4) + 
            scale_color_manual(values = cols[c(1,3)]) +
            guides(color = "none") +
            new_scale_colour() + 
            geom_hline(yintercept = ifelse(l1=="3"|l1=="3.1",truermsts[[floor((i+1)/2)]][3+j],truermsts[[floor((i+1)/2)]][1+j]), size=0.8, color="gray25") + #UPDATE
            geom_hline(data=r_avg,mapping=aes(yintercept=rmst,color=factor(type,labels=c("Cons.","Uncons."))), size=0.8) +
            facet_grid(.~factor(rowname,labels=paste0("Scen. ",1:(nrow(modelparams_sub)/l)))) + labs(color="", y="RMSTs", title=cap, caption="Marginal 40y RMST (reg. standardisation) from each simulation with/without HR constrained to 1. Average of estimates shown by coloured line. \nTrue RMST (from numerical int. of KM of sample size 1,000,000) indicated by black line.") +
            scale_color_manual(values = cols[c(2,4)]) +
            theme_ssparams
          
          g <- g_jitter / g_trtdf / g_cons / g_k95  + 
            plot_layout(heights = c(8, 1,1,1))
          
          ggsave(paste0(aim,"/RMST Plots/dynamicy/DGM",i,"_trt",j,".png"), g, width = 26, height = 17, units = "cm")
          
          g <- (g_jitter+scale_y_continuous(limits=unlist(as.list(minmax)),expand=c(0,0))) / g_trtdf / g_cons / g_k95 + 
            plot_layout(heights = c(8, 1,1,1))
          
          ggsave(paste0(aim,"/RMST Plots/fixy/DGM",i,"_trt",j,".png"), g, width = 26, height = 17, units = "cm")
          
          g_jitter <- ggplot() + 
            geom_point(data=r,mapping=aes(x=as.numeric(sim),y=rmst-ifelse(l1=="3"|l1=="3.1",truermsts[[floor((i+1)/2)]][3+j],truermsts[[floor((i+1)/2)]][1+j]),color=type), alpha=0.4) + 
            scale_color_manual(values = cols[c(1,3)]) +
            guides(color = "none") +
            new_scale_colour() + 
            geom_hline(yintercept = 0, size=0.8, color="gray25") +
            geom_hline(data=r_avg,mapping=aes(yintercept=rmst-ifelse(l1=="3"|l1=="3.1",truermsts[[floor((i+1)/2)]][3+j],truermsts[[floor((i+1)/2)]][1+j]),color=factor(type,labels=c("Cons.","Uncons."))), size=0.8) +
            facet_grid(.~factor(rowname,labels=paste0("Scen. ",1:(nrow(modelparams_sub)/l)))) + labs(color="", y="RMST\nBias", title=cap, caption="Marginal 40y RMST (reg. standardisation) from each simulation with/without HR constrained to 1. Average of estimates shown by coloured line. \nTrue RMST (from numerical int. of KM of sample size 1,000,000) indicated by black line.") +
            scale_color_manual(values = cols[c(2,4)]) +
            theme_ssparams
          
          g <- (g_jitter+scale_y_continuous(limits=unlist(as.list(minmax_bias)),expand=c(0,0))) / g_trtdf / g_cons / g_k95 + 
            plot_layout(heights = c(8, 1,1,1))
          
          ggsave(paste0(aim,"/RMST Plots/bias/DGM",i,"_trt",j,".png"), g, width = 26, height = 17, units = "cm")
          
          
        }
      }
    }
    
  })
  
}

end_time <- Sys.time()
t_rmstplots <- round(end_time - start_time,3)

# RMST TABLE

start_time <- Sys.time()

if(flag_rmsttable){
  #avg rmsts
  rmsts_avg_all <- lapply(seq_along(rmsts_avg), FUN=function(i){
    df <- rmsts_avg[[i]]
    d <- cbind(str_split(df$rowname,"_",simplify=T),df)
    colnames(d)[1:8] <- c("cens", "form", "bldf", "trtdf", "cons", "t1", "k95", "l95")
    d <- cbind(ssparams[i,],d)
    d
  }) %>% list.rbind() %>% mutate(rmsts = paste0(round(trt0_cons,2)," / ",round(trt1_cons,2)," (",round(round(trt1_cons,2)-round(trt0_cons,2),2),")"),
                                 cols = paste0("c_",k95,"_", cons,"_", trtdf)) %>% 
    select(-rowname,-seed,-seed2, -(trt0_uncons:trt1_cons), -k95, -cons, -trtdf)
  
  #append true rmsts to ssparams
  ssparams$rmst40.0 <- rep(sapply(truermsts,FUN=function(tr) tr[1]),each=2)#double for the 2nobs
  ssparams$rmst40.1 <- rep(sapply(truermsts,FUN=function(tr) tr[2]),each=2)
  ssparams$rmst3.0 <- rep(sapply(truermsts,FUN=function(tr) tr[3]),each=2)
  ssparams$rmst3.1 <- rep(sapply(truermsts,FUN=function(tr) tr[4]),each=2)
  
  #rmst biases PBO / TRT (Delta)
  rmsts_bias_all <- lapply(seq_along(rmsts_avg), FUN=function(i){
    df <- rmsts_avg[[i]]
    d <- cbind(str_split(df$rowname,"_",simplify=T),df)
    colnames(d)[1:8] <- c("cens", "form", "bldf", "trtdf", "cons", "t1", "k95", "l95")
    d <- cbind(ssparams[i,],d)
    d
  }) %>% list.rbind() %>% mutate(rmst0=ifelse(cens=="3"|cens=="3.1",rmst3.0,rmst40.0), rmst1=ifelse(cens=="3"|cens=="3.1",rmst3.1,rmst40.1)) %>% 
    mutate(rmsts = paste0(round(trt0_cons-rmst0,2)," / ",round(trt1_cons-rmst1,2)," (",round( round(trt1_cons-trt0_cons,2)-round(rmst1-rmst0,2) ,2),")"),
           cols = paste0("c_",k95,"_", cons,"_", trtdf)) %>% 
    select(-rowname,-seed,-seed2, -(trt0_uncons:trt1_cons), -k95, -cons, -trtdf, -rmst0, -rmst1)
  
  #rmst Delta bias
  rmsts_deltabias_all <- lapply(seq_along(rmsts_avg), FUN=function(i){
    df <- rmsts_avg[[i]]
    d <- cbind(str_split(df$rowname,"_",simplify=T),df)
    colnames(d)[1:8] <- c("cens", "form", "bldf", "trtdf", "cons", "t1", "k95", "l95")
    d <- cbind(ssparams[i,],d)
    d
  }) %>% list.rbind() %>% mutate(rmst0=ifelse(cens=="3"|cens=="3.1",rmst3.0,rmst40.0), rmst1=ifelse(cens=="3"|cens=="3.1",rmst3.1,rmst40.1)) %>% 
    mutate(rmsts = round( (trt1_cons-trt0_cons) - (rmst1-rmst0) ,3),
           cols = paste0("c_",k95,"_", cons,"_", trtdf)) %>% 
    select(-rowname,-seed,-seed2, -(trt0_uncons:trt1_cons), -k95, -cons, -trtdf, -rmst0, -rmst1)
  
  #RMST Delta % bias
  rmsts_deltaPERCbias_all <- lapply(seq_along(rmsts_avg), FUN=function(i){
    df <- rmsts_avg[[i]]
    d <- cbind(str_split(df$rowname,"_",simplify=T),df)
    colnames(d)[1:8] <- c("cens", "form", "bldf", "trtdf", "cons", "t1", "k95", "l95")
    d <- cbind(ssparams[i,],d)
    d
  }) %>% list.rbind() %>% mutate(rmst0=ifelse(cens=="3"|cens=="3.1",rmst3.0,rmst40.0), rmst1=ifelse(cens=="3"|cens=="3.1",rmst3.1,rmst40.1)) %>% 
    mutate(rmsts = round( (((trt1_cons-trt0_cons) - (rmst1-rmst0))/(rmst1-rmst0))*100 ,1),
           cols = paste0("c_",k95,"_", cons,"_", trtdf)) %>% 
    select(-rowname,-seed,-seed2, -(trt0_uncons:trt1_cons), -k95, -cons, -trtdf, -rmst0, -rmst1)
  
  #c = levels(factor(rmsts_avg_all$cens))[1];f = levels(factor(rmsts_avg_all$form))[1]
  for(c in levels(factor(rmsts_avg_all$cens))){
    for(f in levels(factor(rmsts_avg_all$form))){
      cols_order <- paste0("c_",apply(expand.grid(c(3,4,5)+ifelse(c=="20",1,0),c(5,10,20),c("FALSE","TRUE"))[,3:1], MARGIN = 1, FUN=paste0, collapse="_" )) %>% str_remove_all(" ")
      
      rmst_avg_all <- filter(rmsts_avg_all, cens==c & form==f) %>% 
        select(-cens, -form,-bx, -bzk, -bzu, -bldf, -t1, -l95)
      rmst_tbl <- pivot_wider(rmst_avg_all, names_from=cols, values_from = rmsts)
      
      rmst_bias_all <- filter(rmsts_bias_all, cens==c & form==f) %>% 
        select(-cens, -form,-bx, -bzk, -bzu, -bldf, -t1, -l95)
      rmst_tbl_bias <- pivot_wider(rmst_bias_all, names_from=cols, values_from = rmsts)
      
      rmst_deltabias_all <- filter(rmsts_deltabias_all, cens==c & form==f) %>% 
        select(-cens, -form,-bx, -bzk, -bzu, -bldf, -t1, -l95)
      rmst_tbl_deltabias <- pivot_wider(rmst_deltabias_all, names_from=cols, values_from = rmsts)
      
      rmst_deltaPERCbias_all <- filter(rmsts_deltaPERCbias_all, cens==c & form==f) %>% 
        select(-cens, -form,-bx, -bzk, -bzu, -bldf, -t1, -l95)
      rmst_tbl_deltaPERCbias <- pivot_wider(rmst_deltaPERCbias_all, names_from=cols, values_from = rmsts)
      
      #format tables in the right order for headings etc
      if(c %in% c("3","3.1")){
        aim <- "Aim 2"
        rows_order <- paste0(ssparams$haz,"_",ssparams$hazshape) %>% unique
        rmst_tbl <- select(rmst_tbl,-tend) %>% 
          mutate(rows = paste0(haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_bias <- select(rmst_tbl_bias,-tend) %>% 
          mutate(rows = paste0(haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_deltabias <- select(rmst_tbl_deltabias,-tend) %>% 
          mutate(rows = paste0(haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_deltaPERCbias <- select(rmst_tbl_deltaPERCbias,-tend) %>% 
          mutate(rows = paste0(haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
      }else{
        aim <- "Aim 1"
        rows_order <- paste0(ssparams$tend,"_",ssparams$haz,"_",ssparams$hazshape) %>% unique
        rmst_tbl <- mutate(rmst_tbl, rows = paste0(tend,"_",haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_bias <- mutate(rmst_tbl_bias, rows = paste0(tend,"_",haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_deltabias <- mutate(rmst_tbl_deltabias, rows = paste0(tend,"_",haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
        rmst_tbl_deltaPERCbias <- mutate(rmst_tbl_deltaPERCbias, rows = paste0(tend,"_",haz,"_",hazshape)) %>% 
          mutate(rows=factor(rows, levels=rows_order)) %>% 
          arrange(rows, nobs) %>% 
          select(nobs, rows, starts_with("c"))
      }
      
      #reorder columns
      rmst_tbl <- rmst_tbl[,c("nobs","rows",cols_order)]
      rmst_tbl_bias <- rmst_tbl_bias[,c("nobs","rows",cols_order)]
      rmst_tbl_deltabias <- rmst_tbl_deltabias[,c("nobs","rows",cols_order)]
      rmst_tbl_deltaPERCbias <- rmst_tbl_deltaPERCbias[,c("nobs","rows",cols_order)]
      
      #output mean rmst table
      read_docx() %>%
        body_add(paste0(ifelse(c=="3"|c=="3.1",3,20),"y Average RMST: PBO / TRT (Delta)"), style = "heading 1") %>%
        body_add_par(paste0(nsim,"x600 samples, ",
                            "on ", ifelse(c=="3"|c=="3.1",3,20),"y censored data."), style = "Normal") %>%
        body_add_table(rmst_tbl[,c(1:2,3:11)], style = "table_template") %>%
        body_add(paste0(ifelse(c=="3"|c=="3.1",3,20),"y Average RMST: PBO / TRT (Delta)"), style = "heading 1") %>%
        body_add_par(paste0(nsim,"x1000 samples, ",
                            "on ", ifelse(c=="3"|c=="3.1",3,20),"y censored data."), style = "Normal") %>%
        body_add_table(rmst_tbl[,c(1:2,12:20)], style = "table_template") %>%
        print(target =paste0(aim,"/RMST Tables/RMST.docx"))
      #output mean rmst bias table
      read_docx() %>%
        body_add(paste0(ifelse(c=="3"|c=="3.1",3,20),"y Average RMST Bias: PBO / TRT (Delta)"), style = "heading 1") %>%
        body_add_par(paste0(nsim,"x600 samples, ",
                            "on ", ifelse(c=="3"|c=="3.1",3,20),"y censored data."), style = "Normal") %>%
        body_add_table(rmst_tbl_bias[,c(1:2,3:11)], style = "table_template") %>%
        body_add(paste0(ifelse(c=="3"|c=="3.1",3,20),"y Average RMST: PBO / TRT (Delta)"), style = "heading 1") %>%
        body_add_par(paste0(nsim,"x1000 samples, ",
                            "on ", ifelse(c=="3"|c=="3.1",3,20),"y censored data."), style = "Normal") %>%
        body_add_table(rmst_tbl_bias[,c(1:2,12:20)], style = "table_template") %>%
        print(target =paste0(aim,"/RMST Tables/RMST Bias.docx"))
      #output mean delta rmst bias (csv)
      write.csv(rmst_tbl_deltabias,paste0(aim,"/RMST Tables/Delta RMST Bias.csv"))
      #output mean % delta rmst bias (csv)
      write.csv(rmst_tbl_deltaPERCbias,paste0(aim,"/RMST Tables/Delta RMST PercBias.csv"))
      
    }
  }
}

end_time <- Sys.time()
t_rmsttable <- round(end_time - start_time,3)

# TOTAL TIMES

print(t_panelplots)
print(t_rmstplots)
print(t_rmsttable)

