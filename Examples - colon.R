#####################################################
#  consSurv Paper Example                           #
#  Angus Jennings					20Sept23                  #
#####################################################

###############################################
# SETUP       				                        #
###############################################

#### WD/LIBS/SOURCE ####

#update to working directory; where to export outputs
wd <- "WORKING DIRECTORY; OUTPUT"
#ConstrainNR.R
cnr <- "WORKING DIRECTORY; ConstrainNR.R"

#where plots/datasets will output
setwd(wd)

#install required packages
list.of.packages <- c("survival", "pracma", "survPen", "tidyverse", "splines2", "ggthemes", "patchwork", "rlist", "vcd", "ggnewscale")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))suppressPackageStartupMessages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

#import functions from ConstrainNR.R
source(cnr)

#GGTHEME

tt <- theme_base() +
  theme(plot.background =  element_blank(),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey95"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'dashed', colour = "grey95"))

#Reverse Spline Fun

naturalSpline2 <- function (x, log=F, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = NULL, derivs = 0L, integral = FALSE, ...){
  
  if(log){
    x <- log(x)
    if(!(is.null(knots)))knots <- log(knots)
    if(!(is.null(Boundary.knots)))Boundary.knots <- log(Boundary.knots)
  }
  
  x <- -x
  if(!(is.null(knots)))knots <- -knots
  if(!(is.null(Boundary.knots)))Boundary.knots <- -Boundary.knots
  
  naturalSpline(x, df, knots, intercept, Boundary.knots, derivs, integral, ...)
}


###############################################
# METHODS SPLINE PLOT - Fig 1                 #
###############################################

x <- 1:20
y <- c(0,1,1,2,3,5,7,9,9,9,8,5,3,3,2,2,3,3,2,1)/5

knots <- seq(1,20, length.out=8)

k <- knots[-c(1,length(knots))]
bk <- knots[c(1,length(knots))]

t <- seq(0,25,length.out=100)
nS <- naturalSpline2(t, knots=k, Boundary.knots =bk)
nS <- cbind(Int=1,nS)

#PLOT1
nS_gg1 <- as.data.frame(nS)
nS_gg1$t <- t
nS_gg1 <- pivot_longer(as.data.frame(nS_gg1),-t)
nS_gg1$name <- factor(nS_gg1$name, levels = c("Int",1:(length(knots)-1)))

#PLOT2
l1 <- lm(y~naturalSpline2(x, knots=k, Boundary.knots = bk))

l1_c <- coef(l1)
names(l1_c) <- c("Int",1:(length(knots)-1))

nS_gg2 <- nS_gg1
nS_gg2$l1_c <- l1_c[nS_gg2$name]
nS_gg2$value <- nS_gg2$l1_c*nS_gg2$value

l1_tot <- summarise(group_by(nS_gg2,t),value=sum(value))
l1_tot$name <- "Sum"
l1_tot$l1_c <- NA

nS_gg2 <- rbind(nS_gg2, l1_tot)

#PLOT3
l2 <- lm(y~naturalSpline2(x, knots=k, Boundary.knots = bk)[,-1] - 1)

l2_c <- c(0,0,coef(l2))
names(l2_c) <- c("Int",1:(length(knots)-1))

nS_gg3 <- nS_gg1
nS_gg3$l2_c <- l2_c[nS_gg3$name]
nS_gg3$value <- nS_gg3$l2_c*nS_gg3$value

l2_tot <- summarise(group_by(nS_gg3,t),value=sum(value))
l2_tot$name <- "Sum"
l2_tot$l2_c <- NA

nS_gg3 <- rbind(nS_gg3, l2_tot)

#FULL

a.max <- (nS_gg1 %>% group_by(name) %>% summarise(max=max(value)))$max
b.max <- (nS_gg2 %>% filter(!(name=="Sum")) %>% group_by(name) %>% summarise(max=max(ifelse(l1_c>=0,1,-1)*value)))$max
c.max <- (nS_gg3 %>% filter(!(name=="Sum")) %>% group_by(name) %>% summarise(max=max(ifelse(l2_c>=0,1,-1)*value)))$max

weights <- data.frame(name=c("Int",1:(length(knots)-1)),a=1,b=coef(l1),c=c(0,0,coef(l2)),x=rev(knots))
weights$y.a=a.max+0.1
weights$y.b=ifelse(weights$b>=0,1,-1)*(b.max+0.1)
weights$y.c=ifelse(weights$c>=0,1,-1)*(c.max+0.1)

weights$name <- factor(weights$name, levels=c("Int",1:(length(knots)-1)))

cols <- c(RColorBrewer::brewer.pal(length(k)+2,"Set1"),1)
a=0.6

gg1 <- ggplot() + geom_vline(xintercept = knots, color="grey90") + geom_hline(yintercept=0) + 
  geom_line(data=nS_gg1, mapping=aes(x=t,y=value,color=name,linetype=name), size=0.8, alpha=a) + geom_point(aes(x=x,y=y)) + 
  ylim(-0.7,2) + guides(color="none",linetype="none") + scale_color_manual(values=cols) + labs(y="y",title="A. Raw Data and Unweighted Spline Variables") +
  theme_bw()+ theme(axis.title.x = element_blank()) + scale_x_continuous(limits=c(0,25),expand=c(0,0)) + 
  scale_linetype_manual(values=c(rep(1,length(knots)),2)) +
  geom_text(data=weights,mapping=aes(x=x,y=y.a,color=name,label=round(a,1)))

gg2 <- ggplot() + geom_vline(xintercept = knots, color="grey90") + geom_hline(yintercept=0) + 
  geom_line(data=nS_gg2, mapping=aes(x=t,y=value,color=name,linetype=name), size=0.8, alpha=a) + geom_point(aes(x=x,y=y)) + 
  ylim(-0.7,2) + guides(color="none",linetype="none") + scale_color_manual(values=cols) + labs(y="y",title="B. Scaled Spline Variables and their Sum") +
  theme_bw()+ theme(axis.title.x = element_blank()) + scale_x_continuous(limits=c(0,25),expand=c(0,0)) + 
  scale_linetype_manual(values=c(rep(1,length(knots)),2)) +
  geom_text(data=weights,mapping=aes(x=x,y=y.b,color=name,label=round(b,1)))

gg3 <- ggplot() + geom_vline(xintercept = knots, color="grey90") + geom_hline(yintercept=0) + 
  geom_line(data=nS_gg3, mapping=aes(x=t,y=value,color=name,linetype=name), size=0.8, alpha=a) + geom_point(aes(x=x,y=y)) + 
  ylim(-0.7,2) + scale_color_manual(values=cols) + labs(y="y",x="x",color="Spline Var",linetype="Spline Var",title="C. As B but with Spline Int/Variable 1 Constrained to 0") +
  theme_bw() + scale_x_continuous(limits=c(0,25),expand=c(0,0)) + 
  scale_linetype_manual(values=c(rep(1,length(knots)),2)) +
  new_scale_colour()+
  geom_text(data=weights,mapping=aes(x=x,y=y.c,color=name,label=round(c,1))) + 
  scale_color_manual(values=cols) + guides(color="none")

gg <- gg1/gg2/gg3 + plot_layout(guides = "collect") & theme(panel.grid = element_blank(), plot.caption = element_text(hjust=0))

ggsave("Figure 1 - Spline Plot.png", gg, width=22, height=30, units = "cm")


###############################################
# EXAMPLE 1 - Colon Data                      #
###############################################

#format data - colon from package survival
dat <- colon[colon$etype==2 & colon$rx %in% c("Obs","Lev+5FU"),]

dat$yrs <- dat$time/365.25
dat$rx <- as.numeric(dat$rx=="Lev+5FU")

times <- seq(0.01, 10, length.out=100)

#knots
deg <- 5

k <- quantile(dat$yrs[dat$status==1], seq(0,1,length.out=deg))
k95 <- quantile(dat$yrs[dat$status==1], sort(c(seq(0,1,length.out=deg-1),0.95)))
k95[deg] <- 10
# log time models predict 0 pbo haz for low times, making HRs weird, keep to id time scale?
# just predict from 0.3 yrs (~first pbo event)?

maxk <- max(k)

#plot limits
hazmin <- hrmin <- 0
hazmax <- 0.3
hrmax <- 2


## UNCONSTRAINED MODEL ##

#model
sP <- survPen_cons(formula = ~ rx*naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)
#predictions
sP_pred_O <- list.cbind(predict(sP, newdata = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sP_pred_O) <- paste0("O.",colnames(sP_pred_O))
sP_pred_L5FU <- list.cbind(predict(sP, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sP_pred_L5FU) <- paste0("L5FU.",colnames(sP_pred_L5FU))
sP_pred_HR <- list.cbind(predict(sP, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), newdata.ref = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), type = "HR"))

sP_preds <- as.data.frame(cbind(times, sP_pred_O, sP_pred_L5FU, sP_pred_HR))

#plots
g1_1_h <- ggplot(sP_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = O.haz.inf, ymax = O.haz.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = L5FU.haz.inf, ymax = L5FU.haz.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=O.haz), color=2, size=0.9) +
  geom_line(aes(y=L5FU.haz), color=3, size=0.9) +
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) +
  labs(title = "A: Unconstrained Estimates", x="Time (years)", y="Hazard") + 
  tt + theme(axis.title.x = element_blank(), plot.subtitle=element_text(size=10)) + coord_cartesian(ylim=c(hazmin, hazmax))

g1_1_hr <- ggplot(sP_preds, aes(x=times)) + 
  geom_hline(yintercept=1) +
  geom_ribbon(aes(ymin = HR.inf, ymax = HR.sup), alpha=0.2) +
  geom_line(aes(y=HR), size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1,], aes(x=yrs), alpha=0.5) + 
  labs(x="Time (years)", y="HR") + tt + coord_cartesian(ylim=c(hrmin, hrmax))

g1_1_s <- ggplot(sP_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = O.surv.inf, ymax = O.surv.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = L5FU.surv.inf, ymax = L5FU.surv.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=O.surv), color=2, size=0.9) +
  geom_line(aes(y=L5FU.surv), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) +
  lims(y=c(0,1)) + 
  labs(x="Time (years)", y="Survival Probability") + tt + theme(axis.title.x = element_blank())

g1_1  <- g1_1_h / g1_1_s / g1_1_hr


## CONSTRAINED MODEL : HR=1 FROM UPPER OBVS ##

b.ini <- rep(0,length(sP$coefficients))

uncons.keep <- c(1, 3:(deg+1))
cons.zero <- c(2, deg+10)

b.ini[uncons.keep] <- sP$coefficients[uncons.keep]
b.ini[cons.zero] <- 0.000000000001

sPcons <- survPen_cons(formula = ~ rx + naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx:naturalSpline2(yrs, Boundary.knots = k95[c(1,deg)], knots=k95[-c(1,deg)])+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1, beta.ini = b.ini,
                   cons=c(uncons.keep,cons.zero))

#boot cis

boot.fun1 <- function(data, indices){
  
  dat <- data[indices,]
  
  k <- quantile(dat$yrs[dat$status==1], seq(0,1,length.out=deg))
  k95 <- quantile(dat$yrs[dat$status==1], sort(c(seq(0,1,length.out=deg-1),0.95)))
  k95[deg] <- 10
  
  sP <- survPen_cons(formula = ~ rx*naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)])+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)
  
  b.ini <- rep(0,length(sP$coefficients))
  
  uncons.keep <- c(1, 3:(deg+1))
  cons.zero <- c(2, deg+10)
  
  b.ini[uncons.keep] <- sP$coefficients[uncons.keep]
  b.ini[cons.zero] <- 0.000000000001
  
  sPcons <- survPen_cons(formula = ~ rx + naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx:naturalSpline2(yrs, Boundary.knots = k95[c(1,deg)], knots=k95[-c(1,deg)])+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1, beta.ini = b.ini,
                    cons=c(uncons.keep,cons.zero))
  
  p0 <- predict(sPcons, newdata = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0))
  p1 <- predict(sPcons, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0))
  p10 <- predict(sPcons, newdata = data.frame(yrs=seq(0.01, 9.9, length.out=100), rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0),
                 newdata.ref = data.frame(yrs=seq(0.01, 9.9, length.out=100), rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), type = "HR")
  
  return(log(c(p0$haz,p1$haz,p0$surv,p1$surv,p10$HR)))
  
}

set.seed(719831)
bt1 <- boot::boot(dat, boot.fun1, R=1000)

cis1 <- t(
  sapply(1:500, FUN=function(i){
    boot::boot.ci(bt1, conf=0.95, type = c("perc"), index=i)[-(1:3)][[1]]
  })
)

cis1 <- exp( cis1[,4:5] )
colnames(cis1) <- c("lci","uci")

sPcons_pred_O <- list.cbind(predict(sPcons, newdata = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
sPcons_pred_O[,"haz.inf"] <- cis1[,"lci"][1:100]
sPcons_pred_O[,"haz.sup"] <- cis1[,"uci"][1:100]
sPcons_pred_O[,"surv.inf"] <- cis1[,"lci"][201:300]
sPcons_pred_O[,"surv.sup"] <- cis1[,"uci"][201:300]
colnames(sPcons_pred_O) <- paste0("O.",colnames(sPcons_pred_O))
sPcons_pred_L5FU <- list.cbind(predict(sPcons, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
sPcons_pred_L5FU[,"haz.inf"] <- cis1[,"lci"][101:200]
sPcons_pred_L5FU[,"haz.sup"] <- cis1[,"uci"][101:200]
sPcons_pred_L5FU[,"surv.inf"] <- cis1[,"lci"][301:400]
sPcons_pred_L5FU[,"surv.sup"] <- cis1[,"uci"][301:400]
colnames(sPcons_pred_L5FU) <- paste0("L5FU.",colnames(sPcons_pred_L5FU))
sPcons_pred_HR <- list.cbind(predict(sPcons, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), newdata.ref = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), type = "HR"))
sPcons_pred_HR[,"HR.inf"] <- cis1[,"lci"][401:500]
sPcons_pred_HR[,"HR.sup"] <- cis1[,"uci"][401:500]

sPcons_preds <- as.data.frame(cbind(times, timesHR=seq(0.01, 9.9, length.out=100), sPcons_pred_O, sPcons_pred_L5FU, sPcons_pred_HR))

#plot

g1_2_h <- ggplot(sPcons_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = O.haz.inf, ymax = O.haz.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = L5FU.haz.inf, ymax = L5FU.haz.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=O.haz), color=2, size=0.9) +
  geom_line(aes(y=L5FU.haz), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) + 
  labs(title="B: Constrained Treatment Effect Estimates", x="Time (years)", y="Hazard") + 
  tt + theme(axis.title.x = element_blank(), plot.subtitle=element_text(size=10)) + coord_cartesian(ylim=c(hazmin, hazmax))

g1_2_hr <- ggplot(sPcons_preds) + 
  geom_hline(yintercept=1) +
  geom_ribbon(aes(x=timesHR, ymin = HR.inf, ymax = HR.sup), alpha=0.2) +
  geom_line(aes(x=times, y=HR), size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1,], aes(x=yrs), alpha=0.5) + 
  labs(x="Time (years)", y="HR") + tt + coord_cartesian(ylim=c(hrmin, hrmax))

g1_2_s <- ggplot(sPcons_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = O.surv.inf, ymax = O.surv.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = L5FU.surv.inf, ymax = L5FU.surv.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=O.surv), color=2, size=0.9) +
  geom_line(aes(y=L5FU.surv), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) +
  lims(y=c(0,1)) + 
  labs(x="Time (years)", y="Survival Probability") + tt + theme(axis.title.x = element_blank())

g1_2  <- g1_2_h / g1_2_s / g1_2_hr

g1 <- (g1_1 | g1_2)  - grid_legend("center", labels=c("Placebo Hazard/\nSurvival Probability", "Lev+5FU Hazard/\nSurvival Probability"), lty=1, col=c(2,3), frame=F, lwd=2, gp=gpar(cex=1.2), hgap = unit(3, "lines"), vgap = unit(1, "lines")) + plot_layout(widths = c(10,2))

ggsave("Figure 2 - Cons vs Uncons Model Ests.png", g1, height = 10, width=16, units="in")


## CONSTRAINED MODEL : PH THEN WANING ##

sPph <- survPen_cons(formula = ~ naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)
sPphwan <- survPen_cons(formula = ~ naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx:naturalSpline2(yrs, Boundary.knots = c(maxk,10), knots=seq(maxk,10,length.out=4)[-c(1,4)])[,3]+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)

sPph_pred_O <- list.cbind(predict(sPph, newdata = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sPph_pred_O) <- paste0("ph.O.",colnames(sPph_pred_O))
sPph_pred_L5FU <- list.cbind(predict(sPph, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sPph_pred_L5FU) <- paste0("ph.L5FU.",colnames(sPph_pred_L5FU))
sPph_pred_HR <- list.cbind(predict(sPph, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), newdata.ref = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), type = "HR"))
colnames(sPph_pred_HR) <- paste0("ph.",colnames(sPph_pred_HR))

sPphwan_pred_O <- list.cbind(predict(sPphwan, newdata = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sPphwan_pred_O) <- paste0("wan.O.",colnames(sPphwan_pred_O))
sPphwan_pred_L5FU <- list.cbind(predict(sPphwan, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0)))
colnames(sPphwan_pred_L5FU) <- paste0("wan.L5FU.",colnames(sPphwan_pred_L5FU))
sPphwan_pred_HR <- list.cbind(predict(sPphwan, newdata = data.frame(yrs=times, rx=1, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), newdata.ref = data.frame(yrs=times, rx=0, sex=0, age=59.58, obstruct=0, perfor=0, adhere=0, node4=0, extent=3, surg=0), type = "HR"))
colnames(sPphwan_pred_HR) <- paste0("wan.",colnames(sPphwan_pred_HR))

sPph_preds <- as.data.frame(cbind(times, sPph_pred_O, sPphwan_pred_O, sPph_pred_L5FU, sPphwan_pred_L5FU, sPph_pred_HR, sPphwan_pred_HR))

#ph plots
g2_1_h <- ggplot(sPph_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = ph.O.haz.inf, ymax = ph.O.haz.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = ph.L5FU.haz.inf, ymax = ph.L5FU.haz.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=ph.L5FU.haz), color=3, size=0.9) +
  geom_line(aes(y=ph.O.haz), color=2, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) + 
  labs(title="A: Unconstrained Proportional Hazards Estimates", x="Time (years)", y="Hazard") + 
  tt + theme(axis.title.x = element_blank(), plot.subtitle=element_text(size=10)) + coord_cartesian(ylim=c(hazmin, hazmax))

g2_1_hr <- ggplot(sPph_preds, aes(x=times)) + 
  geom_hline(yintercept=1) +
  geom_ribbon(aes(ymin = ph.HR.inf, ymax = ph.HR.sup), alpha=0.2) +
  geom_line(aes(y=ph.HR), size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1,], aes(x=yrs), alpha=0.5) + 
  labs(x="Time (years)", y="HR") + tt + coord_cartesian(ylim=c(hrmin, hrmax))

g2_1_s <- ggplot(sPph_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = ph.O.surv.inf, ymax = ph.O.surv.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = ph.L5FU.surv.inf, ymax = ph.L5FU.surv.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=ph.O.surv), color=2, size=0.9) +
  geom_line(aes(y=ph.L5FU.surv), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) +
  lims(y=c(0,1)) + 
  labs(x="Time (years)", y="Survival Probability") + tt+ theme(axis.title.x = element_blank())

g2_1  <- g2_1_h / g2_1_s / g2_1_hr

#waning plots
g2_2_h <- ggplot(sPph_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = wan.O.haz.inf, ymax = wan.O.haz.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = wan.L5FU.haz.inf, ymax = wan.L5FU.haz.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=wan.O.haz), color=2, size=0.9) +
  geom_line(aes(y=wan.L5FU.haz), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) + 
  labs(title="B: Proportional Hazards + Waning Estimates", x="Time (years)", y="Hazard") + 
  tt + theme(axis.title.x = element_blank(), plot.subtitle=element_text(size=10)) + coord_cartesian(ylim=c(hazmin, hazmax))

g2_2_hr <- ggplot(sPph_preds, aes(x=times)) + 
  geom_hline(yintercept=1) +
  geom_ribbon(aes(ymin = wan.HR.inf, ymax = wan.HR.sup), alpha=0.2) +
  geom_line(aes(y=wan.HR), size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1,], aes(x=yrs), alpha=0.5) + 
  labs(x="Time (years)", y="HR") + tt + coord_cartesian(ylim=c(hrmin, hrmax))

g2_2_s <- ggplot(sPph_preds, aes(x=times)) + 
  geom_ribbon(aes(ymin = wan.O.surv.inf, ymax = wan.O.surv.sup), fill = 2, alpha=0.2) +
  geom_ribbon(aes(ymin = wan.L5FU.surv.inf, ymax = wan.L5FU.surv.sup), fill = 3, alpha=0.2) +
  geom_line(aes(y=wan.O.surv), color=2, size=0.9) +
  geom_line(aes(y=wan.L5FU.surv), color=3, size=0.9) +
  geom_vline(xintercept = maxk, linetype="dashed", color="grey")+
  geom_rug(data=dat[dat$status==1 & dat$rx==0,], aes(x=yrs), alpha=0.5, color=2) +
  geom_rug(data=dat[dat$status==1 & dat$rx==1,], aes(x=yrs), alpha=0.5, color=3) +
  lims(y=c(0,1)) + 
  labs(x="Time (years)", y="Survival Probability") + tt+ theme(axis.title.x = element_blank())

g2_2  <- g2_2_h / g2_2_s / g2_2_hr

g2 <- (g2_1 | g2_2)  - grid_legend("center", labels=c("Placebo Hazard/\nSurvival Probability", "Lev+5FU Hazard/\nSurvival Probability"), lty=1, col=c(2,3), frame=F, lwd=2, gp=gpar(cex=1.2), hgap = unit(3, "lines"), vgap = unit(1, "lines")) + plot_layout(widths = c(10,2))

ggsave("Appendix 1 - Cons vs Uncons PH Model Ests.png", g2, height = 10, width=16, units="in")


## RMSTS ##

#rmsts
sP_ssP <- standSurvPen(object=sP, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)
sPcons_ssP <- standSurvPen(object=sPcons, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)

sPph_ssP <- standSurvPen(object=sPph, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)
sPphwan_ssP <- standSurvPen(object=sPphwan, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)

#rerun with later constraint
k95_ext <- k95
k95_ext[deg] <- 20

sP_2 <- survPen_cons(formula = ~ rx*naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)])+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)

b.ini <- rep(0,length(sP_2$coefficients))
uncons.keep <- c(1, 3:(deg+1))
cons.zero <- c(2, deg+10)
b.ini[uncons.keep] <- sP_2$coefficients[uncons.keep]
b.ini[cons.zero] <- 0.000000000001
sPcons_2 <- survPen_cons(formula = ~ rx + naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx:naturalSpline2(yrs, Boundary.knots = k95_ext[c(1,deg)], knots=k95_ext[-c(1,deg)])+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1, beta.ini = b.ini,
                  cons=c(uncons.keep,cons.zero))

sPph_2 <- survPen_cons(formula = ~ naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)
sPphwan_2 <- survPen_cons(formula = ~ naturalSpline2(yrs, Boundary.knots = k[c(1,deg)], knots=k[-c(1,deg)]) + rx:naturalSpline2(yrs, Boundary.knots = c(maxk,15), knots=seq(maxk,15,length.out=4)[-c(1,4)])[,3]+ sex + age + obstruct + perfor + adhere + node4 + extent + surg, data=dat, t1=yrs, event=status==1)

#recalc rmsts
sP_ssP_2 <- standSurvPen(object=sP_2, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)
sPcons_ssP_2 <- standSurvPen(object=sPcons_2, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)

sPph_ssP_2 <- standSurvPen(object=sPph_2, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)
sPphwan_ssP_2 <- standSurvPen(object=sPphwan_2, data=dat, at=list(O=list(rx=0),L=list(rx=1)), times=30, type=c("rmst"), n.legendre=30)



