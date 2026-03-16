#library(metafor)
#library("mvmeta")
library(mixmeta)
library(MASS)
library(dplyr)

RUNVARS <- commandArgs(T)
RUNMETHOD <- as.integer(RUNVARS[1])
id <- as.integer(RUNVARS[2]) #dim(grid)
ncores <- as.integer(RUNVARS[3])
FILEPATH <- RUNVARS[4]

################################################################################
cond <- NULL
cond <- expand.grid(sz_vec=c(10,20,50,100,300),
                    delta_vec=c(0,0.2,0.5,0.8), #c(0,0.2,0.5),
                    ttau_vec=c(0,0.006,0.02,0.05),
                    r_vec=c(0,0.1,0.3,0.5),#c(0,0.1,0.3,0.5),
                    est.r_vec=c(-0.5,-0.3,-0.1,0,0.1,0.3,0.5),
                    prop_vec=c(0,0.5,0.8))

condition<-do.call(cbind, cond)
index.l<-c(1:nrow(condition))
################################################################################
it<-c(1:1)

grid <- expand.grid(it,index.l)
grid<-as.data.frame(grid)
colnames(grid)<-c("it","index")

args <- grid[id, ]

do_one <- function(it,index){
##############
data<-read.csv("d5.csv",header=F)

  sz=condition[index,1]
  # population parameters
  delta1<-condition[index,2]
  delta2<-condition[index,2]
  ttau1<-condition[index,3]
  ttau2<-condition[index,3]
  rlevel1<-condition[index,4] # level1 y correlation
  rlevel2<-condition[index,4] # level2 true delta correlation
  sigmalevel2<-matrix(c(ttau1,rlevel2*sqrt(ttau1*ttau2),rlevel2*sqrt(ttau1*ttau2),ttau2),2,2)
  e.rlevel1<-condition[index,5]
  prop<-condition[index,6]
  
  set.seed(1)
  nn<-sample(c(1:length(data[,1])),size=sz,replace=T)
  n1<-data[nn,1]
  n2<-data[nn,2]
  
mu1=mu2=ci1=ci2=e.tau1=e.tau2=rho=Qp=mu1p=mu2p=mu1.se=mu2.se=c()
mu1.uni=mu2.uni=ci1.uni=ci2.uni=e.tau1.uni=e.tau2.uni=Qp1.uni=Qp2.uni=mu1p.uni=mu2p.uni=mu1.se.uni=mu2.se.uni=c()

for (i in 1:1000){
  
##### Simulate the data
  set.seed(i)
ttd<-mvrnorm(sz,c(delta1,delta2), sigmalevel2)

d1=d2=c()
for (id in 1:sz){ #for idth study
  sigmalevel1<-matrix(c(1,rlevel1,rlevel1,1),2,2)
  #group1:control
  y.c<-mvrnorm(n1[id],rep(0,2), sigmalevel1)
  y1.c<-y.c[,1]
  y2.c<-y.c[,2]
  #group2
  y.t<-mvrnorm(n2[id],ttd[id,], sigmalevel1)
  y1.t<-y.t[,1]
  y2.t<-y.t[,2]
  
#outcome 1
  y1.tm<-mean(y1.t)
  y1.cm<-mean(y1.c)
  s1.t<-var(y1.t)
  s1.c<-var(y1.c)
  sp<-sqrt( ((n1[id]-1)*s1.t+(n2[id]-1)*s1.c)/(n1[id]+n2[id]-2))
  d1[id]<- (y1.tm- y1.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))

#outcome 1  
  y2.tm<-mean(y2.t)
  y2.cm<-mean(y2.c)
  s2.t<-var(y2.t)
  s2.c<-var(y2.c)
  sp<-sqrt( ((n1[id]-1)*s2.t+(n2[id]-1)*s2.c)/(n1[id]+n2[id]-2))
  d2[id]<- (y2.tm- y2.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))
}

v1i<-(n1+n2)/n1/n2+d1^2/(2*(n1+n2))
v2i<-(n1+n2)/n1/n2+d2^2/(2*(n1+n2))

na_indices <- sample(1:sz, size = floor(sz * prop))
d2[na_indices ]<-NA

### Analyze the data
####multivariate 
cov12i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d1*d2/(2*(n1+n2))
model<-try(mixmeta(cbind(d1,d2), cbind(v1i,cov12i,v2i), method="reml"))

if( class(model)[1] !="try-error"){
res.table<-summary(model)$ coefficients

mu1[i]<-res.table[1,1]
mu2[i]<-res.table[2,1]
mu1.se[i]<-res.table[1,2]
mu2.se[i]<-res.table[2,2]
ci1[i]<-between(delta1,res.table[1,5],res.table[1,6])
ci2[i]<-between(delta2,res.table[2,5],res.table[2,6])
mu1p[i]<-res.table[1,4]
mu2p[i]<-res.table[2,4]
e.tau1[i]<-model$ Psi[1,1] 
e.tau2[i]<-model$ Psi[2,2] 
rho[i]<-summary(model)$corRandom[1,2]
Qp[i]<-summary(model)$ qstat $pvalue[1]} else {
  mu1[i]=mu2[i]=ci1[i]=ci2[i]=e.tau1[i]=e.tau2[i]=rho[i]=Qp[i]=mu1.se[i]=mu2.se[i]=NA
}

#####univarate
model1<-try(mixmeta(d1, v1i, method="reml"))#try(rma(d1, v1i, method="REML"))
model2<-try(mixmeta(d2, v2i, method="reml"))#try(rma(d2, v2i, method="REML"))


if( class(model1)[1] !="try-error"){
  model1<- summary(model1)
mu1.uni[i]<-model1$coefficients[1]
mu1.se.uni[i]<-model1$coefficients[2]
ci1.uni[i]<-between(delta1,model1$coefficients[5],model1$coefficients[6])
mu1p.uni[i]<-model1$coefficients[4]
e.tau1.uni[i]<-model1$ Psi  
Qp1.uni[i]<-model1$ qstat$ pvalue
} else {
  mu1.uni[i]=ci1.uni[i]=mu1p.uni[i]=e.tau1.uni[i]=Qp1.uni[i]=mu1.se.uni[i]=NA
}

if( class(model2)[1] !="try-error"){
  model2<- summary(model2)
  mu2.uni[i]<-model2$coefficients[1]
  mu2.se.uni[i]<-model2$coefficients[2]
  ci2.uni[i]<-between(delta2,model2$coefficients[5],model2$coefficients[6])
  mu2p.uni[i]<-model2$coefficients[4]
  e.tau2.uni[i]<-model2$ Psi  
  Qp2.uni[i]<-model2$ qstat$ pvalue
  # mu2.uni[i]<-model2$b   
  # ci2.uni[i]<-between(delta2,model2$ci.lb,model2$ci.ub)
  # e.tau2.uni[i]<-model2$tau2
  # Qp2.uni[i]<-model2$QEp 
} else {
  mu2.uni[i]=ci2.uni[i]=mu2p.uni[i]=e.tau2.uni[i]=Qp2.uni[i]=mu2.se.uni[i]=NA
}

}

convergence.mul<-length(na.omit(mu1))/1000
convergence.uni1<-length(na.omit(mu1.uni))/1000
convergence.uni2<-length(na.omit(mu2.uni))/1000

### summarize results
if (delta1!=0){
bias.mu=c(mean((na.omit(mu1)-delta1)/delta1),mean((na.omit(mu2)-delta2)/delta2),
          mean((na.omit(mu1.uni)-delta1)/delta1),mean((na.omit(mu2.uni)-delta2)/delta2))} else {
bias.mu=c(mean((na.omit(mu1)-delta1)),mean((na.omit(mu2)-delta2)),
          mean((na.omit(mu1.uni)-delta1)),mean((na.omit(mu2.uni)-delta2)))}

t.mu1.se<-sd(na.omit(mu1))
t.mu2.se<-sd(na.omit(mu2))
bias.muse=c(mean((na.omit(mu1.se)-t.mu1.se)/t.mu1.se),mean((na.omit(mu2.se)-t.mu2.se)/t.mu2.se),
            mean((na.omit(mu1.se.uni)-t.mu1.se)/t.mu1.se),mean((na.omit(mu2.se.uni)-t.mu2.se)/t.mu2.se))

if (ttau1!=0){
bias.tau=c(mean((na.omit(e.tau1)-ttau1)/ttau1),mean((na.omit(e.tau2)-ttau2)/ttau2),
             mean((na.omit(e.tau1.uni)-ttau1)/ttau1),mean((na.omit(e.tau2.uni)-ttau2)/ttau2))} else {
bias.tau=c(mean((na.omit(e.tau1)-ttau1)),mean((na.omit(e.tau2)-ttau2)),
             mean((na.omit(e.tau1.uni)-ttau1)),mean((na.omit(e.tau2.uni)-ttau2)))}

if (rlevel2!=0){
bias.rho<- mean((na.omit(rho)-rlevel2)/rlevel2)} else {
bias.rho<- mean((na.omit(rho)-rlevel2))
}

 coveragemu1<-mean(na.omit(ci1))
 coveragemu2<-mean(na.omit(ci2))
 coveragemu1.uni<-mean(na.omit(ci1.uni))
 coveragemu2.uni<-mean(na.omit(ci2.uni))
 
 rmsemu<-sqrt(c(mean((na.omit(mu1)-delta1)^2),mean((na.omit(mu2)-delta2)^2),
          mean((na.omit(mu1.uni)-delta1)^2),mean((na.omit(mu2.uni)-delta2)^2)))
 
 rmsetau<-sqrt(c(mean((na.omit(e.tau1)-ttau1)^2),mean((na.omit(e.tau2)-ttau2)^2),
            mean((na.omit(e.tau1.uni)-ttau1)^2),mean((na.omit(e.tau2.uni)-ttau2)^2)))
 
 powermu1<-mean(na.omit(mu1p)<=0.05)
 powermu2<-mean(na.omit(mu2p)<=0.05)
 powermu1.uni<-mean(na.omit(mu1p.uni)<=0.05)
 powermu2.uni<-mean(na.omit(mu2p.uni)<=0.05)
 powerQp<-mean(na.omit(Qp)<=0.05)
 powerQp1.uni<-mean(na.omit(Qp1.uni)<=0.05)
 powerQp2.uni<-mean(na.omit(Qp2.uni)<=0.05)
 
 result.table<-c(condition[index,],convergence.mul,convergence.uni1,convergence.uni2,
                 bias.mu,bias.muse,bias.tau,bias.rho, coveragemu1, coveragemu2,
                 coveragemu1.uni, coveragemu2.uni,rmsemu,rmsetau,
                 powermu1, powermu2, powermu1.uni, powermu2.uni,
                 powerQp, powerQp1.uni, powerQp2.uni)
 
 fn <- paste0("result/",index,"result.csv")
 write.table(t(result.table), file = fn,  sep = ",", col.names = FALSE,append = TRUE)
 
}

#####################################################
do.call(do_one, args)
 
 