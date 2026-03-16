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
                    delta_vec=c(0.2), #c(0,0.2,0.5),
                    ttau_vec=c(0,0.006,0.02,0.05),
                    r_vec=c(-0.3,0,0.5),#c(0,0.1,0.3,0.5),
                    est.r_vec=c(-0.3,0,0.5),
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
  delta3<-condition[index,2]
  delta4<-condition[index,2]
  ttau1<-condition[index,3]
  ttau2<-condition[index,3]
  ttau3<-condition[index,3]
  ttau4<-condition[index,3]
  rlevel1<-condition[index,4] # level1 y correlation
  rlevel2<-condition[index,4] # level2 true delta correlation
  taus  <- c(ttau1, ttau2, ttau3, ttau4)
  R     <- matrix(rlevel2, 4, 4)
  diag(R) <- 1
  sigmalevel2 <- diag(sqrt(taus)) %*% R %*% diag(sqrt(taus))
  e.rlevel1<-condition[index,5]
  prop<-condition[index,6]
  
  set.seed(1)
  nn<-sample(c(1:length(data[,1])),size=sz,replace=T)
  n1<-data[nn,1]
  n2<-data[nn,2]
  
  mu1=mu2=mu3=mu4=ci1=ci2=ci3=ci4=e.tau1=e.tau2=e.tau3=e.tau4=c()
  rho12=rho13=rho14=rho23=rho24=rho34=c()
  Qp=mu1p=mu2p=mu3p=mu4p=mu1.se=mu2.se=mu3.se=mu4.se=c()
  mu1.uni=mu2.uni=mu3.uni=mu4.uni=ci1.uni=ci2.uni=ci3.uni=ci4.uni=c()
  e.tau1.uni=e.tau2.uni=e.tau3.uni=e.tau4.uni=c()
  Qp1.uni=Qp2.uni=Qp3.uni=Qp4.uni=mu1p.uni=mu2p.uni=mu3p.uni=mu4p.uni=c()
  mu1.se.uni=mu2.se.uni=mu3.se.uni=mu4.se.uni=c()
  
  for (i in 1:1000){
    
    ##### Simulate the data
    set.seed(i)
    ttd<-mvrnorm(sz,c(delta1,delta2,delta3,delta4), sigmalevel2)
    
    d1=d2=d3=d4=c()
    for (id in 1:sz){ #for idth study
      sigmalevel1<-matrix(rlevel1, 4, 4)
      diag(sigmalevel1) <- 1
      #group1:control
      y.c<-mvrnorm(n1[id],rep(0,4), sigmalevel1)
      y1.c<-y.c[,1]
      y2.c<-y.c[,2]
      y3.c<-y.c[,3]
      y4.c<-y.c[,4]
      #group2
      y.t<-mvrnorm(n2[id],ttd[id,], sigmalevel1)
      y1.t<-y.t[,1]
      y2.t<-y.t[,2]
      y3.t<-y.t[,3]
      y4.t<-y.t[,4]
      
      #outcome 1
      y1.tm<-mean(y1.t)
      y1.cm<-mean(y1.c)
      s1.t<-var(y1.t)
      s1.c<-var(y1.c)
      sp<-sqrt( ((n1[id]-1)*s1.t+(n2[id]-1)*s1.c)/(n1[id]+n2[id]-2))
      d1[id]<- (y1.tm- y1.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))
      
      #outcome 2
      y2.tm<-mean(y2.t)
      y2.cm<-mean(y2.c)
      s2.t<-var(y2.t)
      s2.c<-var(y2.c)
      sp<-sqrt( ((n1[id]-1)*s2.t+(n2[id]-1)*s2.c)/(n1[id]+n2[id]-2))
      d2[id]<- (y2.tm- y2.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))
      
      #outcome 3
      y3.tm<-mean(y3.t)
      y3.cm<-mean(y3.c)
      s3.t<-var(y3.t)
      s3.c<-var(y3.c)
      sp<-sqrt( ((n1[id]-1)*s3.t+(n2[id]-1)*s3.c)/(n1[id]+n2[id]-2))
      d3[id]<- (y3.tm- y3.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))
      
      #outcome 4
      y4.tm<-mean(y4.t)
      y4.cm<-mean(y4.c)
      s4.t<-var(y4.t)
      s4.c<-var(y4.c)
      sp<-sqrt( ((n1[id]-1)*s4.t+(n2[id]-1)*s4.c)/(n1[id]+n2[id]-2))
      d4[id]<- (y4.tm- y4.cm)/sp*(1-3/(4*(n1[id]+n2[id])-9))
    }
    
    v1i<-(n1+n2)/n1/n2+d1^2/(2*(n1+n2))
    v2i<-(n1+n2)/n1/n2+d2^2/(2*(n1+n2))
    v3i<-(n1+n2)/n1/n2+d3^2/(2*(n1+n2))
    v4i<-(n1+n2)/n1/n2+d4^2/(2*(n1+n2))
    
    na_indices1 <- sample(1:sz, size = floor(sz * prop))
    d2[na_indices1 ]<-NA
    
    na_indices2 <- sample(1:sz, size = floor(sz * prop))
    d3[na_indices2 ]<-NA
    
    na_indices3 <- sample(1:sz, size = floor(sz * prop))
    d4[na_indices3 ]<-NA
    
    #mat <- cbind(d1, d2, d3, d4)
    #non_na_count <- rowSums(!is.na(mat))
    #table(non_na_count)
    
    ### Analyze the data
    ####multivariate 
    cov12i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d1*d2/(2*(n1+n2))
    cov13i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d1*d3/(2*(n1+n2))
    cov14i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d1*d4/(2*(n1+n2))
    cov23i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d2*d3/(2*(n1+n2))
    cov24i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d2*d4/(2*(n1+n2))
    cov34i<-e.rlevel1*(n1+n2)/n1/n2+e.rlevel1^2*d3*d4/(2*(n1+n2))
    
    model<-try(mixmeta(cbind(d1,d2,d3,d4), 
                       cbind(v1i,cov12i,cov13i,cov14i,
                             v2i,cov23i,cov24i,
                             v3i,cov34i,v4i), method="reml"))
    
    
    if( class(model)[1] !="try-error" & class(try(summary(model), silent=TRUE))[1] != "try-error"){
      res.table<-summary(model)$ coefficients
      
      mu1[i]<-res.table[1,1]
      mu2[i]<-res.table[2,1]
      mu3[i]<-res.table[3,1]
      mu4[i]<-res.table[4,1]
      mu1.se[i]<-res.table[1,2]
      mu2.se[i]<-res.table[2,2]
      mu3.se[i]<-res.table[3,2]
      mu4.se[i]<-res.table[4,2]
      ci1[i]<-between(delta1,res.table[1,5],res.table[1,6])
      ci2[i]<-between(delta2,res.table[2,5],res.table[2,6])
      ci3[i]<-between(delta3,res.table[3,5],res.table[3,6])
      ci4[i]<-between(delta4,res.table[4,5],res.table[4,6])
      mu1p[i]<-res.table[1,4]
      mu2p[i]<-res.table[2,4]
      mu3p[i]<-res.table[3,4]
      mu4p[i]<-res.table[4,4]
      e.tau1[i]<-model$ Psi[1,1] 
      e.tau2[i]<-model$ Psi[2,2] 
      e.tau3[i]<-model$ Psi[3,3] 
      e.tau4[i]<-model$ Psi[4,4] 
      rho12[i]<-summary(model)$corRandom[1,2]
      rho13[i]<-summary(model)$corRandom[1,3]
      rho14[i]<-summary(model)$corRandom[1,4]
      rho23[i]<-summary(model)$corRandom[2,3]
      rho24[i]<-summary(model)$corRandom[2,4]
      rho34[i]<-summary(model)$corRandom[3,4]
      Qp[i]<-summary(model)$ qstat $pvalue[1]} else {
        mu1[i]=mu2[i]=mu3[i]=mu4[i]=ci1[i]=ci2[i]=ci3[i]=ci4[i]=NA
        e.tau1[i]=e.tau2[i]=e.tau3[i]=e.tau4[i]=NA
        rho12[i]=rho13[i]=rho14[i]=rho23[i]=rho24[i]=rho34[i]=NA
        Qp[i]=mu1.se[i]=mu2.se[i]=mu3.se[i]=mu4.se[i]=NA
      }
    
    #####univarate
    model1<-try(mixmeta(d1, v1i, method="reml"))#try(rma(d1, v1i, method="REML"))
    model2<-try(mixmeta(d2, v2i, method="reml"))#try(rma(d2, v2i, method="REML"))
    model3<-try(mixmeta(d3, v3i, method="reml"))
    model4<-try(mixmeta(d4, v4i, method="reml"))
    
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
    } else {
      mu2.uni[i]=ci2.uni[i]=mu2p.uni[i]=e.tau2.uni[i]=Qp2.uni[i]=mu2.se.uni[i]=NA
    }
    
    if( class(model3)[1] !="try-error"){
      model3<- summary(model3)
      mu3.uni[i]<-model3$coefficients[1]
      mu3.se.uni[i]<-model3$coefficients[2]
      ci3.uni[i]<-between(delta3,model3$coefficients[5],model3$coefficients[6])
      mu3p.uni[i]<-model3$coefficients[4]
      e.tau3.uni[i]<-model3$ Psi  
      Qp3.uni[i]<-model3$ qstat$ pvalue
    } else {
      mu3.uni[i]=ci3.uni[i]=mu3p.uni[i]=e.tau3.uni[i]=Qp3.uni[i]=mu3.se.uni[i]=NA
    }
    
    if( class(model4)[1] !="try-error"){
      model4<- summary(model4)
      mu4.uni[i]<-model4$coefficients[1]
      mu4.se.uni[i]<-model4$coefficients[2]
      ci4.uni[i]<-between(delta4,model4$coefficients[5],model4$coefficients[6])
      mu4p.uni[i]<-model4$coefficients[4]
      e.tau4.uni[i]<-model4$ Psi  
      Qp4.uni[i]<-model4$ qstat$ pvalue
    } else {
      mu4.uni[i]=ci4.uni[i]=mu4p.uni[i]=e.tau4.uni[i]=Qp4.uni[i]=mu4.se.uni[i]=NA
    }
    
  }
  
  convergence.mul<-length(na.omit(mu1))/1000
  convergence.uni1<-length(na.omit(mu1.uni))/1000
  convergence.uni2<-length(na.omit(mu2.uni))/1000
  convergence.uni3<-length(na.omit(mu3.uni))/1000
  convergence.uni4<-length(na.omit(mu4.uni))/1000
  
  ### summarize results
  if (delta1!=0){
    bias.mu=c(mean((na.omit(mu1)-delta1)/delta1),mean((na.omit(mu2)-delta2)/delta2),
              mean((na.omit(mu3)-delta3)/delta3),mean((na.omit(mu4)-delta4)/delta4),
              mean((na.omit(mu1.uni)-delta1)/delta1),mean((na.omit(mu2.uni)-delta2)/delta2),
              mean((na.omit(mu3.uni)-delta3)/delta3),mean((na.omit(mu4.uni)-delta4)/delta4))
  } else {
    bias.mu=c(mean((na.omit(mu1)-delta1)),mean((na.omit(mu2)-delta2)),
              mean((na.omit(mu3)-delta3)),mean((na.omit(mu4)-delta4)),
              mean((na.omit(mu1.uni)-delta1)),mean((na.omit(mu2.uni)-delta2)),
              mean((na.omit(mu3.uni)-delta3)),mean((na.omit(mu4.uni)-delta4)))}
  
  t.mu1.se<-sd(na.omit(mu1))
  t.mu2.se<-sd(na.omit(mu2))
  t.mu3.se<-sd(na.omit(mu3))
  t.mu4.se<-sd(na.omit(mu4))
  t.mu1.uni.se<-sd(na.omit(mu1.uni))
  t.mu2.uni.se<-sd(na.omit(mu2.uni))
  t.mu3.uni.se<-sd(na.omit(mu3.uni))
  t.mu4.uni.se<-sd(na.omit(mu4.uni))
  bias.muse=c(mean((na.omit(mu1.se)-t.mu1.se)/t.mu1.se),mean((na.omit(mu2.se)-t.mu2.se)/t.mu2.se),
              mean((na.omit(mu3.se)-t.mu3.se)/t.mu3.se),mean((na.omit(mu4.se)-t.mu4.se)/t.mu4.se),
              mean((na.omit(mu1.se.uni)-t.mu1.uni.se)/t.mu1.uni.se),mean((na.omit(mu2.se.uni)-t.mu2.uni.se)/t.mu2.uni.se),
              mean((na.omit(mu3.se.uni)-t.mu3.uni.se)/t.mu3.uni.se),mean((na.omit(mu4.se.uni)-t.mu4.uni.se)/t.mu4.uni.se))
  
  if (ttau1!=0){
    bias.tau=c(mean((na.omit(e.tau1)-ttau1)/ttau1),mean((na.omit(e.tau2)-ttau2)/ttau2),
               mean((na.omit(e.tau3)-ttau3)/ttau3),mean((na.omit(e.tau4)-ttau4)/ttau4),
               mean((na.omit(e.tau1.uni)-ttau1)/ttau1),mean((na.omit(e.tau2.uni)-ttau2)/ttau2),
               mean((na.omit(e.tau3.uni)-ttau3)/ttau3),mean((na.omit(e.tau4.uni)-ttau4)/ttau4))} else {
                 bias.tau=c(mean((na.omit(e.tau1)-ttau1)),mean((na.omit(e.tau2)-ttau2)),
                            mean((na.omit(e.tau3)-ttau3)),mean((na.omit(e.tau4)-ttau4)),
                            mean((na.omit(e.tau1.uni)-ttau1)),mean((na.omit(e.tau2.uni)-ttau2)),
                            mean((na.omit(e.tau3.uni)-ttau3)),mean((na.omit(e.tau4.uni)-ttau4)))}
  
  if (rlevel2!=0){
    bias.rho12<- mean((na.omit(rho12)-rlevel2)/rlevel2)
    bias.rho13<- mean((na.omit(rho13)-rlevel2)/rlevel2)
    bias.rho14<- mean((na.omit(rho14)-rlevel2)/rlevel2)
    bias.rho23<- mean((na.omit(rho23)-rlevel2)/rlevel2)
    bias.rho24<- mean((na.omit(rho24)-rlevel2)/rlevel2)
    bias.rho34<- mean((na.omit(rho34)-rlevel2)/rlevel2)
    } else {
      bias.rho12<- mean((na.omit(rho12)-rlevel2))
      bias.rho13<- mean((na.omit(rho13)-rlevel2))
      bias.rho14<- mean((na.omit(rho14)-rlevel2))
      bias.rho23<- mean((na.omit(rho23)-rlevel2))
      bias.rho24<- mean((na.omit(rho24)-rlevel2))
      bias.rho34<- mean((na.omit(rho34)-rlevel2))
    }
  
  coveragemu1<-mean(na.omit(ci1))
  coveragemu2<-mean(na.omit(ci2))
  coveragemu3<-mean(na.omit(ci3))
  coveragemu4<-mean(na.omit(ci4))
  coveragemu1.uni<-mean(na.omit(ci1.uni))
  coveragemu2.uni<-mean(na.omit(ci2.uni))
  coveragemu3.uni<-mean(na.omit(ci3.uni))
  coveragemu4.uni<-mean(na.omit(ci4.uni))
  
  rmsemu<-sqrt(c(mean((na.omit(mu1)-delta1)^2),mean((na.omit(mu2)-delta2)^2),
                 mean((na.omit(mu3)-delta3)^2),mean((na.omit(mu4)-delta4)^2),
                 mean((na.omit(mu1.uni)-delta1)^2),mean((na.omit(mu2.uni)-delta2)^2),
                 mean((na.omit(mu3.uni)-delta3)^2),mean((na.omit(mu4.uni)-delta4)^2)))
  
  rmsetau<-sqrt(c(mean((na.omit(e.tau1)-ttau1)^2),mean((na.omit(e.tau2)-ttau2)^2),
                  mean((na.omit(e.tau3)-ttau3)^2),mean((na.omit(e.tau4)-ttau4)^2),
                  mean((na.omit(e.tau1.uni)-ttau1)^2),mean((na.omit(e.tau2.uni)-ttau2)^2),
                  mean((na.omit(e.tau3.uni)-ttau3)^2),mean((na.omit(e.tau4.uni)-ttau4)^2)))
  
  powermu1<-mean(na.omit(mu1p)<=0.05)
  powermu2<-mean(na.omit(mu2p)<=0.05)
  powermu3<-mean(na.omit(mu3p)<=0.05)
  powermu4<-mean(na.omit(mu4p)<=0.05)
  powermu1.uni<-mean(na.omit(mu1p.uni)<=0.05)
  powermu2.uni<-mean(na.omit(mu2p.uni)<=0.05)
  powermu3.uni<-mean(na.omit(mu3p.uni)<=0.05)
  powermu4.uni<-mean(na.omit(mu4p.uni)<=0.05)
  powerQp<-mean(na.omit(Qp)<=0.05)
  powerQp1.uni<-mean(na.omit(Qp1.uni)<=0.05)
  powerQp2.uni<-mean(na.omit(Qp2.uni)<=0.05)
  powerQp3.uni<-mean(na.omit(Qp3.uni)<=0.05)
  powerQp4.uni<-mean(na.omit(Qp4.uni)<=0.05)
  
  result.table<-c(condition[index,],convergence.mul,convergence.uni1,convergence.uni2,
                  convergence.uni3,convergence.uni4,
                  bias.mu,bias.muse,bias.tau,
                  bias.rho12, bias.rho13,bias.rho14,bias.rho23, bias.rho24, bias.rho34, 
                  coveragemu1, coveragemu2,
                  coveragemu3, coveragemu4,
                  coveragemu1.uni, coveragemu2.uni,coveragemu3.uni, coveragemu4.uni,
                  rmsemu,rmsetau,
                  powermu1, powermu2, powermu3, powermu4, 
                  powermu1.uni, powermu2.uni,powermu3.uni, powermu4.uni,
                  powerQp, powerQp1.uni, powerQp2.uni, powerQp3.uni, powerQp4.uni)
  
  fn <- paste0("result/",index,"result.csv")
  write.table(t(result.table), file = fn,  sep = ",", col.names = FALSE,append = TRUE)
  
}

#####################################################
do.call(do_one, args)

