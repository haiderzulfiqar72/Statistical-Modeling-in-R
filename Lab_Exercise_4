setwd("E:/Study Material/Tampere - Grad/Studies/Year I/Period III/Statistical Modeling 1 - DATA.STAT.740/Lab Excercises/Lab Excercise 4")

#Question 1

data.chromo<-read.table("chromoabnormal.txt", sep="\t", dec=".", header=TRUE)
attach(data.chromo)

#a)
model.chromo<-glm(ca~offset(log(cells))+doseamt*doserate, family=poisson(link="log"), data=data.chromo)
summary(model.chromo)

newdata<-data.frame(doseamt= 4, doserate= 0.75, cells= 64070) 
pred<-predict(model.chromo, newdata=newdata, type="response")

#b)
newd<-data.frame(doseamt= 4, doserate= 0.75, cells= 50000)
mu<-predict(model.chromo, newdata=newd, type="response")

ratio.estimate<-mu/newd$cells

model.matrix(model.chromo)

xf<-t(cbind(1,4,0.75,3))

Var.eYf<-pred*(1+pred*t(xf)%*%vcov(model.chromo)%*%xf)

lower.Yf<-pred-qnorm(0.9)*sqrt(Var.eYf)
upper.Yf<-pred+qnorm(0.9)*sqrt(Var.eYf)

Var.eZf<-((1/newd$cells)^2)*Var.eYf

lower.Zf<-ratio.estimate-qnorm(0.9)*sqrt(Var.eZf)
upper.Zf<-ratio.estimate+qnorm(0.9)*sqrt(Var.eZf)
lower.Zf
upper.Zf


#C)
model.Quasi<-glm(ca~offset(log(cells))+doseamt*doserate, family=quasipoisson(link="log"), data.chromo)
summary(model.Quasi)

model.H0Quasi<-glm(ca~offset(log(cells))+doseamt, family=quasipoisson(link="log"), data.chromo)
summary(model.H0Quasi)

anova(model.H0Quasi, model.Quasi, test="F")
anova(model.H0Quasi, model.Quasi, test="F")$F[2]

#D)
model1<- glm(ca~offset(log(cells))+doseamt*doserate, family=poisson(link="log"), data=data.chromo)
model2<- glm(ca~offset(log(cells))+doseamt*doserate, family=quasipoisson(link="log"), data.chromo)
library(MASS)
model3<- glm.nb(ca~offset(log(cells))+doseamt*doserate,link=log, data=data.chromo)

AIC(model1,model2,model3)
#Here AIC cannot be a good judge because it is always 0 for quasiposson models 

MSE.p<- mean(residuals(model1, type="response")^2)
MSE.qp<- mean(residuals(model2, type="response")^2)
MSE.np<- mean(residuals(model3, type="response")^2)

residuals.1<-residuals(model1, type="pearson")^2
mu.1<-fitted(model1, type="response")
model.pearson1<-lm(residuals.1~mu.1)
summary(model.pearson1)

residuals.2<-residuals(model2, type="pearson")^2
mu.2<-fitted(model2, type="response")
model.pearson2<-lm(residuals.2~mu.2)
summary(model.pearson2)

residuals.3<-residuals(model3, type="pearson")^3
mu.3<-fitted(model1, type="response")
model.pearson3<-lm(residuals.3~mu.3)
summary(model.pearson3)

#Results from pearson residuals are quite similar, taking MSE into account, choosing Poisson model as the best fit

#Question2
data.apple<-read.table("applejuiceCRA7152.txt", sep="\t", dec=".", header=TRUE)
attach(data.apple)

#a)

#As Growth is a binomial variable, so we will use binomial distribution, and so identity, log and inverse link won't be tested

model.l<- glm(Growth~pH+Nisin+Temperature+Brix, family= binomial(link = "logit"), data=data.apple)
model.ql<- glm(Growth~pH+Nisin+Temperature+Brix, family= quasibinomial(link = "logit"), data=data.apple)
summary(model.l)
summary(model.ql)

# Dispersion factor is ~0.9 which is less than 1 so we will use biniomial family.

model.prob<-glm(Growth~pH+Nisin+Temperature+Brix, family= binomial(link = "probit"), data=data.apple)
model.gumb<-glm(Growth~pH+Nisin+Temperature+Brix, family= binomial(link = "cloglog"), data=data.apple)
model.cauch<-glm(Growth~pH+Nisin+Temperature+Brix, family=binomial(link="cauchit"), data=data.apple)

summary(model.prob)
summary(model.gumb)
summary(model.cauch)

residuals.l<-residuals(model.l, type="pearson")^2
mu.l<-fitted(model.l, type="response")
model.pearsonl<-lm(residuals.l~mu.l)
summary(model.pearsonl)

residuals.p<-residuals(model.prob, type="pearson")^2
mu.p<-fitted(model.prob, type="response")
model.pearsonp<-lm(residuals.p~mu.p)
summary(model.pearsonp)

residuals.g<-residuals(model.gumb, type="pearson")^2
mu.g<-fitted(model.gumb, type="response")
model.pearsong<-lm(residuals.g~mu.g)
summary(model.pearsong)

residuals.c<-residuals(model.cauch, type="pearson")^2
mu.c<-fitted(model.cauch, type="response")
model.pearsonc<-lm(residuals.c~mu.c)
summary(model.pearsonc)

MSE.l<- mean(residuals(model.l, type="response")^2)
MSE.prob<- mean(residuals(model.prob, type="response")^2)
MSE.gumb<- mean(residuals(model.gumb, type="response")^2)
MSE.cauch<- mean(residuals(model.cauch, type="response")^2)

AIC(model.l, model.prob, model.gumb, model.cauch)

#Probit model will be used as it seems the better as per above few tests.

#b)
ndata<-data.frame(pH= 4.5, Nisin= 20, Temperature= 30, Brix= 17)
mle<-predict(model.prob, newdata=ndata, type="response")

#c)
eta<-predict(model.prob, newdata=ndata, type="link", se.fit=TRUE)

link.lowerbound<-eta$fit-qnorm(0.975)*eta$se.fit
link.upperbound<-eta$fit+qnorm(0.975)*eta$se.fit

pnorm(eta$fit-qnorm(0.975)*eta$se.fit)
pnorm(eta$fit+qnorm(0.975)*eta$se.fit)

#d)
mu.prob<-predict(model.prob, newdata=ndata, type="response")
YS.pred<- 100 * mu.prob


data<-data.apple
mu.hat<-predict(model.prob, newdata=data, type="response")
index<-dim(data)[1]
n<-100

e.b<-numeric()

for(b in 1:1000){
  
sum.yb<-numeric()
for(i in 1:index){
    
    sum.yb[i]<-sum(sample(0:1, 1, replace=TRUE,prob=c(1-mu.hat[i],mu.hat[i])))
    
}
  
model.c<-glm(sum.yb~pH+Nisin+Temperature+Brix, family=binomial(link="probit"), data=data)
ndata<-data.frame(pH= 4.5, Nisin= 20, Temperature= 30, Brix= 17)
mu.fB<-predict(model.c, newdata=ndata, type="response")

e.b[b]<-sum(YS.pred)-sum(mu.fB*100)
}

var.error<-var(e.b)

z<-qnorm(c(0.9))

lower.bound<-YS.pred-z*sqrt(var.error)
upper.bound<-YS.pred+z*sqrt(var.error)
