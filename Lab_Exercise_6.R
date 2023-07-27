setwd("E:/Study Material/Tampere - Grad/Studies/Year I/Period III/Statistical Modeling 1 - DATA.STAT.740/Lab Excercises/Lab Excercise 6")

#Question 1

data.tire<-read.table("tirereliability.txt", sep="\t", dec=".", header=TRUE)
attach(data.tire)
library(survival)
library(eha)

#a) 
model.cox<- coxph(Surv(survival, complete)~wedge, data=data.tire)
coef(model.cox)

#b)
newdata<-data.frame(wedge=c(0.6))
sf<-survfit(model.cox, newdata=newdata)
summary(sf, times=1)
summary(sf, times=1)$surv


plot(sf)

#c)
newdata1<-data.frame(wedge=c(0.6,1.6))
risk<-predict(model.cox, newdata=newdata1, type="risk")
risk[1]/risk[2]

#d)
model.full<-coxph(Surv(survival, complete)~wedge+peelForce+interBelt+wedge:peelForce, data=data.tire)
summary(model.full)

model.H0<-coxph(Surv(survival, complete)~peelForce+interBelt, data=data.tire)
summary(model.H0)
anova(model.H0, model.full)

#e)
newdata2<-data.frame(wedge=c(0.6), peelForce=c(0.8), interBelt=c(0.7))
sf1<-survfit(model.full, newdata=newdata2,conf.type="plain")
summary(sf1, times=1)
summary(sf1, times=1)$lower
summary(sf1, times=1)$upper
plot(sf,conf.int=c("none"))
plot(sf)
bhest <- basehaz(model.full)


#Question 2 

#a)
model.wph<-phreg(Surv(survival, complete)~wedge, data=data.tire, dist="weibull")
summary(model.wph)
coef(model.wph)

plot(model.wph, fn = "sur")

p<-exp(coef(model.wph)[3])
lambda<-exp(coef(model.wph)[2])
beta<-coef(model.wph)[1]

x1<-0.6
x2<- 1.6
t<- 100

hazard<- ((p/lambda)*(t/lambda)^(p-1) * exp(beta*x1)) / ((p/lambda)*(t/lambda)^(p-1) * exp(beta*x2))
hazard

exp(beta*(x1-x2))

#b)
lambda.star<-lambda/exp((x2*beta)/p)
mu<-lambda.star*gamma(1+(1/p))
mu

#c)
t.star<-rweibull(10000, shape=p, scale=lambda.star)
lowerbound<-quantile(t.star, c(0.1))
upperbound<-quantile(t.star, c(0.9))
lowerbound
upperbound

# qweibull(0.1, shape=p, scale=lambda.star)
# qweibull(0.9, shape=p, scale=lambda.star)

#d) 
model.complex<- phreg(Surv(survival, complete)~wedge+peelForce+interBelt+wedge*peelForce, data=data.tire, dist="weibull")
coef(model.complex)

p_new<-exp(coef(model.complex)[6])
beta_new<-coef(model.complex)[1:4]
lambda_new<-exp(coef(model.complex)[5])

x_new<-c(wedge=c(0.6),peelForce=c(0.8), interBelt=c(0.7), interse= c(0.6*0.8))
lambda.star_n<-lambda_new/exp((x_new%*%beta_new)/p_new)

survival1<-1-pweibull(1.00,shape=p_new, scale=lambda.star_n)
survival1
