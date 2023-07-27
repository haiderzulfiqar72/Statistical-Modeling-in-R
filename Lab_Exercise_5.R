setwd("E:/Study Material/Tampere - Grad/Studies/Year I/Period III/Statistical Modeling 1 - DATA.STAT.740/Lab Excercises/Lab Excercise 5")

#Question 1

data.retino<-read.table("retinopathy.txt", sep="\t", dec=".", header=TRUE)
attach(data.retino)
library(nnet)
library(mvtnorm)
library(ordinal)
library(nlme)

#a)
data.retino$RET <- as.factor(data.retino$RET)
model.retino<-multinom(RET~DIAB, data=data.retino)
summary(model.retino)

newdata<-data.frame(DIAB=20)
predict(model.retino, newdata=newdata, type="probs")

#b)
model.cum<- clm(RET~DIAB, data=data.retino)
summary(model.cum)

newdata<-data.frame(DIAB=20)
pred.cum<- predict(model.cum, newdata=newdata, type="prob")
data.frame(newdata,pred.cum)

#c)
model.2<- clm(RET~DIAB+GH+BP+SM, data=data.retino)
summary(model.2)

y0<-RET==0
y1<-RET==1
y2<-RET==2

Y<-cbind(as.numeric(y0),as.numeric(y1),as.numeric(y2))
         
ndata <-data.retino
ndata$RET<- NULL

fitted1<-predict(model.cum, newdata=ndata, type="prob")$fit
eM1<- Y-fitted1
MSE.1<- mean(diag(eM1%*%t(eM1)))

fitted2<-predict(model.2, newdata=ndata, type="prob")$fit
eC2<-Y-fitted2
MSE.2<- mean(diag(eC2%*%t(eC2)))

#Hence, Full Model is better

#Question 2
data.yield<-read.table("NitrogenYield.txt", sep="\t", dec=".", header=TRUE)
attach(data.yield)

#a)
model.poly2<-lm(Yield~Nitrogen+I(Nitrogen^2))
coef(model.poly2)

#b)
model.polyexp<- glm(Yield~log(Nitrogen), family=gaussian(link="log"), data=data.yield) 
summary(model.polyexp)

nd<-data.frame(Nitrogen=150)
predict(model.polyexp, newdata=nd, type="response")

#c)
model.asymp<-nls(Yield~SSasymp(Nitrogen, Asym,R0,lrc), data=data.yield)
coef(model.asymp)

#d)
model.ssm<- nls(Yield~SSmicmen(Nitrogen, Vm, K), data=data.yield)
summary(model.ssm)

predict(model.ssm, newdata=nd, type="probs")

#e)
model.gom<-nls(Yield~SSasymp(Nitrogen, Asym, b2, b3), data=data.yield)
summary(model.gom)


beta<- coef(model.gom)
cov.beta<- vcov(model.gom)

beta.star<-rmvnorm(100000, mean = beta, sigma = cov.beta)
newd<-data.frame(Nitrogen=150)

Asym<-beta.star[,1]
R0<-beta.star[,2]
lrc<-beta.star[,3]

mu.star<- Asym+(R0-Asym)*exp(-exp(lrc)*newd$Nitrogen)

sigma2<-sigma(model.gom)^2
yf.star<-rnorm(100000, mean=mu.star, sd=sqrt(sigma2))
pred.lowerbound<-quantile(yf.star, c(0.1))
pred.upperbound<-quantile(yf.star, c(0.9))
pred.lowerbound
pred.upperbound
