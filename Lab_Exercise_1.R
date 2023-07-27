setwd("E:/Study Material/Tampere - Grad/Studies/Year I/Period III/Statistical Modeling 1 - DATA.STAT.740/Lab Excercises/Lab Excercise 2")

#Question 1

data<-read.table("ratstime.txt", sep="\t", dec=".", header=TRUE)
attach(data)

model<-lm(time~factor(poison)+factor(treat), data=data)
summary(model)

#a)

newdata<-expand.grid(poison= c('I','II','III'), treat=c('A','B','C','D'))

Xm<-model.matrix(~factor(poison)+factor(treat), data=newdata)
betahat<- coef(model)

x1mean<-(Xm[1,]+Xm[2,]+Xm[3,])/3
x2mean<-(Xm[4,]+Xm[5,]+Xm[6,])/3

x1mean%*%betahat
x2mean%*%betahat

kt<- x1mean-x2mean
kt%*%betahat

K<-cbind(kt)

q<-1
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#Other Pairwise Average Difference
#A and C
x3mean<-(Xm[1,]+Xm[2,]+Xm[3,])/3
x4mean<-(Xm[7,]+Xm[8,]+Xm[9,])/3 
kt<- x3mean-x4mean
K<-cbind(kt)
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#A and D
x5mean<-(Xm[1,]+Xm[2,]+Xm[3,])/3
x6mean<-(Xm[10,]+Xm[11,]+Xm[12,])/3 
kt<- x5mean-x6mean
K<-cbind(kt)
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#B and C
x7mean<-(Xm[4,]+Xm[5,]+Xm[6,])/3
x8mean<-(Xm[7,]+Xm[8,]+Xm[9,])/3 
kt<- x7mean-x8mean
K<-cbind(kt)
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#B and D
x9mean<-(Xm[4,]+Xm[5,]+Xm[6,])/3
x10mean<-(Xm[10,]+Xm[11,]+Xm[12,])/3 
kt<- x9mean-x10mean
K<-cbind(kt)
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#C and D
x11mean<-(Xm[7,]+Xm[8,]+Xm[9,])/3
x12mean<-(Xm[10,]+Xm[11,]+Xm[12,])/3 
kt<- x11mean-x12mean
K<-cbind(kt)
Wald<-(t(t(K)%*%betahat)%*%solve(t(K)%*%vcov(model)%*%K)%*%t(K)%*%betahat)/q
Wald
p.value<-pf(Wald, 1, 42, lower.tail = FALSE)
p.value

#b)

pred<-kt%*%betahat
sigma2<-sigma(model)^2

# predictive hypothesis testing

T<-pred/sqrt(sigma2*(2+(kt)%*%solve(t(model.matrix(model))%*%model.matrix(model))%*%(kt)))

d.value<-2*pt(abs(T),df=42, lower.tail = FALSE)
d.value


#c)

model.identity<- glm(time~factor(poison)+factor(treat), family=gaussian(link="identity"), data=data)
model.log<-glm(time~factor(poison)+factor(treat), family=gaussian(link="log"), data=data)
model.inverse<-glm(time~factor(poison)+factor(treat), family=gaussian(link="inverse"), data=data)
coef(model.identity)
coef(model.log)
coef(model.inverse)

# MSE values

mean(residuals(model.identity, type="response")^2)
mean(residuals(model.log, type="response")^2)
mean(residuals(model.inverse, type="response")^2)

#Inverse looks the best model wrt to smaller mse

#Question 2

#a)

data1<-read.table("Alba.txt", sep="\t", dec=".", header=TRUE)
attach(data1)

model<-glm(DryMatter~factor(Herbicide)+Dose+ factor(Herbicide):Dose,family=Gamma(link="inverse"),data=data1)
summary(model)

newdata<-data.frame(Dose=c(50), Herbicide="Glyphosate")
mu<- predict(model, newdata= newdata, type="response")

#b)
newdata<-data.frame(Dose=c(50), Herbicide="Glyphosate")
pred<-predict(model, newdata=newdata, type="response")
pred

model.matrix(model)

xf<-c(1,1,50,50)
xf


phi<-summary(model)$dispersion
Var.Yf<-phi
D.f<- (-1/pred^2)
Var.ef<-Var.Yf+(D.f^2)*t(xf)%*%vcov(model)%*%xf

lower.yf<-pred-qnorm(0.9)*sqrt(Var.ef)
upper.yf<-pred+qnorm(0.9)*sqrt(Var.ef)

lower.yf
upper.yf


#c)

modelci<-glm(DryMatter~factor(Herbicide)+Dose+ factor(Herbicide):Dose,family=inverse.gaussian(link="1/mu^2"),data=data1)
summary(modelci)

newdata<-data.frame(Dose=c(50), Herbicide="Glyphosate")

eta<- predict(modelci,newdata=newdata, type='link', se.fit=TRUE)
link.lowerbound<- eta$fit-qnorm(0.975)*eta$se.fit
link.upperbound<- eta$fit+qnorm(0.975)*eta$se.fit

lower.mu<-sqrt(1/link.lowerbound)
upper.mu<-sqrt(1/link.upperbound)

lower.mu
upper.mu
