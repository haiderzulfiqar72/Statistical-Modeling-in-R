#Question 1

##a)

data<-read.table("ozone.txt", header=TRUE, sep="\t", dec=".")
data

model.identity<- glm(ozone~rad+temp, family=gaussian(link="identity"),data=data)

model.inverse<-glm(ozone~rad+temp+wind, family=gaussian(link="inverse"), data=data)
model.log<-glm(ozone~rad+temp+wind, family=gaussian(link="log"), data=data)
model.exp<-glm(ozone~log(rad)+log(temp)+log(wind), family=gaussian(link="log"), data=data)

AIC(model.identity, model.inverse, model.log, model.exp)

#Based on AIC values we get, logarithm model (model.log) seems to be the best fit.

##b)
model.n<-glm(ozone~rad+temp+wind, family=gaussian(link="log"), data=data)

residuals.n<-residuals(model.n, type="pearson")^2
mu.n<-fitted(model.n, type="response")
model.pearsonN<-lm(residuals.n~mu.n)
summary(model.pearsonN)

model.gamma<-glm(ozone~rad+temp+wind, family=Gamma(link="log"), data=data)
residuals.gamma<-residuals(model.gamma, type="pearson")^2
mu.gamma<-fitted(model.gamma, type="response")
model.pearsonG<-lm(residuals.gamma~mu.gamma)
summary(model.pearsonG)

model.invgamma<-glm(ozone~rad+temp+wind, family=inverse.gaussian(link="log"), data=data)
residuals.invgamma<-residuals(model.invgamma, type="pearson")^2
mu.invgamma<-fitted(model.invgamma, type="response")
model.pearsonIG<-lm(residuals.invgamma~mu.invgamma)
summary(model.pearsonIG)

plot(fitted(model.n, type="response"), residuals(model.n, type="response")^2)
plot(fitted(model.n, type="response"), residuals(model.n, type="pearson")^2)
plot(fitted(model.gamma, type="response"), residuals(model.gamma, type="pearson")^2)
plot(fitted(model.invgamma, type="response"), residuals(model.invgamma, type="pearson")^2)

#Inverse Gaussian is the most suitable distribution

##c)
windtemp = data$wind-data$temp
df2 <- cbind(data, windtemp)

model.f<- glm(ozone~rad+temp+wind, family=Gamma(link="log"),data=data)
model.H0<- glm(ozone~rad+windtemp, family=Gamma(link="log"),data=df2)

anova(model.f,model.H0, test='F')


#Question 2

##a)
data<-read.table("weld.txt", header=TRUE, sep="\t", dec=".")
data

model1<- glm(Strength~factor(Drying)+factor(Material),family= gaussian(link='identity'),data=data)
model2<- glm(Strength~factor(Drying)+factor(Material),family= gaussian(link='log'),data=data)
model3<- glm(Strength~factor(Drying)+factor(Material),family= gaussian(link='inverse'),data=data)
model4<- glm(Strength~factor(Drying)+factor(Material),family= Gamma(link='log'),data=data)
model5<- glm(Strength~factor(Drying)+factor(Material),family= Gamma(link='inverse'),data=data)
model6<- glm(Strength~factor(Drying)+factor(Material),family= Gamma(link='identity'),data=data)
model7<- glm(Strength~factor(Drying)+factor(Material),family=inverse.gaussian(link="inverse"),data=data)

data.frame(model= c("model1", "model2", "model3","model4","model5","model6","model7"), AIC= sapply(list(model1, model2, model3, model4, model5, model6,model7), AIC))

#Choosing model7 

betahat<- cbind(coef(model7))

X<- model.matrix(model7)
summary(model7)

x1<- c(1,0,0)      
x2<- c(1,1,1)     

pred<-(t(x2)-t(x1))%*%betahat
sigma2<-sigma(model7)^2

# predictive hypothesis testing
T<-pred/sqrt(sigma2*(2+(t(x2)-t(x1))%*%solve(t(X)%*%X)%*%(x2-x1)))
#Predictive hypothesis give the test statistic value as 0.1763087

d<-2*pt(abs(T),df=13, lower.tail = FALSE)
#d-value comes out to be 0.8627683

##b)
model7<- glm(Strength~factor(Drying)+factor(Material),family=inverse.gaussian(link="inverse"),data=data)

newdata<-expand.grid(Drying=c(0,1), Material=c(0,1))
X<-model.matrix(~factor(Drying)+factor(Material), data=newdata)

kt<- X[2,]-X[1,]

kt<-t(kt)
result<-glht(model7, linfct = kt)

K<-cbind(X[1,]-t(X[-1,,drop=FALSE]),
         X[2,]-t(X[-(1:2),,drop=FALSE]),
        X[3,]-t(X[-(1:3),,drop=FALSE]),
        X[4,]-t(X[-(1:4),,drop=FALSE]))

library(multcomp)
pairwise<-glht(model7, linfct = t(K))
summary(pairwise)

#test statistic value of -14.310 is the largest pairwise difference

# data.frame(t(K),summary(pairwise)$test$pvalues)
# summary(pairwise, Chisqtest())
# 
# model0<-glm(Strength~1, family=inverse.gaussian(link="inverse"), data=data)
# anova(model0,model7, test="F")

##c)
model.invg<- glm(Strength~factor(Drying)*factor(Material)*factor(Thickness)*factor(Angle)*factor(Opening)*factor(Preheating),family=inverse.gaussian(link="inverse"),data=data)
model.nidentity<- glm(Strength~factor(Drying)*factor(Material)*factor(Thickness)*factor(Angle)*factor(Opening)*factor(Preheating),family=gaussian(link="identity"),data=data)
model.nlog<- glm(Strength~factor(Drying)*factor(Material)*factor(Thickness)*factor(Angle)*factor(Opening)*factor(Preheating),family=gaussian(link="log"),data=data)
model.ninverse<- glm(Strength~factor(Drying)*factor(Material)*factor(Thickness)*factor(Angle)*factor(Opening)*factor(Preheating),family=gaussian(link="inverse"),data=data)

AIC(model.invg,model.nidentity,model.nlog,model.ninverse)
#Choosing inverse gaussian model further and performing stepwise regression

full<- model.invg
reduced<- glm(Strength~1, family=inverse.gaussian(link="inverse"), data=data)
step.lm <- step(reduced, scope= list(lower=reduced, upper =full), direction= "both")
summary(step.lm)

#Based on the results of stepwise regression, the best model seems to be following, Strength ~ factor(Material) + factor(Drying) + factor(Opening) + factor(Preheating) + factor(Material):factor(Opening) + factor(Opening):factor(Preheating) + factor(Drying):factor(Preheating)

model.chosen<- glm(Strength~factor(Material) + factor(Drying) + factor(Opening) + factor(Preheating) + factor(Material):factor(Opening) + 
                     factor(Opening):factor(Preheating) + factor(Drying):factor(Preheating), 
                   family = inverse.gaussian(link = "inverse"), data = data)

fitted(model.chosen)[1]


