#Q1

data<-read.table("paper.txt", sep="\t", dec=".", header=TRUE)
attach(data)

model<-lm(strength~hardwood+pressure)

#a)
coef(model)
#Maximum likelihood estimate of B1 is -0.175595  

#b)
sigma(model)^2
#unbisaed estimate for sigma^2 is 1.328534

#c)
fitted(model)[36]
#fitted value for last observation i.e. 36 is 198.3235 

#d)
new.data<- data.frame(hardwood=c(8),pressure=c(550))
mu<- predict(model, newdata= new.data)
#MLE for expected value under the given data points is 197.6835 

#e)
mu<- predict(model, newdata= new.data, interval = 'prediction', level= .80)
#Estimate for the lower bound is 196.1169

#f)
model.H0<- lm(strength~hardwood)
anova(model,model.H0, test='F')
  #Since P-value here is 0.001674 which is < 0.05, we would reject the null hypothesis, i.e. it is better to include B1 (pressure)


#Q2

data1<-read.table("makiwaraboard.txt", sep="\t", dec=".", header=TRUE)
attach(data1)

model1<- lm(Deflection~factor(WoodType)+factor(BoardType))

#a)
exp.data<- data.frame(WoodType= "Oak", BoardType= "Tapered")
mle<- predict(model1, newdata= exp.data)
#mle under given scenario is 51.3127

#b)
model12<- lm(Deflection~factor(WoodType)+factor(BoardType)+factor(WoodType):factor(BoardType))
betahat<- coef(model12)

k1<- c(0,0,0,0,0,1,0,0)
k2<- c(0,0,0,0,0,0,1,0)
k3<- c(0,0,0,0,0,0,0,1)

K<- cbind(k1, k2, k3)
q=3

Wald<- (t(t(K)%*%betahat) %*% solve(t(K)%*%vcov(model12)%*%K)%*%t(K)%*%betahat)/q
#Wald Test Statistic gives the value of 0.1522429

p.value<- pf(Wald, 3, 328, lower.tail= FALSE)
#p-value for the wald test statistc is 0.9282088

anova(model1,model12)
#Verified as well that P-value is same coming from Wald test as well as from Anova test

#c)
betahat1<- cbind(coef(model12))

X<- model.matrix(model12)

x1<- c(1,1,0,0,0,0,0,0)      
x2<- c(1,0,0,1,1,0,0,1)     

pred<-(t(x2)-t(x1))%*%betahat1
sigma2<-sigma(model12)^2

# predictive hypothesis testing
T<-pred/sqrt(sigma2*(2+(t(x2)-t(x1))%*%solve(t(X)%*%X)%*%(x2-x1)))
#Predictive hypothesis give the test statistic value as -0.710013

d<-2*pt(abs(T),df=328, lower.tail = FALSE)
#d-value comes out to be 0.398245


