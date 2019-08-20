library(survival)
library(casebase)
library(glmnet)
#process for creating a manual absolute risk function
data=ERSPC


u=Surv(time = data$Follow.Up.Time, event = data$DeadOfPrCa)

xa=as.data.frame(data[,1])
xa$normal=rnorm(length(xa[,1]),mean=0,sd=6)
xa$wernormal=rnorm(length(xa[,1]),mean=5,sd=5)
x=as.matrix(xa)
coxfit=cv.glmnet(x=x,y=u, family="cox",alpha=0)

coxcoef=coef(coxfit)
betaHat=coxcoef@x
tab <- data.frame(table(data[data$DeadOfPrCa == 1, "Follow.Up.Time"])) 
y <- as.numeric(levels(tab[, 1]))[tab[, 1]] #ordered distinct event times
d <- tab[, 2] 


h0 <- rep(NA, length(y))
for(l in 1:length(y))
{
  h0[l] <- d[l] / sum(exp(x[data$Follow.Up.Time>=y[l],] %*% betaHat))
}
plot(h0,type="l")
