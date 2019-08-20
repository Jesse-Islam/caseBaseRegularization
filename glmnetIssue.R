library(casebase)
library(glmnet)
#devtools::install_github("tidyverse/reprex")
library(reprex)
library(h2o)
h2o.init()




#####################################################################################
#####################################################################################

#There appears to be an issue with glmnet running at all, for any dataset.

#when either fitSmoothHazard (or what I need, fitSmoothHazard.fit)
# is used with family=c("glmnet") I run into this issue where the
# code never finishes running... For example: 
######Following Code will infinite loop######
#y=as.matrix(ERSPC[,c(2,3)])
#x=as.matrix(ERSPC[,c(1)])
#glmnetFit=fitSmoothHazard.fit(x,y,time="Follow.Up.Time",event="DeadOfPrCa",ratio=10,family=c("glmnet"),alpha=0.5)

#If you try to run the code above, it ends up never finishing.

#using just glm works fine:
y=as.matrix(ERSPC[,c(2,3)])
x=as.matrix(ERSPC[,c(1)])
glmFit=fitSmoothHazard.fit(x,y,time="Follow.Up.Time",event="DeadOfPrCa",ratio=10)

#When I dissect the main steps myself (I don't think I'm skipping a crucial step)
cbERSPC=sampleCaseBase(ERSPC,event="DeadOfPrCa",time="Follow.Up.Time",ratio=10)
y=as.matrix(cbERSPC[,c(3)])
x=as.matrix(cbERSPC[,c(1,2,4)])
cvGlmFit=cv.glmnet(x,y,family=c("binomial"))


#In that case, I'm leaving the offset as a variable. If I were to separate the offset:
######following code will infinite loop######
x=as.matrix(cbERSPC[,c(1,2)])
offset=as.matrix(cbERSPC[,c(4)])
#cvGlmFitOffset=cv.glmnet(x,y,family=c("binomial"),offset=offset)
#With glmfit,same issue:
#glmFitOffset=glmnet(x,y,family=c("binomial"),offset=x[,3])

#So there appears to be something off with the current cv.glmnet function's implementation
#should we include the offset parameter.


# It appears H2o's implementation works as expected with the offset:
cbERSPC.Split <- h2o.splitFrame(data = as.h2o(cbERSPC),ratios = 0.8, seed = 1234)
train <- cbERSPC.Split[[1]]
valid <- cbERSPC.Split[[2]]
x_names=colnames(x)
colnames(y)=c("DeadOfPrCa")
y_names=colnames(y)
h2oGlmFit=h2o.glm(x_names,y_names,training_frame = train,validation_frame=valid,nfolds=10,family=c("binomial"),solver=c("L_BFGS"),alpha=0.5,offset_column ="offset" )
#the only issue being that I need to reformat their output to match what absoluteRisk()
#Is expecting

#For Monday, I'll just use the offset as a variable. Like the following:
cbERSPC=sampleCaseBase(ERSPC,event="DeadOfPrCa",time="Follow.Up.Time",ratio=100)
y=as.matrix(cbERSPC[,c(3)])
x=as.matrix(cbERSPC[,c(1,2,4)])
cvGlmFitThatIWillUse=cv.glmnet(x,y,family=c("binomial"))

sessionInfo()