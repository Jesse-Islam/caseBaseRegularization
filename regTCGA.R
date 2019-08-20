#these packages were taken from github on 23 July 2019
library(casebase)
library(TCGA2STAT)
#these packages were pulled from bioconductor
library(BiocManager)
library(glmnet)
library(doParallel)



fitSmoothHazard.fitted <- function(x, y, formula_time, time, event, family = c("glm", "gbm", "glmnet"),
                                   censored.indicator, ratio = 100, ...) {
  family <- match.arg(family)
  if (family == "gam") stop("The matrix interface is not available for gam")
  if (family == "gbm" && !requireNamespace("gbm", quietly = TRUE)) {
    stop("Pkg gbm needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (family == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Pkg glmnet needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Default to linear term
  if (missing(formula_time)) {
    formula_time <- formula(paste("~", time))
    timeVar <- time
  } else {
    timeVar <- if (length(formula_time) == 3) all.vars(formula_time[[3]]) else all.vars(formula_time)
  }
  # There should only be one time variable
  stopifnot(length(timeVar) == 1)
  
  # Try to infer event from
  if (missing(event)) {
    varNames <- checkArgsTimeEvent(data = as.data.frame(y), time = timeVar)
    eventVar <- varNames$event
  } else eventVar <- event
  
  typeEvents <- sort(unique(y[,eventVar]))
  # Call sampleCaseBase
  originalData <- list("x" = x,
                       "y" = y)
  class(originalData) <- c(class(originalData), "data.fit")
  if (missing(censored.indicator)) {
    sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
                                 timeVar, eventVar,
                                 comprisk = (length(typeEvents) > 2),
                                 ratio)
  } else {
    sampleData <- sampleCaseBase(as.data.frame(cbind(y, x)),
                                 timeVar, eventVar,
                                 comprisk = (length(typeEvents) > 2),
                                 censored.indicator, ratio)
  }
  sample_event <- as.matrix(sampleData[,eventVar])
  sample_time_x <- cbind(as.matrix(sampleData[,!names(sampleData) %in% c(eventVar, timeVar, "offset")]),
                         model.matrix(update(formula_time, ~ . -1), sampleData))
  sample_offset <- sampleData$offset
  
  # Fit a binomial model if there are no competing risks
  if (length(typeEvents) == 2) {
    out <- switch(family,
                  "glm" = glm.fit(sample_time_x, sample_event,
                                  family = binomial(),
                                  offset = sample_offset),
                  "glmnet" = glmnet::cv.glmnet(sample_time_x, sample_event,
                                               family = "binomial",
                                               ...),
                  "gbm" = gbm::gbm.fit(sample_time_x, sample_event,
                                       distribution = "bernoulli",
                                       offset = sample_offset,
                                       verbose = FALSE, ...))
    
    out$originalData <- originalData
    out$typeEvents <- typeEvents
    out$timeVar <- timeVar
    out$eventVar <- eventVar
    out$matrix.fit <- TRUE
    out$formula_time <- formula_time
    out$offset<- sample_offset
  } else {
    stop("Not implemented yet")
    # Otherwise fit a multinomial regression
    # withCallingHandlers(model <- vglm(formula, family = multinomial(refLevel = 1),
    #                                   data = sampleData),
    #                     warning = handler_fitter)
    #
    # out <- new("CompRisk", model,
    #            originalData = originalData,
    #            typeEvents = typeEvents,
    #            timeVar = timeVar,
    #            eventVar = eventVar)
  }
  return(out)
}

#lusc.rnaseq2 <- getTCGA(disease="LUSC", data.type="RNASeq2", clinical=TRUE)
lusc.rnaseq2 <-readRDS('lusc.rnaseq2.rds')
highDimSurvData=na.omit(lusc.rnaseq2$merged.dat)
highDimNames=colnames(highDimSurvData)
fmla=as.formula(paste("status~ bs(OS) +",paste(highDimNames[20200:length(highDimNames)],collapse = "+")))
#highDimSurvData[, -c(1:3)]

cbTCGA=sampleCaseBase(highDimSurvData,event="status",time="OS",ratio=10)
y=as.matrix(cbTCGA[,c(2)])
x=as.matrix(cbTCGA[,c(3:100,20505)])
timeTest=cv.glmnet(x,y,family=c("binomial"))

new_data=as.data.frame(t(x[5,]))
ab=absoluteRisk(timeTest,time = seq(0,300, 1),newdata = new_data)
#hard coded fixes
timeTest$matrix.fit=1
#defaulting to matrix version, for cv.glmnet makes sense, as it doesn't have a native formula interface


y=as.matrix(highDimSurvData[,c(2,3)])
x=as.matrix(highDimSurvData[,c(4:length(highDimSurvData[1,]))])
uppers=rep(Inf, each=length(x[1,]))
lowers=rep(-Inf,each=length(x[1,]))
uppers=c(uppers,)
lowers=c(lowers,0)
registerDoParallel(2)
fit=fitSmoothHazard.fit(x,y,time="OS",event="status",family=c("glmnet"),ratio=10,lower.limits=lowers,upper.limits=uppers,parallel=TRUE,nfold=3)
wholeFit=fitSmoothHazard.fit(x,y,time="OS",event="status",family=c("glmnet"),ratio=10,parallel=FALSE,nfold=8,alpha=0)

new_data=as.data.frame(t(x[5,]))
#new_data$offset=fit$offset[1]
ab=absoluteRisk(wholeFit,time = seq(0,10000, 100),newdata = new_data)
plot(ab,type = "l")

