rm(list=ls())
library(caret)
library(evtree)
library(credittree)
library(proftree)
library(imbalance)
library(EMP)
library(timeR)
library(foreach)
library(doParallel)

UK <- read.csv2("UK.csv")
UK$ï..Censor <- ifelse(UK$ï..Censor == 0, 1, 0)
UK$ï..Censor <- as.factor(UK$ï..Censor)
UK$Censore <- NULL
UK$Open <- NULL
UK$Custgend <- NULL
UK$ï..Censor <- as.factor(UK$ï..Censor)

Results <- data.frame(Seed=integer(),
                      NrOfStartTrees=integer(),
                      EMP=factor(),
                      stringsAsFactors=FALSE)
fullUK <- UK

test <- UK[UK$ï..Censor == 1, ]
test <- test[c(1:4800),]
UK <- rbind(test,UK[UK$ï..Censor ==0,])
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
a <- foreach(i= 5:6, combine = 'c') %dopar%{
  seed <- i
  j <- 50#,150,200,250,300,500,1000,2000
    print(j)
    set.seed(seed)
    train_indices <- caret::createDataPartition(y = UK$ï..Censor, times = 1, p = 0.80, list = FALSE)
    train<- UK[train_indices,]
    test <- UK[-train_indices, ]
    set.seed(seed)
    model <- credittree::credittree(formula = train$ï..Censor ~ ., #Age + Amount + Curradd + Curremp + Depchild + Freqpaid + Insprem + Loantype + Marstat + Term + Homeowns + Purpose
                        data = train,
                        control = credittree::credittree.control(minbucket = 3,
                                                     minsplit = 10,
                                                     maxdepth = 5,
                                                     miniterations = 1,
                                                     niterations = 10,
                                                     ntrees = j,
                                                     verbose = F,
                                                     errortol = -0.5,))
    print(i)
    scores <- predict(model, newdata = test, type="prob")
    #predicted_labels <- predict
    
    #Results <- rbind(Results,c(seed=i,NrOfIterations = j,empCreditScoring(scores = scores[,2], classes = test$ï..Censor )))
    c(seed=i,NrOfIterations = j,EMP::empCreditScoring(scores = scores[,2], classes = test$ï..Censor ))
}
parallel::stopCluster(cl)