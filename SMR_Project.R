library(ISLR2)
data(Hitters)
# Replace NA with column mean 
Hitters_updated <- Hitters
for (col in names(Hitters_updated)) {
  if (is.numeric(Hitters_updated[[col]])) {
    col_mean <- mean(Hitters_updated[[col]], na.rm = TRUE)
    Hitters_updated[[col]][is.na(Hitters_updated[[col]])] <- col_mean
  }
}
############################################################################
#####################Best Subset selection#####################
library(leaps)
best.fitted.model=regsubsets(Salary~.,Hitters_updated, nvmax=19)
summary=summary(best.fitted.model)
summary
coef(best.fitted.model,1)
coef(best.fitted.model,2)
coef(best.fitted.model,3)
coef(best.fitted.model,4)
coef(best.fitted.model,5)
coef(best.fitted.model,6)
coef(best.fitted.model,7) 
coef(best.fitted.model,8)
coef(best.fitted.model,9)
coef(best.fitted.model,10)
coef(best.fitted.model,11)
coef(best.fitted.model,12)
coef(best.fitted.model,13)
coef(best.fitted.model,14)
coef(best.fitted.model,15)
coef(best.fitted.model,16)
coef(best.fitted.model,17)
coef(best.fitted.model,18)
coef(best.fitted.model,19)
names(summary)

##################### Best model based on R^2 Evaluation ###############
summary$adjr2
plot(summary$adjr2 , xlab = "Number of Variables",
     ylab = "Adjusted RSq")
best.r2 <- which.max(summary$adjr2)
coef(best.fitted.model,best.r2)
points(best.r2, summary$adjr2[best.r2], col = "red", cex = 2, pch = 20)

##################### Best model based on CP Evaluation ###############
summary$cp
plot(summary$cp, xlab = "Number of Variables",
     ylab = "CP Value")
best.cp <- which.min(summary$cp)
points(best.cp,summary$cp[best.cp],col='red',cex=2,pch=20)

##################### Best model based on BIC Evaluation ###############
summary$bic
plot(summary$bic, xlab = "Number of Variables",
     ylab = "BIC Value")
best.bic <- which.min(summary$bic)
points(best.bic,summary$bic[best.bic],col='red',cex=2,pch=20)

############################## Best model based on AIC Evaluation ###############

library(MASS)
full_model <- lm(Salary ~ ., data = Hitters_updated)
best_aic_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(best_aic_model)

######################################################
############################## Based ob best R^2 fitting linear model ###################

model_data <- model.matrix(Salary ~ ., data = Hitters_updated)
best.vars <- names(coef(best.fitted.model, best.r2))[-1]  
X_best <- model_data[, best.vars]
best.lm <- lm(Hitters_updated$Salary ~ X_best)
summary(best.lm)

################## Residual for the above linear model ###############

################## PRESS Residual ##################
x=model.matrix(best.lm)
PRESS_res=summary(best.lm)$res/(1-hat(x))
PRESS_res
plot(PRESS_res)
abline(h=0)

PRESS=sum(PRESS_res^2)
TSS=sum(anova(best.lm)$'Sum Sq')
predrsq=1-PRESS/(TSS)
predrsq

################## standard and student residual ##################
rstudent(best.lm)
plot(rstudent(best.lm),ylab="Studentized residual")
rstandard(best.lm)
plot(rstandard(best.lm),ylab="Standarized residual")

############################### Finding outliers ####################################

      ############################ Leverage point  #####################
                ############.     using hat values.    #####################
h=hatvalues(best.lm)
h
k=11
n=322
t=2*k/n
for (i in 1:322)
  if(h[i]>t)
  {print(h[i])}
num_outliers <- sum(h[1:322] > t)
cat("Number of high-leverage points in first 322 rows:", num_outliers, "\n")


############################# influential point ###################
#####################.      using COOK'S Distance.   ##################### 
cook=cooks.distance(best.lm)
cook
n=322
k=11
cutoff = 4/(n-k-1)
cutoff
for (i in 1:322)
  if(cook[i]>cutoff)
  {print(cook[i])} 
num_outliers <- sum(cook[1:322] > cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")

########################.     Using DFFITS.     ##################

dff=dffits(best.lm)
abs(dff)
n=322
k=11
cutoff = 2*sqrt(k/n)
cutoff
for (i in 1:322)
  if(abs(dff)[i]>cutoff)
  {print(abs(dff)[i])} 
num_outliers <- sum(abs(dff)[1:322]>cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")


########################     Using DFBEATs     ##################

dfb=dfbetas(best.lm)
abs(dfb)
n=322
k=11
cutoff=2/sqrt(n)
for (i in 1:322) 
  if(abs(dfb)[i] > cutoff)
  {print(abs(dfb)[i])}
num_outliers <- sum(abs(dfb)[1:322]>cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")


########################     Using CovRatio     ##################

cv=covratio(best.lm)
cv
n=322
k=11
for (i in 1:322)
  if(cv[i]>1+3*k/n|cv[i]<1-3*k/n)
  {print(cv[i])}
num_outliers <- sum(cv[1:322]>1+3*k/n|cv[1:322]<1-3*k/n)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")

####################### Durbin-watson test ##############
library(DescTools)
DurbinWatsonTest(best.lm, alternative ="less")
DurbinWatsonTest(best.lm, alternative ="greater")
DurbinWatsonTest(best.lm, alternative ="two.sided")

################ Checking Heteroscedasticity  #############

plot(best.lm$fit,best.lm$res,xlab="fitted", ylab="Residual")
abline(h=0)

########################Checking Normality #############################
########p-p Plot %%%%%%%%%%%%%
library(faraway)
par(mar=c(2,2,2,0))
probDist <- pnorm(rstudent(best.lm))
plot(ppoints(length(rstudent(best.lm))), sort(probDist),
     main = "PP Plot_studentized", xlab = "Observed Probability", 
     ylab = "Expected Probability") 
abline(0,1)

##################  Q-Q Plot  ###########

qqnorm(rstudent(best.lm))
qqline(rstudent(best.lm), col='red')
  

##################   BOX-COX Transformation on best subsset using R^2 #####################################
library(car)
max(Hitters_updated$Salary)/min(Hitters_updated$Salary)
best.vars <- names(coef(best.fitted.model, best.r2))[-1] 
model_data <- model.matrix(Salary ~ ., data = Hitters_updated)
X_best <- model_data[, best.vars]
p1 <- powerTransform(Salary ~ X_best, family = "bcPower", data = Hitters_updated)
summary(p1)

################Fitting linear model using Box Cox transform ######

lambda <- p1$lambda
Salary_bc <- bcPower(Hitters_updated$Salary, lambda, jacobian.adjusted = TRUE)
df_model <- as.data.frame(X_best)
df_model$Salary_bc <- Salary_bc
bc_model <- lm(Salary_bc ~ ., data = df_model)
summary(bc_model)


######## since BC value > 10 hence it is helpful

# Apply Box-Cox transformation to Salary using the fitted lambda


################## Residual for the above linear model ###############

################## PRESS Residual ##################
x=model.matrix(bc_model)
PRESS_res=summary(bc_model)$res/(1-hat(x))
PRESS_res
plot(PRESS_res)
abline(h=0)

PRESS=sum(PRESS_res^2)
TSS=sum(anova(bc_model)$'Sum Sq')
predrsq=1-PRESS/(TSS)
predrsq

################## standard and student residual ##################
rstudent(bc_model)
plot(rstudent(bc_model),ylab="Studentized residual")
rstandard(bc_model)
plot(rstandard(bc_model),ylab="Standarized residual")

############################### Finding outliers ####################################

############################ Leverage point  #####################
############.     using hat values.    #####################
h=hatvalues(bc_model)
h
k=11
n=322
t=2*k/n
for (i in 1:322)
  if(h[i]>t)
  {print(h[i])}
num_outliers <- sum(h[1:322] > t)
cat("Number of high-leverage points in first 322 rows:", num_outliers, "\n")


############################# influential point ###################
#####################.      using COOK'S Distance.   ##################### 
cook=cooks.distance(bc_model)
cook
n=322
k=11
cutoff = 4/(n-k-1)
cutoff
for (i in 1:322)
  if(cook[i]>cutoff)
  {print(cook[i])} 
num_outliers <- sum(cook[1:322] > cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")

########################.     Using DFFITS.     ##################

dff=dffits(bc_model)
abs(dff)
n=322
k=19
cutoff = 2*sqrt(k/n)
cutoff
for (i in 1:322)
  if(abs(dff)[i]>cutoff)
  {print(abs(dff)[i])} 
num_outliers <- sum(abs(dff)[1:322]>cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")


########################     Using DFBEATs     ##################

dfb=dfbetas(bc_model)
abs(dfb)
n=322
k=19
cutoff=2/sqrt(n)
for (i in 1:322) 
  if(abs(dfb)[i] > cutoff)
  {print(abs(dfb)[i])}
num_outliers <- sum(abs(dfb)[1:322]>cutoff)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")


########################     Using CovRatio     ##################

cv=covratio(bc_model)
cv
n=322
k=19
for (i in 1:322)
  if(cv[i]>1+3*k/n|cv[i]<1-3*k/n)
  {print(cv[i])}
num_outliers <- sum(cv[1:322]>1+3*k/n|cv[1:322]<1-3*k/n)
cat("Number of Influential points in first 322 rows:", num_outliers, "\n")

####################### Durbin-watson test ##############
library(DescTools)
DurbinWatsonTest(bc_model, alternative ="less")
DurbinWatsonTest(bc_model, alternative ="greater")
DurbinWatsonTest(bc_model, alternative ="two.sided")

################ Checking Heteroscedasticity  #############

plot(bc_model$fit,bc_model$res,xlab="fitted", ylab="Residual")
abline(h=0)

########################Checking Normality #############################
########p-p Plot %%%%%%%%%%%%%
library(faraway)
par(mar=c(2,2,2,0))
probDist <- pnorm(rstudent(bc_model))
plot(ppoints(length(rstudent(bc_model))), sort(probDist),
     main = "PP Plot_studentized", xlab = "Observed Probability", 
     ylab = "Expected Probability") 
abline(0,1)

##################  Q-Q Plot  ###########

qqnorm(rstudent(bc_model))
qqline(rstudent(bc_model), col='red')


######################################################################################
##################### Validation set approach ########################
set.seed(1)
train=sample(c(TRUE,FALSE), nrow(Hitters_updated),rep=T)
train
test=(!train)
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
train.mat
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
test.mat

##############best subset for train data
library(leaps)
regfit.best=regsubsets(Salary~.,data=Hitters_updated[train, ],nvmax=19, method="forward")
test.mse=rep(NA, 19)
for (i in 1:19){
  coefi=coef(regfit.best, id=i)
  pred=test.mat[ ,names(coefi)]%*%coefi
  test.mse[i]=mean((Hitters_updated$Salary[test]-pred)^2)
}
test.mse

############# Finding Best model using minumum validation error ######
which.min(test.mse)
coef(regfit.best, which.min(test.mse))##best model according to validation set
par(mfrow=c(1,1))
plot(test.mse)


###################### K-Fold CV approach ########################## 
##### Predictor function to used in calculating CV error #########      


predict.regsubsets=function(object,newdata,id){ 
  form=as.formula(object$call[[2]])
  mat=model.matrix(form,newdata)
  coefi=coef(object,id=id)
  xvars=names(coefi)
  mat[,xvars]%*%coefi
}
sm=summary(best.fit)
sm$obj$call[[3]]
######  Compute the test MSE using cross-validation set approach #############
k=10
n=nrow(Hitters_updated)
set.seed(1)
folds=sample(rep(1:k, length=n))
folds
cv.errors=matrix(NA,k, 19, dimnames = list(NULL,paste(1:19)))
cv.errors
#######
for(j in 1:k){ 
  best.fit=regsubsets(Salary~.,
                      data=Hitters_updated[folds!=j,],
                      nvmax=19)
  for(i in 1:19){ 
    pred=predict(best.fit,Hitters_updated[folds==j,],id=i) 
    cv.errors[j,i]=mean( (Hitters_updated$Salary[folds==j]-pred)^2)
  }
}
##############
mean.cv.errors=apply(cv.errors,2,mean)
mean.cv.errors
which.min(mean.cv.errors)
par(mfrow=c(1,1))
plot(mean.cv.errors,type="b")

#######Best Kfold ##########
reg.best=regsubsets(Salary~.,data=Hitters_updated, nvmax=19)
coef(reg.best,9)  #since  for 9 it is min




################################################################################
######################## Ridge Regression ####################
######### Mking grids and ridge model coefficints using any grid value ########
grid=10^seq(10,-2,length=100)
grid
library(glmnet)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)
ridge.mod$lambda[45] 
coef(ridge.mod)[, 45] 

##################################################################
#####  Estimation of the TEST MSE using Validation set ##########
###################################################################
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
grid=10^seq(10,-4,length.out=99)
grid=c(grid,0)
grid
set.seed(18)
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
test=(-train)
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
test.mat
y.test=y[test]
##Fit the regression model, for all lambda, for the train data
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=grid, thresh=1e-12)
plot(ridge.mod,xvar="lambda")

########## Test MSE for lambda
ridge.pred3=predict(ridge.mod,s=ridge.mod$lambda[46],newx=x[test,], exact=T) 
mean((ridge.pred3-y.test)^2)
predict(ridge.mod,s=ridge.mod$lambda[46],type="coefficients")[1:19,]


###########Cross validation error for Ridge Regression Model ################

set.seed(16)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
test=(-train)
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
y.test=y[test]

library(glmnet)
cv.out=cv.glmnet(x[train,],y[train],alpha=0, nfolds=5)
plot(cv.out)
#######  Choosing best Lambda  ###########

best.lambda=cv.out$lambda.min
best.lambda

##############  TEST MSE corresponding to best lambda
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=best.lambda, thresh=1e-12)
ridge.pred4=predict(ridge.mod,s=best.lambda,newx=x[test,]) 
mean((ridge.pred4-y.test)^2)




###########Fitting using best lambda
Final_ridge=glmnet(x,y,alpha=0)
predict(Final_ridge,type="coefficients",s=best.lambda)[1:19,]



######### THE LASSO ############################################



library(glmnet)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
grid=10^seq(10,-2,length=100)
lasso.mod=glmnet(x,y,alpha=1,lambda=grid, thresh=1e-12)
lasso.mod$lambda[95] 
dim(coef(lasso.mod))


#####  Estimation of the TEST MSE using Validation set ##########

set.seed(19)
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
test=(-train)
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
y.test=y[test]
lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=grid, thresh=1e-12)
plot(lasso.mod)
##### TEST MSE for lambda[1]##########
lasso.pred1=predict(lasso.mod,s=lasso.mod$lambda[94] ,newx=x[test,])
mean((lasso.pred1-y.test)^2)#Test MSE for lambda[94]
predict(lasso.mod,s=lasso.mod$lambda[94] ,type="coefficients")[1:19,]


###############################Analaysis based on Train Data ##########################################
set.seed(11)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
test=(-train)
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
y.test=y[test]

library(glmnet)
cv.out=cv.glmnet(x[train,],y[train],alpha=1, nfolds=5)
plot(cv.out)
best.lambda=cv.out$lambda.min
best.lambda
Lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=best.lambda, thresh=1e-12)
Lasso.pred4=predict(Lasso.mod,newx=x[test,]) 
mean((Lasso.pred4-y.test)^2)
########### Refit the model using the full data by using the lambda chosen
#by cross-validation technique
Full.mod=glmnet(x,y,alpha=1, lambda=best.lambda)
coef(Full.mod)



################ Choosing best lambda using Cross validation on full data  #################
set.seed(18)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
grid=10^seq(10,-2,length=100)
grid
cv.full=cv.glmnet(x,y,alpha=1, lambda=grid, nfolds=5)
plot(cv.full)
best.lambda2=cv.full$lambda.min
best.lambda2
lasso.full=glmnet(x,y,alpha=1, lambda=best.lambda2)
coef(lasso.full)

################  Fitting Elastic NET on a train data #####################
################  and finding test mse on the test data  ##################
###########################################################################
set.seed(12)
x=model.matrix(Salary~.,Hitters_updated)[,-1]
y=Hitters_updated$Salary
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(Salary~., data=Hitters_updated[train, ])
test=(-train)
test.mat=model.matrix(Salary~., data=Hitters_updated[test, ])
y.test=y[test]
list.fit=list() ####Creating empty list
for(i in 0:10){
  mod.alpha=paste("alpha", i/10)
  list.fit[[mod.alpha]]=cv.glmnet(x[train,],y[train],alpha=i/10, nfolds=5)
}
########## TEST MSE calculation ##########
All_mod=data.frame()
for(i in 0:10){
  mod.alpha=paste("alpha", i/10)
  pred.net= predict(list.fit[[mod.alpha]],
                    s=list.fit[[mod.alpha]]$lambda.min, newx=x[test,]) 
  mse.net=mean((y.test-pred.net)^2)
  mse.list=data.frame(alpha=i/10, mse=mse.net, model.name=mod.alpha)
  All_mod=rbind(All_mod,mse.list)
}
All_mod
min.mse=min(All_mod[,"mse"])
min.mse
All_mod[All_mod$mse==min.mse,]
####  Final Model in elastic net based on the full data ####################
Final_model=glmnet(x,y,alpha=0.2, lambda=list.fit[["alpha 0.2"]]$lambda.min)
coef(Final_model) 


#############################################

#############################################

#############################################

#####################Logistic Regression########################

full.model=glm(NewLeague~., data=Hitters_updated, family=binomial)
summary(full.model)
full.prob <- predict(full.model, Hitters_updated, type = "response")
full.prob
full.pred <- rep("A", nrow(Hitters_updated))
full.pred[full.prob > 0.5] <- "N"
full.pred <- factor(full.pred, levels = levels(Hitters_updated$NewLeague))
# Confusion matrix
table(Predicted = full.pred, Actual = Hitters_updated$NewLeague)
mean(full.pred == Hitters_updated$NewLeague)

#############################################

###  Fitting logistic model on Train set of the data ############
#####################################################################
train=sample(1:nrow(Hitters_updated), nrow(Hitters_updated)*0.70) 
train.mat=Hitters_updated[train, ]## Train data
test=(-train)
test.mat=Hitters_updated[test, ]## Test data
#### Fitting the model on train set
logis.fit=glm(NewLeague~.,data=Hitters_updated,
              family=binomial, subset=train)
summary(logis.fit)
pred.prob=predict(logis.fit,test.mat,type="response")
pred.prob
library(pROC)
roc_score=roc(test.mat[,20], pred.prob)
roc_score
par(pty="s")##Margin set
plot(roc_score ,main ="ROC curve -- Logistic Regression ", legacy.axes=T)
plot(roc_score ,main ="ROC curve -- Logistic Regression ", legacy.axes=T,
     xlab="False positive rate", ylab="True positive rate")

###############################################################
############## Subset Selection  ############################
#############################################################
library(mlbench)
#######  Fitting model using stepwise, forward and backward  ###############
full.model=glm(NewLeague~., data=Hitters_updated, family=binomial)
summary(full.model)
null.model=glm(NewLeague~1, data=Hitters_updated, family=binomial)
summary(null.model)
library(MASS)
stepAIC(full.model, direction="both")
stepAIC(full.model, direction="backward")
stepAIC(null.model, scope=list(lower=null.model, upper=full.model), direction="forward")

#######################





















##################  Shrinkage Method ##########%%%%%%%%%%%%
######### THE LASSO ############################################
#################################################################
##############################################################

set.seed(11)
x=model.matrix(NewLeague~.,Hitters_updated,)[,-1]
y=Hitters_updated$NewLeague
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(NewLeague~., data=Hitters_updated[train, ])## Train data without response
test=(-train)
x.test=Hitters_updated[test,] ##test data with response
x.test
test.mat=model.matrix(NewLeague~., data=x.test)[,-1]## Test data without response
test.mat ##Test data without response
y.test=y[test]
y.test ##Test responses
##cv.glmnet has some dafault lambda values; 
#we may also set our own lambda values through grid
library(glmnet)
cv.out=cv.glmnet(x[train,],y[train],alpha=1, family="binomial", nfolds=5)
plot(cv.out)
#######  Chossing best Lambda  ###########
best.lambda=cv.out$lambda.min
best.lambda
########## best Lasso model corresponding to min lambda #####
Lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=best.lambda, family="binomial")
coef(Lasso.mod)
##############  ROC for best model  #########
##For model (fitted with glmnet function) test data should be without response
pred.prob=predict(Lasso.mod,test.mat,type="response")##test matrix should be withouit data
pred.prob
library(pROC)
roc_dia=roc(y.test, pred.prob) ## ROC score
roc_dia
par(pty="s")##Margin set
plot(roc_dia ,main ="ROC curve -- Logistic Regression ", legacy.axes=T,
     xlab="False positive rate", ylab="True positive rate") ##ROC curve
########### Refit the model using the full data  ##########
Lasso.full=glmnet(x,y,alpha=1, lambda=best.lambda, family="binomial")
coef(Lasso.full)
###############################################
#######  Comparing it with full model (i.e., with MLE model)  ###########
#######################################################
full.mod <- glm(NewLeague~.,data=Hitters_updated,
                family = binomial,
                subset =train)
summary(full.mod)
##For model (fitted with glm function) test data should be with response
prob.full=predict(full.mod,x.test,type="response")
prob.full
library(pROC)
roc_full=roc(y.test, prob.full) ## ROC score
roc_full
par(pty="s")##Margin set
plot(roc_dia ,main ="ROC curve -- Logistic Regression ", legacy.axes=T,
     xlab="False positive rate", ylab="True positive rate") ##ROC curve




################  Fitting Elastic NET on a train data #####################
###########################################################################

set.seed(12)
x=model.matrix(NewLeague~.,Hitters_updated,)[,-1]
y=Hitters_updated$NewLeague
train=sample(1:nrow(x), nrow(x)*(2/3)) 
train.mat=model.matrix(NewLeague~., data=Hitters_updated[train, ])## Train data without response
test=(-train)
x.test=Hitters_updated[test,] ##test data with response
test.mat=model.matrix(NewLeague~., data=x.test)[,-1]## Test data without response
test.mat ##Test data without response
y.test=y[test] ##Test responses
##############################################################################
##5-fold cross validation based on train data; we may set our own lambda=grid
##For each alpha, we get the best model corresponding to best lambda 
list.fit=list() ####Creating empty list
for(i in 0:10){
  mod.alpha=paste("alpha", i/10)
  list.fit[[mod.alpha]]=cv.glmnet(x[train,],y[train],alpha=i/10, family="binomial", nfolds=5)
}
########## Binomial deviation calculation ##########
##Fitting the elastic net corresponding to each alpha (and its corresponding best lambda)
elastic.mod=list()
for(i in 0:10){
  mod.alpha=paste("alpha", i/10)
  elastic.mod[[mod.alpha]]=glmnet(x[train,],y[train],
                                  alpha=i/10,lambda=list.fit[[mod.alpha]]$lambda.min, family="binomial") 
}
###################### #######
##ROC measure corresponding to each alpha (and its corresponding best lambda)
roc.list=data.frame()
All_mod=data.frame()
for(i in 0:10){
  mod.alpha=paste("alpha", i/10)
  pred.net= predict(elastic.mod[[mod.alpha]], test.mat, type="response") 
  roc.net=auc(roc(y.test, pred.net))
  roc.list=data.frame(alpha=i/10, ROC=roc.net, model.name=mod.alpha)
  All_mod=rbind(All_mod,roc.list)
}
All_mod
max.ROC=max(All_mod[,"ROC"])
max.ROC
All_mod[All_mod$ROC==max.ROC,]
####  Final Model in elastic net based on the full data ####################
Final_model=glmnet(x,y,alpha=0.2, lambda=list.fit[["alpha 0.2"]]$lambda.min, family="binomial")
coef(Final_model)


############  Poisson regression ###############################
################################################################
## Poission regression on fulldata #######
mod.pois <- glm(bikers~ mnth + hr + workingday + temp + weathersit, 
                data = Bikeshare, family = poisson)
summary(mod.pois)
########################
plot(coef(mod.pois)[2:12], xlab="Month", ylab="coefficient",
     xaxt="n", col="blue", pch=19, type="o")
axis(side = 1, at = 1:12, labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"))
plot(coef(mod.pois)[1:23], xlab = "Hour", ylab = "Coefficient", 
     col = "blue", pch = 19, type = "o")


###############################

mod.pois <- glm(AtBat ~ League + Division + NewLeague + Hits + HmRun + Runs + 
                  RBI + Walks + Years + CAtBat + CHits + CHmRun + CRuns + 
                  CRBI + CWalks + PutOuts + Assists + Errors,
                data = Hitters_updated, family = poisson)

# 4. Model summary
summary(mod.pois)

# 5. Plot coefficients for the first 15 numeric predictors (adjust if needed)
coef_vals <- coef(mod.pois)[2:16]  # Skipping Intercept
plot(coef_vals, xlab = "Predictors", ylab = "Coefficient", 
     xaxt = "n", col = "blue", pch = 19, type = "o")
axis(side = 1, at = 1:length(coef_vals), labels = names(coef_vals), las=2, cex.axis=0.7)

######its result shows atbat highest at NewLeagueN






