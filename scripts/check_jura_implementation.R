source('revision/R/jura_common.R')
source('revision/R/jura_profile_estimator.R')
source('revision/R/jura_profile_fit.R')
datasets <- list(jura=readRDS('revision/data/jura-data.rds'))
bundle <- readRDS('revision/results/application-jura.rds')
checks <- list()
for(dataset in names(datasets)) {
  target <- datasets[[dataset]]$primary_target
  ans <- Filter(function(x) x$results$target[1]==target && x$results$repeat_id[1]==1L,bundle$details)[[1]]
  for(model in c('ordinary','trend')) {
    fit <- ans$full_models[[model]]
    tr <- data.frame(x=fit$xy[,1],y=fit$xy[,2],
                     response=datasets[[dataset]]$data$response[ans$geometry$train_rows])
    te <- data.frame(x=ans$geometry$target_xy[,1],y=ans$geometry$target_xy[,2])
    dep <- fit$dependence
    formula <- if(model=='ordinary') response~1 else response~x+y
    ref <- suppressWarnings(gstat::krige(formula,locations=~x+y,data=tr,newdata=te,
      model=gstat::vgm(dep$psill,'Exp',dep$range,dep$nugget),debug.level=0))$var1.pred
    actual <- predict_spatial(fit,as.matrix(te))
    error <- max(abs(ref-actual))
    stopifnot(error<1e-7*max(1,max(abs(ref))))
    checks[[length(checks)+1]] <- data.frame(check='independent gstat predictions',dataset,model,max_error=error)
  }
  # Change only target outcomes by a deliberately large amount. Neither model,
  # covariance diagnostics, CV predictions nor target predictions may change.
  perturbed <- datasets
  ix <- ans$geometry$target_rows
  perturbed[[dataset]]$data$response[ix] <- 1000+seq_along(ix)*100
  again <- analyse_application(dataset,target,1L,perturbed)
  p <- ans$predictions; q <- again$predictions
  stopifnot(identical(p$predicted,q$predicted),
    identical(ans$results$cv_mse,again$results$cv_mse),
    identical(ans$results$J2,again$results$J2),
    identical(ans$results$W,again$results$W),
    all(again$results$target_mse>ans$results$target_mse))
  checks[[length(checks)+1]] <- data.frame(check='target outcome exclusion',dataset,model='both',max_error=0)
}
# Recompute the original package-fitter configuration as a separate sensitivity check.
source('revision/R/jura_common.R')
original <- analyse_application('jura','published',1L,datasets)$results
reference <- subset(read.csv('revision/results/application-jura-original.csv'),
                    target=='published' & repeat_id==1L)
reference <- reference[match(paste(original$model,original$scheme),paste(reference$model,reference$scheme)),]
columns <- c('A','W','J1','J2','U1','U2','cv_mse','target_mse')
error <- max(abs(as.matrix(original[columns])-as.matrix(reference[columns])))
stopifnot(is.finite(error), error<1e-7)
checks[[length(checks)+1L]] <- data.frame(check='original fitter recomputation',dataset='jura',model='both',max_error=error)
write.csv(do.call(rbind,checks),'revision/results/supplementary/jura/implementation-checks.csv',row.names=FALSE)
print(do.call(rbind,checks))
