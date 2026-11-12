source('revision/R/local_libraries.R')
source('revision/R/common.R')
pkgload::load_all('package', quiet=TRUE)
# Jura-only manuscript analysis; routines preserved from the checked application comparison.

# Unknown-mean Gaussian prediction: ordinary or universal kriging. Every call
# receives only its own training observations, including covariance estimation.
fit_spatial <- function(xy, response, model) {
  xy <- as.matrix(xy)
  center <- colMeans(xy); scale <- apply(xy,2,sd)
  stopifnot(all(scale>0), length(response)==nrow(xy), nrow(xy)>=8)
  design <- function(z) {
    if(model=='ordinary') return(matrix(1,nrow(z),1,dimnames=list(NULL,'Intercept')))
    cbind(Intercept=1, sweep(sweep(as.matrix(z),2,center,'-'),2,scale,'/'))
  }
  F <- design(xy)
  z <- if(model=='ordinary') response else lm.fit(F,response)$residuals
  messages <- character()
  dep <- withCallingHandlers(spLeakage::estimate_dependence(
    data.frame(x=xy[,1],y=xy[,2],response=z),'response',coords=c('x','y')),
    warning=function(w) { messages <<- c(messages,conditionMessage(w)); invokeRestart('muffleWarning') })
  variance <- dep$psill+dep$nugget
  stopifnot(is.finite(variance),variance>0,is.finite(dep$range),dep$range>0)
  D <- pair_distance(xy)
  S <- dep$psill*exp(-D/dep$range)
  diag(S) <- variance
  jitter <- 0
  U <- tryCatch(chol(S),error=function(e) NULL)
  if(is.null(U)) {
    jitter <- variance*1e-10
    U <- chol(S+diag(jitter,nrow(S)))
  }
  solve_S <- function(b) backsolve(U,forwardsolve(t(U),b))
  SiF <- solve_S(F)
  information <- crossprod(F,SiF)
  beta <- as.numeric(solve(information,crossprod(F,solve_S(response))))
  alpha <- as.numeric(solve_S(response-as.numeric(F%*%beta)))
  vg <- dep$variogram
  fitted_vg <- dep$nugget+dep$psill*(1-exp(-vg$dist/dep$range))
  fit_record <- data.frame(model=model,n=nrow(xy),range_km=dep$range,nugget=dep$nugget,
    psill=dep$psill,signal=dep$signal_prop,variance=variance,
    zero_nugget=dep$nugget/variance<1e-7,negligible_signal=dep$signal_prop<1e-6,
    long_range=dep$range>=4.999*max(vg$dist),
    normalised_variogram_rmse=sqrt(weighted.mean((vg$gamma-fitted_vg)^2,vg$np))/variance,
    jitter=jitter,warnings=paste(unique(messages),collapse=' | '))
  list(xy=xy,model=model,center=center,scale=scale,beta=beta,alpha=alpha,
       dependence=dep,record=fit_record)
}

predict_spatial <- function(fit,xy) {
  xy <- as.matrix(xy)
  F <- if(fit$model=='ordinary') matrix(1,nrow(xy),1) else
    cbind(1,sweep(sweep(xy,2,fit$center,'-'),2,fit$scale,'/'))
  B <- fit$dependence$psill*exp(-pair_distance(xy,fit$xy)/fit$dependence$range)
  as.numeric(F%*%fit$beta+B%*%fit$alpha)
}

make_application_geometry <- function(dataset,target_type,repeat_id) {
  d <- dataset$data
  dataset_id <- if(dataset$primary_target=='published') 1L else 2L
  target_id <- match(target_type,c('published','uniform','north'))
  split_seed <- 1200000L+10000L*dataset_id+1000L*target_id+repeat_id
  fold_seed <- 1300000L+10000L*dataset_id+1000L*target_id+repeat_id
  set.seed(split_seed)
  target_rows <- switch(target_type,
    published=which(d$published_role=='validation'),
    uniform=sort(sample(seq_len(nrow(d)),dataset$n_holdout)),
    north=tail(order(d$y_m,d$site_id),dataset$n_holdout))
  train_rows <- setdiff(seq_len(nrow(d)),target_rows)
  xy_m <- as.matrix(d[train_rows,c('x_m','y_m')])
  target_m <- as.matrix(d[target_rows,c('x_m','y_m')])
  origin <- colMeans(xy_m)
  xy <- sweep(xy_m,2,origin,'-')/1000
  target_xy <- sweep(target_m,2,origin,'-')/1000
  D <- pair_distance(xy); diag(D) <- Inf
  buffer <- median(apply(D,1,min))
  set.seed(fold_seed)
  folds <- make_candidates(xy,k=5L,buffer=buffer)
  points <- function(z) sf::st_as_sf(data.frame(x=z[,1],y=z[,2]),coords=c('x','y'),crs=dataset$crs)
  set.seed(fold_seed+100000L)
  kn_messages <- character()
  kn <- withCallingHandlers(suppressMessages(CAST::knndm(points(xy_m),
    predpoints=points(target_m),k=5L)),warning=function(w) {
      kn_messages <<- c(kn_messages,conditionMessage(w)); invokeRestart('muffleWarning') })
  folds$knndm <- Map(function(te,tr) list(test=te,train=tr),kn$indx_test,kn$indx_train)
  D <- pair_distance(xy); Dq <- pair_distance(target_xy,xy)
  dp <- apply(Dq,1,min)
  dc <- list()
  for(scheme in names(folds)) {
    fs <- folds[[scheme]]
    stopifnot(length(fs)==5L,
      identical(sort(as.integer(unlist(lapply(fs,`[[`,'test')))),seq_along(train_rows)))
    dc[[scheme]] <- rep(NA_real_,length(train_rows))
    for(f in fs) {
      stopifnot(length(f$train)>=8,length(f$test)>0,
                length(intersect(f$train,f$test))==0,
                all(c(f$train,f$test)%in%seq_along(train_rows)))
      dc[[scheme]][f$test] <- apply(D[f$test,f$train,drop=FALSE],1,min)
    }
    stopifnot(all(is.finite(dc[[scheme]])))
  }
  # CAST computes in metres; all diagnostics here use kilometres.
  W <- distance_summary(dc$knndm,dp,function(h) rep(0,length(h)))[['W']]
  stopifnot(abs(W*1000-kn$W)<max(1e-4,abs(kn$W)*1e-7),min(dp)>0)
  list(train_rows=train_rows,target_rows=target_rows,xy=xy,target_xy=target_xy,
       folds=folds,dc=dc,dp=dp,buffer_km=buffer,split_seed=split_seed,fold_seed=fold_seed,
       knndm_warnings=paste(unique(kn_messages),collapse=' | '))
}

analyse_application <- function(dataset_name,target_type,repeat_id,datasets) {
  ds <- datasets[[dataset_name]]
  geo <- make_application_geometry(ds,target_type,repeat_id)
  y <- ds$data$response[geo$train_rows]
  # Target outcomes are read for scoring only after folds/geometry are fixed.
  target_y <- ds$data$response[geo$target_rows]
  rows <- fit_rows <- prediction_rows <- list()
  full_models <- list()
  for(model in c('ordinary','trend')) {
    full_fit <- fit_spatial(geo$xy,y,model)
    full_models[[model]] <- full_fit
    full_record <- full_fit$record; full_record$fit <- 'full'
    fit_rows[[length(fit_rows)+1L]] <- full_record
    target_prediction <- predict_spatial(full_fit,geo$target_xy)
    target_mse <- mean((target_prediction-target_y)^2)
    # Same marginal response covariance for both prediction models and all
    # schemes; trend-residual covariance is a separately labelled sensitivity.
    marginal <- full_models$ordinary$dependence
    for(scheme in names(geo$folds)) {
      folds <- geo$folds[[scheme]]
      cv_prediction <- rep(NA_real_,length(y))
      for(fi in seq_along(folds)) {
        f <- folds[[fi]]
        fit <- fit_spatial(geo$xy[f$train,,drop=FALSE],y[f$train],model)
        cv_prediction[f$test] <- predict_spatial(fit,geo$xy[f$test,,drop=FALSE])
        rec <- fit$record; rec$fit <- paste(scheme,fi,sep=':')
        fit_rows[[length(fit_rows)+1L]] <- rec
      }
      stopifnot(all(is.finite(cv_prediction)),all(is.finite(target_prediction)))
      diagnostics <- distance_summary(geo$dc[[scheme]],geo$dp,marginal$rho)
      residual <- distance_summary(geo$dc[[scheme]],geo$dp,full_fit$dependence$rho)
      score <- mean((y-cv_prediction)^2)
      rows[[length(rows)+1L]] <- data.frame(dataset=dataset_name,target=target_type,
        repeat_id=repeat_id,model=model,scheme=scheme,n_train=length(y),n_target=length(target_y),
        cv_mse=score,target_mse=target_mse,gap=target_mse-score,absolute_error=abs(target_mse-score),
        A=diagnostics[['A']],W=diagnostics[['W']],J1=diagnostics[['J1']],J2=diagnostics[['J2']],
        U1=diagnostics[['U1']],U2=diagnostics[['U2']],crosses=diagnostics[['crosses']],
        J2_residual=residual[['J2']],variance_J2=(marginal$psill+marginal$nugget)*diagnostics[['J2']],
        buffer_km=geo$buffer_km,range_km=marginal$range,signal=marginal$signal_prop,
        residual_range_km=full_fit$dependence$range,residual_signal=full_fit$dependence$signal_prop,
        min_train=min(vapply(folds,function(f) length(f$train),integer(1))),
        mean_train=mean(vapply(folds,function(f) length(f$train),integer(1))),
        split_seed=geo$split_seed,fold_seed=geo$fold_seed)
      prediction_rows[[length(prediction_rows)+1L]] <- data.frame(dataset=dataset_name,target=target_type,
        repeat_id=repeat_id,model=model,scheme=scheme,role='cv',
        site_id=ds$data$site_id[geo$train_rows],observed=y,predicted=cv_prediction,
        squared_error=(y-cv_prediction)^2)
    }
    prediction_rows[[length(prediction_rows)+1L]] <- data.frame(dataset=dataset_name,target=target_type,
      repeat_id=repeat_id,model=model,scheme='full',role='target',
      site_id=ds$data$site_id[geo$target_rows],observed=target_y,predicted=target_prediction,
      squared_error=(target_y-target_prediction)^2)
  }
  fits <- do.call(rbind,fit_rows)
  fits$dataset <- dataset_name; fits$target <- target_type; fits$repeat_id <- repeat_id
  list(results=do.call(rbind,rows),fits=fits,predictions=do.call(rbind,prediction_rows),
       geometry=geo,full_models=full_models)
}
