# Exact copy of primary fit_spatial, changing only covariance estimation and adding audit fields.
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
  dep <- withCallingHandlers(profile_dependence(
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
  for(name in names(dep$profile_audit)) fit_record[[name]] <- dep$profile_audit[[name]]
  list(xy=xy,model=model,center=center,scale=scale,beta=beta,alpha=alpha,
       dependence=dep,record=fit_record)
}

