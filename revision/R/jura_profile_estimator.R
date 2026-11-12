# For each range the variogram is linear in nugget and partial sill. Enumerate
# the constrained quadratic optimum and each active boundary analytically.
profile_dependence <- function(data,response,coords) {
  original <- spLeakage::estimate_dependence(data,response,coords=coords)
  vg <- original$variogram
  h_scale <- max(vg$dist); g_scale <- max(vg$gamma)
  if(g_scale<=0) return(original)
  h <- vg$dist/h_scale; y <- vg$gamma/g_scale; w <- vg$np/sum(vg$np)
  ey <- sum(w*y)
  linear_solution <- function(log_range,coefficients=FALSE) {
    x <- -expm1(-h/exp(log_range))
    ex <- sum(w*x); exx <- sum(w*x*x); exy <- sum(w*x*y)
    clip <- function(z,lo,hi) pmin(pmax(z,lo),hi)
    candidates <- rbind(c(0,clip(exy/exx,1e-8,2)),
                        c(1,clip((exy-ex)/exx,1e-8,2)),
                        c(clip(ey-1e-8*ex,0,1),1e-8),
                        c(clip(ey-2*ex,0,1),2))
    if(exx-ex*ex>1e-15) {
      b <- (exy-ex*ey)/(exx-ex*ex); a <- ey-b*ex
      if(a>=0&&a<=1&&b>=1e-8&&b<=2) candidates <- rbind(candidates,c(a,b))
    }
    loss <- apply(candidates,1,function(p) sum(w*(p[1]+p[2]*x-y)^2))
    best <- which.min(loss)
    if(coefficients) return(c(nugget=candidates[best,1],psill=candidates[best,2],
                             range=exp(log_range),loss=loss[best]))
    loss[best]
  }
  grid <- seq(log(1e-6),log(5),length.out=81L)
  loss <- vapply(grid,linear_solution,numeric(1))
  local <- which(loss[2:80]<=loss[1:79]&loss[2:80]<=loss[3:81])+1L
  # Flat plateaus at unidentifiably short scales need only their best grid value;
  # refining every equal-valued point adds work without additional information.
  local <- local[vapply(local,function(i) loss[i]<max(loss[i-1],loss[i+1])-1e-14,logical(1))]
  trials <- c(grid[which.min(loss)],grid[1],tail(grid,1))
  for(i in local) trials <- c(trials,optimize(linear_solution,c(grid[i-1],grid[i+1]),tol=1e-10)$minimum)
  solutions <- lapply(trials,linear_solution,coefficients=TRUE)
  best <- solutions[[which.min(vapply(solutions,function(p) p[['loss']],numeric(1)))]]
  old_loss <- sum(w*((original$nugget+original$psill*(1-exp(-vg$dist/original$range)))/g_scale-y)^2)
  changed <- best[['loss']]<old_loss-1e-12
  out <- original
  if(changed) {
    out$nugget <- best[['nugget']]*g_scale
    out$psill <- best[['psill']]*g_scale
    out$range <- best[['range']]*h_scale
    out$practical_range <- 3*out$range
    out$signal_prop <- out$psill/(out$psill+out$nugget)
    out$coef <- list(nugget=out$nugget,psill=out$psill,range=out$range)
    out$vcov <- NULL
    signal <- out$signal_prop; range <- out$range
    out$rho <- function(h) signal*exp(-h/range)
  }
  new_loss <- sum(w*((out$nugget+out$psill*(1-exp(-vg$dist/out$range)))/g_scale-y)^2)
  stopifnot(new_loss<=old_loss+1e-12)
  out$profile_audit <- c(original_objective=old_loss,profile_objective=new_loss,
    improved=as.numeric(changed),original_range=original$range,original_signal=original$signal_prop)
  out
}
