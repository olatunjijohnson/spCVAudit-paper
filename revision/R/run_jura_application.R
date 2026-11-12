source('revision/R/jura_common.R')
original_fit <- '--original-fit' %in% commandArgs(trailingOnly=TRUE)
if (!original_fit) {
  source('revision/R/jura_profile_estimator.R')
  source('revision/R/jura_profile_fit.R')
}
if(!file.exists('revision/data/jura-data.rds')) source('revision/R/prepare_jura_data.R')
datasets <- list(jura=readRDS('revision/data/jura-data.rds'))
grid <- expand.grid(target=c('published','uniform','north'),repeat_id=seq_len(40L),stringsAsFactors=FALSE)
grid <- grid[order(grid$repeat_id,grid$target),]; rownames(grid)<-NULL
grid$id <- seq_len(nrow(grid))
grid$seed <- 1210000L+1000L*match(grid$target,c('published','uniform','north'))+grid$repeat_id
stem <- if (original_fit) 'application-jura-original' else 'application-jura'
workers <- as.integer(Sys.getenv('SLURM_NTASKS','4'))
sources <- c('revision/R/run_jura_application.R','revision/R/jura_common.R',
             'revision/R/jura_profile_estimator.R','revision/R/jura_profile_fit.R',
             'revision/R/common.R','revision/data/jura-data.rds',list.files('package/R',full.names=TRUE,pattern='[.]R$'))
hashes <- tools::md5sum(sources)
outdir <- file.path('revision/results',paste0(stem,'-batches'))
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
write.csv(grid,file.path('revision/results',paste0(stem,'-configuration.csv')),row.names=FALSE)
one <- function(i) {
  cfg <- grid[i,]; path <- file.path(outdir,sprintf('case-%03d.rds',cfg$id))
  if(file.exists(path)) {
    old <- readRDS(path); stopifnot(identical(old$source_md5,hashes)); return(old)
  }
  ans <- tryCatch(analyse_application('jura',cfg$target,cfg$repeat_id,datasets),
    error=function(e) list(error=conditionMessage(e),configuration=cfg))
  if(is.null(ans$error)) {
    for(nm in c('results','fits','predictions')) ans[[nm]]$id <- cfg$id
  }
  ans$source_md5 <- hashes
  saveRDS(ans,path); ans
}
records <- list(); started <- proc.time()[[3]]
for(batch in split(seq_len(nrow(grid)),ceiling(seq_len(nrow(grid))/8))) {
  out <- parallel::mclapply(batch,one,mc.cores=workers,mc.set.seed=FALSE)
  records <- c(records,out)
  cat('Completed',max(batch),'/',nrow(grid),'in',round(proc.time()[[3]]-started,1),'seconds\n');flush.console()
}
errors <- Filter(function(x) !is.null(x$error),records)
saveRDS(errors,file.path('revision/results',paste0(stem,'-errors.rds')))
if(length(errors)) stop('Jura configurations failed; inspect retained errors')
combine <- function(name) do.call(rbind,lapply(records,`[[`,name))
bundle <- list(results=combine('results'),fits=combine('fits'),predictions=combine('predictions'),details=records,
  metadata=list(time_utc=format(Sys.time(),tz='UTC',usetz=TRUE),config=grid,session=sessionInfo(),
    package_version=as.character(packageVersion('spLeakage')),executed_source_md5=hashes,
    fitter_body=paste(deparse(body(if (original_fit) spLeakage::estimate_dependence else profile_dependence)),collapse='\n'),
    scope='Repeated analysis of one Jura survey; observed holdout MSE, not expected repeated-field risk.'))
stopifnot(nrow(bundle$results)==960L,nrow(bundle$fits)==5040L)
saveRDS(bundle,file.path('revision/results',paste0(stem,'.rds')))
write.csv(bundle$results,file.path('revision/results',paste0(stem,'.csv')),row.names=FALSE)
cat('Jura complete:',nrow(bundle$results),'evaluations and',nrow(bundle$fits),'covariance fits.\n')
