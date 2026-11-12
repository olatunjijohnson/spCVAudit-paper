source('revision/R/local_libraries.R')
b <- readRDS('revision/results/application-jura.rds'); z <- b$results
stopifnot(nrow(z)==960L,all(table(z$id)==8L),all(is.finite(z$cv_mse)),all(is.finite(z$target_mse)))
outdir <- 'revision/results/summaries';dir.create(outdir,showWarnings=FALSE)
grouped <- split(seq_len(nrow(z)),interaction(z[c('target','model','scheme')],drop=TRUE,lex.order=TRUE))
candidate <- do.call(rbind,lapply(grouped,function(ii) {
  d <- z[ii,];r <- d[1,c('target','model','scheme','n_train','n_target')]; r$n_repeats<-nrow(d)
  for(v in c('cv_mse','target_mse','gap','absolute_error','W','J2','buffer_km','range_km','signal')) {
    r[[paste0(v,'_mean')]] <- mean(d[[v]])
    for(q in c(25,50,75)) r[[paste0(v,'_q',q)]] <- unname(quantile(d[[v]],q/100))
  }
  r
}))
write.csv(candidate,file.path(outdir,'application-jura-candidates.csv'),row.names=FALSE)
selection <- list()
for(ii in split(seq_len(nrow(z)),interaction(z$id,z$model,drop=TRUE))) {
  d <- z[ii,]
  for(rule in c('W','J2','A','J1','U1','U2','J2_residual','fixed_random','fixed_knndm','best_observed')) {
    criterion <- switch(rule,fixed_random=as.numeric(d$scheme!='random'),
      fixed_knndm=as.numeric(d$scheme!='knndm'),best_observed=d$absolute_error,abs(d[[rule]]))
    take <- which(criterion<=min(criterion)+1e-12)
    selection[[length(selection)+1L]] <- data.frame(d[1,c('id','target','repeat_id','model')],rule,
      absolute_error=mean(d$absolute_error[take]),cv_mse=mean(d$cv_mse[take]),target_mse=d$target_mse[1],
      selected=paste(sort(d$scheme[take]),collapse=';'),n_tied=length(take))
  }
}
selection <- do.call(rbind,selection)
write.csv(selection,file.path(outdir,'application-jura-selection-repetitions.csv'),row.names=FALSE)
groups <- split(seq_len(nrow(selection)),interaction(selection[c('target','model','rule')],drop=TRUE))
summary <- do.call(rbind,lapply(groups,function(ii) {
  d <- selection[ii,];r <- d[1,c('target','model','rule')]
  r$mae <- mean(d$absolute_error);r$mean_cv_mse<-mean(d$cv_mse);r$mean_target_mse<-mean(d$target_mse)
  r$q25_absolute_error<-unname(quantile(d$absolute_error,.25));r$q75_absolute_error<-unname(quantile(d$absolute_error,.75))
  r$n_repeats<-nrow(d);r
}))
write.csv(summary,file.path(outdir,'application-jura-selection.csv'),row.names=FALSE)
details <- Filter(function(a) a$results$repeat_id[1]==1L,b$details)
geometry <- lapply(details,function(a) list(target=a$results$target[1],xy=a$geometry$xy,target_xy=a$geometry$target_xy,
  cv_distances=a$geometry$dc,target_distances=a$geometry$dp,buffer_km=a$geometry$buffer_km))
jsonlite::write_json(geometry,file.path(outdir,'application-jura-geometry.json'),auto_unbox=TRUE,digits=16)
lines <- c('% Generated from application-jura-selection.csv.',
  '\\begin{tabular}{llrrrrr}','\\toprule',
  'Target & Mean model & $W$ & $|J_2|$ & Random & kNNDM & Best of four \\\\','\\midrule')
for(target in c('published','uniform','north')) for(model in c('ordinary','trend')) {
  values <- vapply(c('W','J2','fixed_random','fixed_knndm','best_observed'),function(rule)
    summary$mae[summary$target==target&summary$model==model&summary$rule==rule],numeric(1))
  lines <- c(lines,sprintf('%s & %s & %.4f & %.4f & %.4f & %.4f & %.4f \\\\',
    c(published='Original split',uniform='Random',north='Northern')[[target]],
    c(ordinary='Constant',trend='Linear trend')[[model]],values[1],values[2],values[3],values[4],values[5]))
}
writeLines(c(lines,'\\bottomrule','\\end{tabular}'),'revision/manuscript/table-jura.tex')
cat('Generated Jura candidates, selection comparisons, geometry and supplementary table.\n')
