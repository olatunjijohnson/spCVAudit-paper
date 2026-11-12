# Export the documented 259/100 Jura split without any target-dependent choices.
source('revision/R/local_libraries.R')
e <- new.env(); data('jura',package='gstat',envir=e)
j <- rbind(e$jura.pred,e$jura.val)
j$site_id <- c(sprintf('train-%03d',seq_len(nrow(e$jura.pred))),
               sprintf('validation-%03d',seq_len(nrow(e$jura.val))))
j$published_role <- rep(c('training','validation'),c(nrow(e$jura.pred),nrow(e$jura.val)))
points <- sf::st_transform(sf::st_as_sf(j,coords=c('long','lat'),crs=4326,remove=FALSE),2056)
coords <- sf::st_coordinates(points)
j$x_m <- coords[,1]; j$y_m <- coords[,2]; j$response <- log(j$Cd)
stopifnot(nrow(j)==359L,all(is.finite(j$response)),all(j$Cd>0),
          !anyDuplicated(j[c('x_m','y_m')]))
jura <- list(data=j,crs=2056,units='squared natural-log cadmium',response='Natural log cadmium (mg/kg)',
             source='gstat::jura; georeferenced jura.pred and jura.val; EPSG:2056',
             primary_target='published',n_holdout=100L)
saveRDS(jura,'revision/data/jura-data.rds')
write.csv(j,'revision/data/jura-sites.csv',row.names=FALSE)
writeLines(c('Jura topsoil data: cadmium concentration and coordinates.',
  'Source: gstat::jura, jura.pred (259 sites) and jura.val (100 sites).',
  paste('Exported with gstat',packageVersion('gstat')),
  'Documentation: https://r-spatial.github.io/gstat/reference/jura.html',
  'Atteia, Dubois and Webster (1994), Environmental Pollution 86:315-327.',
  'Goovaerts (1997), Geostatistics for Natural Resources Evaluation, Appendix C.',
  'Geographic coordinates transformed from WGS84 to EPSG:2056.',
  'response is natural-log Cd; Cd is recorded in mg/kg of topsoil.',
  'Original training/validation roles and all original data columns are retained.'),
  'revision/data/JURA-SOURCE.txt')
