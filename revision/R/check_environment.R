source("revision/R/local_libraries.R")
cat("R:", R.version.string, "\n")
cat("Host:", Sys.info()[["nodename"]], "\n")
for (p in c("CAST", "sf", "ranger", "gstat", "pkgload", "jsonlite", "ggplot2", "patchwork", "digest")) {
  loaded <- tryCatch(requireNamespace(p, quietly = TRUE), error = function(e) FALSE)
  cat(p, if (loaded) as.character(utils::packageVersion(p)) else "UNAVAILABLE", "\n")
}
print(.libPaths())
print(sessionInfo())
