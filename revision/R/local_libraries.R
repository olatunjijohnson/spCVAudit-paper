# Use installed libraries from earlier R minor versions only as a fallback.
# The analysis never downloads or modifies these packages. Versions are recorded
# in each result's sessionInfo; smoke tests check the actual functions we use.
extra_libraries <- list.dirs(dirname(.libPaths()[1]), recursive = FALSE, full.names = TRUE)
extra_libraries <- extra_libraries[grepl("^[0-9]+[.][0-9]+$", basename(extra_libraries))]
.libPaths(unique(c(.libPaths(), rev(sort(extra_libraries)))))
# A dedicated CSF library takes precedence over user-wide installations.
project_library <- file.path("revision", "csf", "library")
if (dir.exists(project_library)) .libPaths(c(normalizePath(project_library), .libPaths()))
