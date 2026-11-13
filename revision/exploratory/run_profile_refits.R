# Post-hoc fitting sensitivity; declared in revision/audit/profile-refit-protocol.json.
# Runs an unmodified study runner with every variogram fit replaced by the Jura
# profile search. Usage: Rscript revision/exploratory/run_profile_refits.R <study> [--pilot]
study <- commandArgs(trailingOnly = TRUE)[1]
runners <- list(
  sensitivity = list(file = "revision/R/run_sensitivity.R", fits = 1L,
    stem = c(), suffix = "-rng-profile"),
  nested = list(file = "revision/R/run_nested_estimation.R", fits = 1L,
    stem = c(), suffix = "-rng-profile"),
  `nested-k10` = list(file = "revision/exploratory/run_nested_k10.R", fits = 1L,
    stem = c(`"nested-estimation-k10"` = '"nested-estimation-k10-profile"',
             `"nested-k10-pilot"` = '"nested-k10-profile-pilot"'), suffix = ""))
if (is.na(study) || !study %in% names(runners)) stop("Study must be one of: ", paste(names(runners), collapse = ", "))
spec <- runners[[study]]
source("revision/R/local_libraries.R")
source("revision/R/jura_profile_estimator.R")
Sys.setenv(REVISION_RUN_SUFFIX = spec$suffix)

code <- readLines(spec$file)
call <- "spLeakage::estimate_dependence("
replaced <- sum(lengths(regmatches(code, gregexpr(call, code, fixed = TRUE))))
if (replaced != spec$fits) stop("Expected ", spec$fits, " fitter call(s) in ", spec$file, ", found ", replaced)
code <- gsub(call, "profile_dependence(", code, fixed = TRUE)
for (old in names(spec$stem)) {
  if (sum(grepl(old, code, fixed = TRUE)) != 1L) stop("Expected one stem ", old, " in ", spec$file)
  code <- gsub(old, spec$stem[[old]], code, fixed = TRUE)
}

eval(parse(text = code, keep.source = FALSE), envir = globalenv())

# The runner has saved its bundle under `stem`; record the substitution.
path <- file.path("revision/results", paste0(stem, ".rds"))
bundle <- readRDS(path)
bundle$metadata$profile_refit <- list(
  protocol = "revision/audit/profile-refit-protocol.json",
  runner = spec$file, runner_md5 = unname(tools::md5sum(spec$file)),
  wrapper_md5 = unname(tools::md5sum("revision/exploratory/run_profile_refits.R")),
  estimator_md5 = unname(tools::md5sum("revision/R/jura_profile_estimator.R")),
  protocol_md5 = unname(tools::md5sum("revision/audit/profile-refit-protocol.json")),
  replaced_calls = replaced,
  profile_fitter_body = paste(deparse(body(profile_dependence)), collapse = "\n"),
  note = "fitter_body/dependence_implementation above record the package fitter used as the profile starting point.")
saveRDS(bundle, path)
cat("Profile refit", study, "saved to", path, "\n")
