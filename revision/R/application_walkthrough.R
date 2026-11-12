# Public-API walkthrough for the Jura published target, repetition 1.
source("revision/R/local_libraries.R")
pkgload::load_all("package", quiet=TRUE)
bundle <- readRDS("revision/results/application-jura.rds")
detail <- Filter(function(d) d$results$target[1]=="published" && d$results$repeat_id[1]==1L,bundle$details)[[1]]
geo <- detail$geometry
target <- spLeakage::prediction_target(grid=geo$target_xy)
audits <- lapply(geo$folds,function(folds) spLeakage::detect_leakage(geo$xy,folds,target,
  dependence=detail$full_models$ordinary$dependence))
table <- do.call(rbind,lapply(names(audits),function(scheme) {
  x<-audits[[scheme]]
  data.frame(scheme,A_km=x$A,W_km=x$W,J1=x$SLI_rho,J2=x$SLI_rho2,U1=x$W_rho,crosses=x$crosses)
}))
reference <- subset(bundle$results,target=="published" & repeat_id==1L & model=="ordinary")
reference <- reference[match(table$scheme,reference$scheme),]
stopifnot(max(abs(as.matrix(table[c("A_km","W_km","J1","J2","U1")])-
  as.matrix(reference[c("A","W","J1","J2","U1")])))<1e-7)
write.csv(table,"revision/results/application-public-workflow.csv",row.names=FALSE)
pdf("revision/figures/application-public-workflow.pdf",width=8,height=6)
par(mfrow=c(2,2),mar=c(4,4,3,1))
for(scheme in names(audits)) {plot(audits[[scheme]]);title(main=scheme)}
invisible(dev.off())
print(table,row.names=FALSE)
cat("Distances use kilometres. These diagnostics describe nearest-neighbour distributions;\n",
    "observed target MSE is evaluated separately using the 100 reserved soil measurements.\n",sep="")
