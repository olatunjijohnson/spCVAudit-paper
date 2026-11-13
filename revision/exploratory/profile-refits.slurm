#!/bin/bash --login
#SBATCH --job-name=spcv-profile
#SBATCH --partition=multicore
#SBATCH --ntasks=8
#SBATCH --time=08:00:00
#SBATCH --output=/net/scratch/x27905oj/spLeakage-profile-refits-20261113/profile-refits-%j.log
# Profile-search refits of the estimated-covariance studies
# (revision/audit/profile-refit-protocol.json), in the reference environment
# used for the retained fitted studies.
set -euo pipefail
module purge
module load apps/gcc/R/4.5.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TZ=Etc/UTC
export R_LIBS_USER=/net/scratch/x27905oj/spLeakage-revision-20261028/revision/csf/library:/mnt/iusers01/maths01/x27905oj/R/x86_64-pc-linux-gnu-library/4.5:/opt/apps/apps/gcc/R/4.5.0/lib64/R/library:/mnt/iusers01/maths01/x27905oj/R/x86_64-pc-linux-gnu-library/4.4:/mnt/iusers01/maths01/x27905oj/R/x86_64-pc-linux-gnu-library/4.3:/mnt/iusers01/maths01/x27905oj/R/x86_64-pc-linux-gnu-library/4.1:/mnt/iusers01/maths01/x27905oj/R/x86_64-pc-linux-gnu-library/3.6
cd /net/scratch/x27905oj/spLeakage-profile-refits-20261113/source
git rev-parse HEAD
git status --short
# Same numerical environment as the retained runs, or stop before any study.
Rscript --vanilla revision/R/check_reference_fitter.R .
for study in sensitivity nested nested-k10; do
  echo "== $study $(date -u +%FT%TZ)"
  Rscript --vanilla revision/exploratory/run_profile_refits.R "$study"
done
echo "== compare $(date -u +%FT%TZ)"
Rscript --vanilla revision/exploratory/compare_profile_refits.R
echo "== done $(date -u +%FT%TZ)"
