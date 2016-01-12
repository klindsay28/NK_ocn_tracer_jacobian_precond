#!/bin/csh -f

#BSUB -P P93300070                  # project number
#BSUB -W 0:45                       # wall clock time (in minutes)
#BSUB -a poe                        # select the poe elim
#BSUB -n 144                        # number of tasks
#BSUB -R "span[ptile=2]"            # tasks per node
#BSUB -J jacobian_precond           # job name
#BSUB -oo jacobian_precond.%J.out   # output filename
#BSUB -q small                      # queue

# setenv TARGET_CPU_LIST "-1"
setenv MP_LABELIO yes
setenv MP_STDOUTMODE unordered

echo entering test_solve_ABdist.csh
date

source ../scripts/shared_vars

setenv in_fname     $WDIR/fcn_eval_000.nc
setenv out_fname    /glade/scratch/klindsay/B_dist.nc
setenv matrix_fname /glade/scratch/klindsay/matrix.nc

rm -f $out_fname
cp $in_fname $out_fname

source /etc/profile.d/modules.csh
module load job_memusage

mpirun.lsf job_memusage.exe ./bin/solve_ABdist -D1 -n12,12 \
  -v IAGE_RESTORE_1DAY_CUR $matrix_fname $out_fname
set ret_code = $status
if ($ret_code) then
  echo solve_ABdist returned with error code $ret_code
  echo exiting test_solve_ABdist.csh
  date
  exit 1
endif
date

ncdiff -A -v $varlist_perturb_cur $out_fname $in_fname $out_fname

echo returning from test_solve_ABdist.csh
date

