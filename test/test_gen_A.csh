#!/bin/csh -f

echo entering test_gen_A.csh
date

source /glade/apps/opt/lmod/lmod/init/csh

source ../scripts/shared_vars

setenv out_fname /glade/scratch/klindsay/matrix.nc
rm -f $out_fname

set circ_fname = $WDIR/c.e12.C.T62_g37.pre.001.pop.h.0001.nc
set circ_fname = $WDIR/c.e12.C.T62_g16.no_ovf.001.pop.h.0151.nc
set circ_fname = $WDIR/c.e12.C.T62_g16.ovf.001.pop.h.0151.nc

setenv day_cnt    365.0    # number of days in forcing

module load job_memusage

job_memusage.exe ./bin/gen_A -D1 \
  -o adv_type,upwind3 -o hmix_type,isop_file \
  -o vmix_type,file -o sink_type,const_shallow,365.0,10.0e2 \
  -c $day_cnt -o reg_fname,$regfile $circ_fname $out_fname
set ret_code = $status
if ($ret_code) then
  echo gen_A returned with error code $ret_code
  echo exiting test_gen_A.csh
  date
  exit 1
endif
date

ncdump -h $out_fname

echo returning from test_gen_A.csh
date

