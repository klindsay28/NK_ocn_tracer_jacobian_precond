#!/bin/csh -f

rm -f Depends.new
touch Depends.new

echo 'OBJS =' > tempfile
foreach file ( `grep -c '^main' *.c | grep -v :1 | cut -f1 -d:` )
   echo ${file:r}.o >> tempfile
end
cat tempfile | paste -s -d' ' >> Depends.new

foreach file ( *.c )
   echo ${file:r}.o: $file >! tempfile
   grep 'include "' $file | cut -f2 -d\" | grep -v -E 'netcdf|superlu_ddefs|mpi' >> tempfile
   cat tempfile | paste -s -d' ' >> Depends.new
   rm -f tempfile
end

diff Depends Depends.new

mv Depends.new Depends

