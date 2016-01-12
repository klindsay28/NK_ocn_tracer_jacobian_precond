#!/bin/csh -f

foreach file ( *.h *.c )
   echo $file
   indent $file
   diff ${file}~ $file
   if ! $status mv ${file}~ $file
end
