#!/bin/sh

NTHREAD=16
NAIAD=/Users/williamdawn/work/naiad/src/naiad.x

for CASE in $( cat manifest.txt )
do
  cd ${CASE}
  OMP_NUM_THREAD=${NTHREAD} ${NAIAD} ${CASE}.inp
  cd ..
done
