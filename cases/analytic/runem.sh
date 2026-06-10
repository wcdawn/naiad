#!/bin/sh

NTHREAD=16
NAIAD="$(pwd)/../../src/naiad.x"

for CASE in $( cat manifest.txt )
do
  cd ${CASE}
  OMP_NUM_THREAD=${NTHREAD} ${NAIAD} ${CASE}.inp || exit 1
  cd ..
done
