#!/bin/bash
set -
FILES=$(ls inputs/*.bin)

for F in $FILES ; do
  B=$(basename $F .bin)
  echo $B
  grib2/wgrib2/wgrib2 $F -netcdf inputs/$B.nc
done

