#!/bin/bash
set -o nounset
set -o errexit

if [ $# -ne 1 ]; then 
  echo "usage: grib_to_netcdf_rt.sh <directory>"
  exit 1
fi

FILES=$(ls $1/*.bin)

for F in $FILES ; do
  B=$(basename $F .bin)
  echo $B
  grib2/wgrib2/wgrib2 $F -netcdf $1/$B.nc
done

