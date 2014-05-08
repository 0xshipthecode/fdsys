#!/bin/bash
set -o nounset
set -o errexit

if [ $# -ne 1 ]; then 
  echo "usage: grib_to_netcdf_archive.sh <directory>"
  exit 1
fi

FILES=$(ls $1/L*)

for F in $FILES ; do
  B=$(basename $F)
  echo $B
  grib2/wgrib2/wgrib2 -not error $F -netcdf $1/$B.nc
done

