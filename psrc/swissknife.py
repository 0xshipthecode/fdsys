

import netCDF4
import sys

from spatial_model_utilities import find_closest_grid_point,great_circle_distance


def get_grid_values_from_file(ncpath,varname,lat,lon):
  d = netCDF4.Dataset(ncpath)
  glat, glon = d.variables['Lat'][:,:], d.variables['Lon'][:,:]
  V = d.variables[varname]

  i, j = find_closest_grid_point(lon, lat, glon, glat)
  dist = great_circle_distance(glon[i,j],glat[i,j],lon,lat)
  vals = V[i,j,...]
  d.close()

  return dist, i, j, vals


def main(args):

  if len(args) < 3:
    print("usage: %s <cmd> [args] <nc-file>\n")
    print("current commands:")
    print("display var lat lon <nc-file> ->  shows all values of var nearest to lat/lon + distance to point")
    sys.exit(1)


  if args[1] == 'display':
    varname, lat, lon = args[2], float(args[3]), float(args[4])
    ncpath = args[5]
    dist,i,j,vals = get_grid_values_from_file(ncpath,varname,lat,lon)
    print("Distance from [%g,%g] to nearest grid point [%d,%d] is %g km." % (lat,lon,i,j,dist))
    print("Values of [%s]: " % varname)
    print(vals)



if __name__ == "__main__":
  main(sys.argv)
