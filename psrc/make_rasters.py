

import netCDF4
import numpy as np
import sys

from raster_renderer import basemap_raster_mercator


def render_rasters(vs,ncpath):

  d = netCDF4.Dataset(ncpath)

  lat,lon = d.variables['Lat'][:,:],d.variables['Lon'][:,:]
  bounds = None

  for v in vs:

    vals = None
    if len(v)==4:
      vn,ndx,oname,longname = v
      vals = d.variables[vn][:,:,ndx]
    else:
      vn,oname,longname = v
      vals = d.variables[vn][:,:]

    raster,bounds = basemap_raster_mercator(lon,lat,vals)

    # dump binary data to file
    with open('outputs/' + oname + '.png','w') as f:
      f.write(raster)

  # all bounds will be exactly the same
  return bounds


if __name__ == '__main__':

  if len(sys.argv) < 2:
    print('usage: %s <state-file>')

  # generate the rasters
  ncpath = sys.argv[1]
  fields = [('FMC_GC',0,'fmc-1hr','1-hr'),
            ('FMC_GC',1,'fmc-10hr','10-hr'),
            ('FMC_GC',2,'fmc-100hr','100-hr'),
            ('Ed','drying-fme','Drying FME'),
            ('Ew','wetting-fme','Wetting FME'),
            ('RELH','rel-humidity','Relative Humidity'),
            ('T2','temp-2m','Temperature 2m')]
  Nf = len(fields)

  bounds = render_rasters(fields,ncpath)
  minlat,maxlat = bounds[0][1],bounds[2][1]
  minlon,maxlon = bounds[0][0],bounds[2][0]

  with open('data/rasters.json','w') as f:
    f.write('[\n')
    for i in range(Nf):
      fldinfo = fields[i]
      N = len(fldinfo)
      f.write('  {\n')
      f.write('    "name": "'+fldinfo[N-1]+'",\n')
      f.write('    "bounds": [ [ ' + str(minlat) + ', ' + str(minlon) + '], [ ' + str(maxlat) + ', ' + str(maxlon) + ' ] ],\n')
      f.write('    "url": "rasters/' + fldinfo[N-2] + '.png"\n')
      if i < Nf-1:
        f.write('  },\n')
      else:
        f.write('  }\n')
    f.write(']\n')
