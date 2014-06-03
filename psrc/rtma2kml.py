#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import simplekml as kml
from mpl_toolkits.basemap import Basemap
import numpy as np
import netCDF4 as nc4
import sys
import os
import StringIO

from raster_renderer import make_colorbar, basemap_raster_mercator


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print('Usage: %s <rtma-file> <varname> <units> <target-file>' % sys.argv[0])
        sys.exit(1)

    # open the netCDF dataset
    d = nc4.Dataset(sys.argv[1])
    varname = sys.argv[2]
    cb_units = sys.argv[3]
    outf = sys.argv[4]

    # extract variables
    fa = d.variables[varname][:,:,0]
    lon = d.variables['Lon'][:,:]
    lat = d.variables['Lat'][:,:]

    # construct kml file
    doc = kml.Kml(name = varname)
    fa_min,fa_max = np.nanmin(fa),np.nanmax(fa)

    print("[rtma2kml] rendering colorbar and colorbar as a screen overlay ...")
    cb_png_data = make_colorbar([fa_min,fa_max],'vertical',2,mpl.cm.jet,cb_units,varname)
    with open('colorbar.png', 'w') as f:
        f.write(cb_png_data)
    doc.addfile('colorbar.png')

    cbo = doc.newscreenoverlay(name='colorbar')
    cbo.overlayxy = kml.OverlayXY(x=0,y=1,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
    cbo.screenxy = kml.ScreenXY(x=0.02,y=0.98,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
    cbo.size = kml.Size(x=0,y=0,xunits=kml.Units.fraction,yunits=kml.Units.fraction)
    cbo.color = kml.Color.rgb(255,255,255,a=150)
    cbo.visibility = 1
    cbo.icon.href="colorbar.png"

    print("[rtma2kml] raster range [%g, %g], mean %g" % (np.amin(fa),np.amax(fa),np.mean(fa)))
    print("[rtma2kml] rendering raster from variable %s (Mercator projection) ..." % varname)
    ground = doc.newgroundoverlay(name=varname,color="80ffffff")
    raster_png_data,corner_coords = basemap_raster_mercator(lon,lat,fa)
    with open('raster.png','w') as f:
      f.write(raster_png_data)
    doc.addfile('raster.png')
    ground.icon.href="raster.png"
    ground.gxlatlonquad.coords = corner_coords

    doc.savekmz(outf)

    # cleanup
    print("[rtma2kml] cleaning up temp files ...")
    os.remove('colorbar.png')
    os.remove('raster.png')

    print("[rtma2kml] done.")