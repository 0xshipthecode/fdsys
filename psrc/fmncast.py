# -*- coding: utf-8 -*-
"""
Created on Sun May 3rd 2014 - Martin Vejmelka

Performs a nowcasting step using current weather conditions, observations and
background guess.  Results in an analysis and new covariance.

@author: martin
"""

from trend_surface_model import fit_tsm
from grid_moisture_model import GridMoistureModel
from spatial_model_utilities import great_circle_distance, find_closest_grid_point
from observation import Observation

import numpy as np
import os
import sys
import string
from datetime import datetime
import pytz
import netCDF4


def total_seconds(tdelta):
    """
    Utility function for python < 2.7, 2.7 and above have total_seconds()
    as a function of timedelta.
    """
    return tdelta.microseconds / 1e6 + (tdelta.seconds + tdelta.days * 24 * 3600)


def load_raws_observations(obs_file,glat,glon):
  """
  Loads all of the RAWS observations valid at the time in question
  and converts them to Observation objects.
  """
  # load observations & register them to grid
  orig_obs = np.loadtxt(obs_file,delimiter=',')
  obss = []
  for oo in orig_obs:
      i, j = find_closest_grid_point(oo[6],oo[7],glat,glon)
      obst = datetime(int(oo[0]),int(oo[1]),int(oo[2]),int(oo[3]),int(oo[4]),int(oo[5]))
      obss.append(Observation(obst,oo[6],oo[7],oo[8],oo[9],(i,j)))
  return obss


def load_rtma_data(in_dir):
  """
  This function expects the files ds.temp.nc, ds.td.nc, ds.precipa.nc
  to be in the directory inputs and the file rtma_elevation.nc to be in
  the directory static.  All of the data is loaded and stored in a dictionary.
  """
  file_list = [("inputs/ds.temp.nc", "TMP_2maboveground", "T2"),
               ("inputs/ds.td.nc", "DPT_2maboveground", "DPT"),
               ("inputs/ds.precipa.nc", "APCP_surface", "RAIN")]

  rtma_time = None
  data = {}
  for rtmaf, rtmav, name in file_list:
      d = netCDF4.Dataset(rtmaf)
      data[name] = d.variables[rtmav][0,:,:]
      rtma_time = d.variables['time']
      d.close()

  # converted accumulated precip into rain intensity
  data['RAIN'] = 0.01 * data['RAIN']

  d = netCDF4.Dataset("static/rtma_elevation.nc")
  data['HGT'] = d.variables['HGT_surface'][0,:,:]
  Lat = d.variables['latitude'][:,:]
  Lon = d.variables['longitude'][:,:] - 360
  data['Lat'] = Lat
  data['Lon'] = Lon

  # compute the relative humidity according to NOAA formula [http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf]
  T2c = data['T2'] - 273.15
  DTc = data['DPT'] - 273.15
  e  = 6.11 * 10**(7.5 * DTc / (DTc + 237.3))
  es = 6.11 * 10**(7.5 * T2c / (T2c + 237.3))
  relh = e / es * 100
  data['RH'] = relh

  # chop all datasets to Colorado coordinates [36.9 <-> 41.1, -109.1 <-> -101.9]
  i1, i2, j1, j2 = 0, 0, Lat.shape[0], Lat.shape[1]
  done = False
  while not done:
    done = True
    tmp = np.where(np.amax(Lat[:, j1:j2],axis=1) < 37)[0][-1]
    if i1 != tmp:
      i1 = tmp
      done = False
    tmp = np.where(np.amin(Lat[:, j1:j2],axis=1) > 41)[0][0]
    if i2 != tmp:
      i2 = tmp
      done = False
    tmp = np.where(np.amax(Lon[i1:i2,:],axis=0) < -109)[0][-1]
    if j1 != tmp:
      j1 = tmp
      done = False
    tmp = np.where(np.amin(Lon[i1:i2,:],axis=0) > -102)[0][0]
    if j2 != tmp:
      j2 = tmp
      done = False

  for k,v in data.iteritems():
    data[k] = v[i1:i2,j1:j2]

  print("Chopping RTMA surface data to lat ndx [%d <-> %d], lon ndx [%d <-> %d]" % (i1,i2,j1,j2))

  # add time after slicing
  data['Time'] = datetime.fromtimestamp(rtma_time[0])

  return data


def compute_equilibria(T,H):
    d = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    w = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
    d *= 0.01
    w *= 0.01
    return d, w


def parse_datetime(s):
    gmt_tz = pytz.timezone('GMT')
    dt = datetime.strptime(s, "%Y/%m/%d  %H:%M:%S")
    return dt.replace(tzinfo=gmt_tz)


def run_module():

    # load RTMA data
    print("Loading RTMA data ...")
    data = load_rtma_data("inputs")

    # retrieve variables from RTMA
    lat, lon, hgt = data['Lat'], data['Lon'], data['HGT']
    t2, relh, rain = data['T2'], data['RH'], data['RAIN']
    ed, ew = compute_equilibria(t2,relh)
    tm = data['Time']

    dom_shape = lat.shape
    print('INFO: domain size is %d x %d grid points.' % dom_shape)
    print('INFO: domain extent is lats (%g to %g) lons (%g to %g).' % (np.amin(lat),np.amax(lat),np.amin(lon),np.amax(lon)))

    # initialize output file
    out_file = netCDF4.Dataset('outputs/fm.nc', 'w')
    out_file.createDimension('fuel_moisture_classes_stag', 5)
    out_file.createDimension('south_north', dom_shape[0])
    out_file.createDimension('west_east', dom_shape[1])
    ncfmc = out_file.createVariable('FMC_GC', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
    ncfmc_cov = out_file.createVariable('FMC_COV', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag', 'fuel_moisture_classes_stag'))
    print('INFO: opened outputs/fm.nc for writing assimilated data.')

    ### Load observation data from the stations

    # compute the diagonal distance between grid points
    grid_dist_km = great_circle_distance(lon[0,0], lat[0,0], lon[1,1], lat[1,1])
    print('INFO: diagonal distance in grid is %g' % grid_dist_km)

    obss = load_raws_observations('inputs/raws_ingest_201405041900.csv',lat,lon)
    print('Loaded %d observations.' % (len(obss)))
    print obss

    # set up parameters
    Nk = 3  # we simulate 4 types of fuel
    Q = np.diag([1e-4, 1e-5, 1e-6, 1e-6, 1e-6])
    P0 = np.diag([0.01, 0.01, 0.01, 0.001, 0.001])
    Tk = np.array([1.0, 10.0, 100.0]) * 3600
    dt = 3600
    print("INFO: Time step is %d seconds." % dt)

    # preprocess all static covariates
    X = np.zeros((dom_shape[0], dom_shape[1], 4))
    X[:,:,1] = 1.0
    X[:,:,2] = hgt / 2000.0
    if np.any(rain) > 0.05:
      X[:,:,3] = rain
    else:
      X = X[:,:,:3]

    # load current state (or initialize from equilibrium if not found)
    fm_init = None
    fm_cov_init = None
    try:
      in_file = netCDF4.Dataset('inputs/fm.nc', 'r')
      fm_init = np.zeros(10)
      fm_cov_init = np.zeros(100)
    except Exception:
      print('INFO: input file not found, initializing from equilibrium')
      fm_init = 0.5 * (ed + ew)
      fm_cov_init = P0

    models = GridMoistureModel(fm_init[:,:,np.newaxis][:,:,np.zeros((3,),dtype=np.int)], Tk, P0)

    print('INFO: performing forecast at: [time=%s].' % str(tm))

    # compute the FORECAST
    models.advance_model(ed, ew, rain, dt, Q)
    f = models.get_state()

    # examine the assimilated fields (if assimilation is activated)
    for i in range(3):
      print('INFO [%d]: [min %g, mean %g, max %g]' % (i, np.amin(f[:,:,i]), np.mean(f[:,:,i]), np.amax(f[:,:,i])))
      if np.any(f[:,:,i] < 0.0):
        print("WARN: in field %d there were %d negative moisture values !" % (i, np.count_nonzero(f[:,:,i] < 0.0)))
      if np.any(f[:,:,i] > 2.5):
        print("WARN: in field %d there were %d moisture values above 2.5!" % (i, np.count_nonzero(f[:,:,i] > 2.5)))

    # fit the trend surface model to data
    tsm, tsm_var = fit_tsm(obss, X)
    if np.count_nonzero(tsm > 2.5) > 0:
        print('WARN: in TSM found %d values over 2.5, %d of those had rain, clamped to 2.5' %
                (np.count_nonzero(tsm > 2.5),
                 np.count_nonzero(np.logical_and(tsm > 2.5, rain > 0.0))))
        tsm[tsm > 2.5] = 2.5
    if np.count_nonzero(tsm < 0.0) > 0:
        print('WARN: in TSM found %d values under 0.0, clamped to 0.0' % np.count_nonzero(tsm < 0.0))
        tsm[tsm < 0.0] = 0.0

    # if there were any observations, run the kalman update step
    Kg = np.zeros((dom_shape[0], dom_shape[1], len(Tk)+2))
    models.kalman_update_single2(tsm, tsm_var, 1, Kg)

    # check post-assimilation results
    for i in range(3):
        if np.any(models.get_state()[:,:,i] < 0.0):
            print("WARN: in field %d there were %d negative moisture values !" % (i, np.count_nonzero(models.get_state()[:,:,i] < 0.0)))

    # store post-assimilation (or forecast depending on whether observations were available) FM-10 state and variance
    ncfmc[:,:,:] = models.get_state()
    ncfmc_cov[:,:,:,:] = models.P

    # store the diagnostics in a binary file when done
    diagnostics().dump_store(os.path.join(cfg['output_dir'], 'diagnostics.bin'))

    # close the netCDF file (relevant if we did write into FMC_GC)
    if out_file is not None:
        out_file.close()


if __name__ == '__main__':
#    import profile
#    import pstats
#    profile.run('run_module(); print', 'spatial_model.stats')

#    stats = pstats.Stats('spatial_model.stats')
#    stats.strip_dirs()
#    stats.sort_stats('cumulative')
#    stats.print_stats()

    run_module()
    sys.exit(0)

