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
      # check & remove non-sense zero observations
      if oo[8] > 0:
        obss.append(Observation(obst,oo[6],oo[7],oo[8],oo[9],(i,j)))
  return obss


def check_values_in_range(name,vals,vmin,vmax):
  if np.any(vals[:] < vmin):
    print('WARNING: Var [%s] has %d values below minimum %g, clamping ...' % (name, np.sum(vals[:]<vmin), vmin))
    vals[vals < vmin] = vmin
  if np.any(vals[:] > vmax):
    print('WARNING: Var [%s] has %d values above maximum %g, clamping ...' % (name, np.sum(vals[:]>vmax), vmax))
    vals[vals > vmax] = vmax


def time_from_dir(dirname):
  ts = dirname[dirname.index("/")+1:]
  return datetime(int(ts[:4]),int(ts[4:6]),int(ts[6:8]),int(ts[9:11]), 0, 0)


def load_rtma_data(in_dir):
  """
  This function expects the files ds.temp.nc, ds.td.nc, ds.precipa.nc
  to be in the directory inputs and the file rtma_elevation.nc to be in
  the directory static.  All of the data is loaded and stored in a dictionary.
  """
  # decode the time from the input directory name, convention is YYYYMMDD-HH00
  ts = in_dir[in_dir.index("/")+1:]
  suffix = ts[:8] + ts[9:11] + "00.nc"
  print suffix
  file_list = [('ds.precipa.nc', 'LEIA98_KWBR_' + suffix, 'APCP_surface',      'RAIN'),
               ('ds.td.nc',      'LRIA98_KWBR_' + suffix, 'DPT_2maboveground', 'DPT'),
               ('ds.temp.nc',    'LTIA98_KWBR_' + suffix, 'TMP_2maboveground', 'T2')]

  # FIXME: why does the precipitation have the incorrect timestamp?

  rtma_time = None
  data = {}
  for rtmaf, rtmaf2, rtmav, name in file_list:
      d = None
      path = os.path.join(in_dir, rtmaf)
      path2 = os.path.join(in_dir, rtmaf2)
      if os.path.isfile(path):
        d = netCDF4.Dataset(path)
      elif os.path.isfile(path2):
        d = netCDF4.Dataset(path2)
      else:
        print('ERROR: no input file found for variable %s' % name)
        sys.exit(2)

      data[name] = d.variables[rtmav][0,:,:]
      rtma_time = d.variables['time'][0]
      d.close()

  # converted accumulated precip into rain intensity in mm / hr ?
  data['RAIN'] = 0.01 * data['RAIN']

  d = netCDF4.Dataset("static/rtma_elevation.nc")
  data['HGT'] = d.variables['HGT_surface'][0,:,:]
  Lat = d.variables['latitude'][:,:]
  Lon = d.variables['longitude'][:,:] - 360
  data['Lat'] = Lat
  data['Lon'] = Lon

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

  # compute the relative humidity according to NOAA formula [http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf]
  T2c = data['T2'] - 273.15
  DTc = data['DPT'] - 273.15
  e  = 6.11 * 10**(7.5 * DTc / (DTc + 237.3))
  es = 6.11 * 10**(7.5 * T2c / (T2c + 237.3))
  data['RH'] = e / es * 100

  # check all input values
  check_values_in_range('RH',data['RH'],0,100)
  check_values_in_range('T2',data['T2'],250,320)
  check_values_in_range('DPT',data['DPT'],250,310)

  # add time after slicing
  data['Time'] = datetime.fromtimestamp(rtma_time,tz=pytz.utc)

  print('Loaded RTMA fields for %s.' % str(data['Time']))

  return data


def compute_equilibria(T,H):
  """
  Computes atmospheric drying/wetting moisture equilibria from the temperature [K]
  and relative humidity [%].
  """
  d = 0.924*H**0.679 + 0.000499*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
  w = 0.618*H**0.753 + 0.000454*np.exp(0.1*H) + 0.18*(21.1 + 273.15 - T)*(1 - np.exp(-0.115*H))
  d *= 0.01
  w *= 0.01
  return d, w


def esmf_time(dt):
  """
  Format a datetime in ESMF format.
  """
  return '%04d-%02d-%02d_%02d:%02d:%02d' % (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)


def run_data_assimilation(in_dir0, in_dir1, fm_dir):

  # load RTMA data for previous
  print("INFO: loading RTMA data for time t-1 from [%s] ..." % in_dir0)
  tm0 = time_from_dir(in_dir0)
  tm = time_from_dir(in_dir1)
  max_back = 6
  data0 = None
  while max_back > 0:
    in_dir0 = 'inputs/%04d%02d%02d-%02d00' % (tm0.year, tm0.month, tm0.day, tm0.hour)
    print('INFO: searching for RTMA in directory %s' % in_dir0)
    data0 = load_rtma_data(in_dir0)
    if data0 is not None:
      break
    print("WARN: cannot find RTMA data for time %s in directory %s, going back one hour" % (str(tm0),in_dir0))
    max_back -= 1
    tm0 = tm0 - timedelta(0,3600)

  if data0 is None:
    print("FATAL: cannot find a suitable previous RTMA analysis fror time %s." % str(tm))
    return

  print("INFO: loading RTMA data for time t from [%s] ..." % in_dir1)
  data1 = load_rtma_data(in_dir1)

  if data1 is None:
    print('FATAL: insufficient environmnetal data for time %s, skipping ...' % tm)
    return

  # retrieve variables from RTMA
  lat, lon, hgt = data0['Lat'], data0['Lon'], data0['HGT']

  t20, relh0 = data0['T2'], data0['RH']
  t21, relh1, rain = data1['T2'], data1['RH'], data1['RAIN']
  ed0, ew0 = compute_equilibria(t20,relh0)
  ed1, ew1 = compute_equilibria(t21,relh1)
  tm0, tm = data0['Time'], data1['Time']
  tm_str = tm.strftime('%Y%m%d-%H00')
  tm_str0 = tm0.strftime('%Y%m%d-%H00')

  # compute mean values for the Equilibria at t-1 and at t
  ed = 0.5 * (ed0 + ed1)
  ew = 0.5 * (ew0 + ew1)

  dom_shape = lat.shape
  print('INFO: domain size is %d x %d grid points.' % dom_shape)
  print('INFO: domain extent is lats (%g to %g) lons (%g to %g).' % (np.amin(lat),np.amax(lat),np.amin(lon),np.amax(lon)))
  print('INFO: stepping from time %s to time %s' % (tm0, tm))

  # initialize output file
  out_fm_file = os.path.join(fm_dir, 'fm-%s.nc' % tm_str)
  out_file = netCDF4.Dataset(out_fm_file, 'w')
  out_file.createDimension('fuel_moisture_classes_stag', 5)
  out_file.createDimension('south_north', dom_shape[0])
  out_file.createDimension('west_east', dom_shape[1])
  ncfmc_fc = out_file.createVariable('FMC_GC_FC', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
  ncfmc_an = out_file.createVariable('FMC_GC', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
  nckg = out_file.createVariable('K', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag'))
  ncfmc_cov = out_file.createVariable('FMC_COV', 'f4', ('south_north', 'west_east','fuel_moisture_classes_stag', 'fuel_moisture_classes_stag'))
  ncrelh = out_file.createVariable('RELH','f4', ('south_north', 'west_east'))
  ncrelh[:,:] = relh1
  nctemp = out_file.createVariable('T2','f4', ('south_north', 'west_east'))
  nctemp[:,:] = t21
  nclat = out_file.createVariable('Lat', 'f4', ('south_north', 'west_east'))
  nclat[:,:] = lat
  nclon = out_file.createVariable('Lon', 'f4', ('south_north', 'west_east'))
  nclon[:,:] = lon

  print('INFO: opened %s and wrote XLAT,XLONG,RELH,T2 fields.' % out_fm_file)

  ### Load observation data from the stations

  # compute the diagonal distance between grid points
  grid_dist_km = great_circle_distance(lon[0,0], lat[0,0], lon[1,1], lat[1,1])
  print('INFO: diagonal distance in grid is %g' % grid_dist_km)

  raws_path = os.path.join(in_dir1, 'raws_ingest_%4d%02d%02d-%02d%02d.csv' % (tm.year,tm.month,tm.day,tm.hour,tm.minute))
  obss = load_raws_observations(raws_path,lat,lon)
  print('INFO: Loaded %d observations.' % (len(obss)))

  # set up parameters
  Nk = 3  # we simulate 4 types of fuel
  Q = np.diag([1e-4, 1e-5, 1e-6, 1e-6, 1e-6])
  P0 = np.diag([0.01, 0.01, 0.01, 0.001, 0.001])
  Tk = np.array([1.0, 10.0, 100.0])
  dt = 3600
  print("INFO: Time step is %d seconds." % dt)

  # preprocess all covariates
  X = np.zeros((dom_shape[0], dom_shape[1], 4))
  X[:,:,1] = 1.0
  X[:,:,2] = hgt / 2000.0
  if np.any(rain) > 0.01:
    X[:,:,3] = rain
  else:
    X = X[:,:,:3]

  # load current state (or initialize from equilibrium if not found)
  fm0 = None
  fm_cov0 = None
  in_fm_file = os.path.join(fm_dir, 'fm-%s.nc' % tm_str0)
  if os.path.isfile(in_fm_file):
    in_file = netCDF4.Dataset(in_fm_file, 'r')
    fm0 = in_file.variables['FMC_GC'][:,:,:]
    fm_cov0 = in_file.variables['FMC_COV'][:,:,:,:]
    print('INFO: found input file %s, initializing from it [fm is %dx%dx%d, fm_cov is %dx%dx%dx%d]' %
      (in_fm_file,fm0.shape[0],fm0.shape[1],fm0.shape[2],fm_cov0.shape[0],fm_cov0.shape[1],
       fm_cov0.shape[2],fm_cov0.shape[3]))
    in_file.close()
  else:
    print('INFO: input file %s not found, initializing from equilibrium' % in_fm_file)
    fm0 = 0.5 * (ed + ew)
    fm0 = fm0[:,:,np.newaxis][:,:,np.zeros((5,),dtype=np.int)]
    fm0[:,:,3] = -0.04
    fm0[:,:,4] = 0
    fm_cov0 = P0

  models = GridMoistureModel(fm0, Tk, 0.08, 2, 0.6, 7, fm_cov0)

  print('INFO: performing forecast at: [time=%s].' % str(tm))

  # compute the FORECAST
  models.advance_model(ed, ew, rain, dt, Q)
  f = models.get_state()
  ncfmc_fc[:,:,:] = f

  # examine the assimilated fields (if assimilation is activated)
  for i in range(3):
    print('INFO [%d]: [min %g, mean %g, max %g]' % (i, np.amin(f[:,:,i]), np.mean(f[:,:,i]), np.amax(f[:,:,i])))
    if np.any(f[:,:,i] < 0.0):
      print("WARN: in field %d there were %d negative moisture values !" % (i, np.count_nonzero(f[:,:,i] < 0.0)))
    if np.any(f[:,:,i] > 2.5):
      print("WARN: in field %d there were %d moisture values above 2.5!" % (i, np.count_nonzero(f[:,:,i] > 2.5)))

  if len(obss) > 0:

    print('INFO: running trend surface model ...')

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

    print('INFO: running KF ...')

    # run the kalman update step
    Kg = np.zeros((dom_shape[0], dom_shape[1], len(Tk)+2))
    models.kalman_update_single2(tsm[:,:,np.newaxis], tsm_var[:,:,np.newaxis,np.newaxis], 1, Kg)

    # check post-assimilation results
    for i in range(3):
        if np.any(models.get_state()[:,:,i] < 0.0):
            print("WARN: in field %d there were %d negative moisture values !" % (i, np.count_nonzero(models.get_state()[:,:,i] < 0.0)))

    nckg[:,:,:] = Kg

  print('INFO: storing results in netCDF file %s.' % out_fm_file)

  # store post-assimilation (or forecast depending on whether observations were available) FM-10 state and variance
  ncfmc_an[:,:,:] = models.get_state()
  ncfmc_cov[:,:,:,:] = models.get_state_covar()

  # close the netCDF file (relevant if we did write into FMC_GC)
  out_file.close()

  print('INFO: SUCCESS')


if __name__ == '__main__':
#    import profile
#    import pstats
#    profile.run('run_module(); print', 'spatial_model.stats')

#    stats = pstats.Stats('spatial_model.stats')
#    stats.strip_dirs()
#    stats.sort_stats('cumulative')
#    stats.print_stats()

  if len(sys.argv) != 4:
    print('Usage: %s <in_dir0> <in_dir1> <fm_dir>' % sys.argv[0])
    sys.exit(1)

  run_data_assimilation(sys.argv[1],sys.argv[2],sys.argv[3])
  sys.exit(0)
