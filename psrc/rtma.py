
import netCDF4
import numpy as np
import os

from datetime import datetime
import pytz
import glob


def check_values_in_range(name,vals,vmin,vmax):
  if np.any(vals[:] < vmin):
    print('WARNING: Var [%s] has %d values below minimum %g, clamping ...' % (name, np.sum(vals[:]<vmin), vmin))
    vals[vals < vmin] = vmin
  if np.any(vals[:] > vmax):
    print('WARNING: Var [%s] has %d values above maximum %g, clamping ...' % (name, np.sum(vals[:]>vmax), vmax))
    vals[vals > vmax] = vmax


def find_region_indices(glat,glon,minlat,maxlat,minlon,maxlon):
  """
  Find the indices i1:i2 (lat dimension) and j1:j2 (lon dimension)
  that contain the desired region (minlat-maxlat,minlon-maxlon).
  """
  i1, i2, j1, j2 = 0, glat.shape[0], 0, glat.shape[1]
  done = False
  while not done:
    done = True
    tmp = np.where(np.amax(glat[:, j1:j2],axis=1) < minlat)[0][-1]
    if i1 != tmp:
      i1 = tmp
      done = False
    tmp = np.where(np.amin(glat[:, j1:j2],axis=1) > maxlat)[0][0]
    if i2 != tmp:
      i2 = tmp
      done = False
    tmp = np.where(np.amax(glon[i1:i2,:],axis=0) < minlon)[0][-1]
    if j1 != tmp:
      j1 = tmp
      done = False
    tmp = np.where(np.amin(glon[i1:i2,:],axis=0) > maxlon)[0][0]
    if j2 != tmp:
      j2 = tmp
      done = False
  return i1,i2,j1,j2


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
        return None

      data[name] = d.variables[rtmav][0,:,:]
      rtma_time = d.variables['time'][0]
      d.close()

  d = netCDF4.Dataset("static/rtma_elevation.nc")
  data['HGT'] = d.variables['HGT_surface'][0,:,:]
  Lat = d.variables['latitude'][:,:]
  Lon = d.variables['longitude'][:,:] - 360
  data['Lat'] = Lat
  data['Lon'] = Lon

  # chop all datasets to Colorado coordinates [36.9 <-> 41.1, -109.1 <-> -101.9]
  i1, i2, j1, j2 = find_region_indices(Lat,Lon,37,41,-109,-102)
  print('INFO: chopping all grids data to %d-%d x %d-%d' % (i1,i2,j1,j2))

  for k,v in data.iteritems():
    data[k] = v[i1:i2,j1:j2]

  # compute the relative humidity according to NOAA formula [http://www.srh.noaa.gov/images/epz/wxcalc/vaporPressure.pdf]
  t2 = data['T2'] - 273.15
  td = data['DPT'] - 273.15
  data['RH'] = 100*np.exp(17.625*243.04*(td - t2) / (243.04 + t2) / (243.0 + td))

  # check all input values
  check_values_in_range('RH',data['RH'],0,100)
  check_values_in_range('T2',data['T2'],200,350)
  check_values_in_range('DPT',data['DPT'],200,350)

  # add time after slicing
  data['Time'] = datetime.fromtimestamp(rtma_time,tz=pytz.utc)

  print('Loaded RTMA fields for %s.' % str(data['Time']))

  return data



def subset_netcdf_file(ncpath0,ncpath1,minlat,maxlat,minlon,maxlon):
  """
  Subset the netCDF file <ncpath0> to contain data in the (region minlat-maxlat,minlon-maxlon).
  Write the resulting file into <ncpath1>.
  """

  d0 = netCDF4.Dataset(ncpath0)
  d1 = netCDF4.Dataset(ncpath1,'w',format='NETCDF4')

  lat = d0.variables['latitude'][:,:]
  lon = d0.variables['longitude'][:,:] - 360

  i1,i2,j1,j2 = find_region_indices(lat,lon,37,41,-109,-102)

  d1.createDimension('time', None)
  d1.createDimension('x', j2-j1)
  d1.createDimension('y', i2-i1)

  for vname,v0 in d0.variables.iteritems():
    v1 = None
    if hasattr(v0,'_FillValue'):
      v1 = d1.createVariable(vname, v0.dtype, v0.dimensions,fill_value=getattr(v0,'_FillValue'))
    else:
      v1 = d1.createVariable(vname, v0.dtype, v0.dimensions)

    # copy all attributes except fillvalue
    for va in v0.ncattrs():
      if va != '_FillValue':
        setattr(v1,va,getattr(v0,va))

    if len(v0.dimensions) == 2:
      v1[:,:] = v0[i1:i2,j1:j2]
    elif len(v0.dimensions) == 3:
      v1[:,:,:] = v0[:,i1:i2,j1:j2]
    else:
      if vname=='x':
        v1[:] = v0[j1:j2]
      elif vname=='y':
        v1[:] = v0[i1:i2]
      else:
        v1[:] = v0[:]

  d0.close()
  d1.close()


def subset_rtma_analysis(path_pairs):
  for f0,f1 in path_pairs:
    subset_netcdf_file(f0,f1,37,41,-109,-102)


def subset_archive_rtma_analysis(rtma_dir):
  file_list = glob.glob(rtma_dir + "/L*.nc")
  pairs = []

  print file_list

  for f in file_list:
    if 'LEIA' in f:
      pairs.append((f, os.path.join(rtma_dir,'ds.precipa.nc')))
    elif 'LRIA' in f:
      pairs.append((f, os.path.join(rtma_dir,'ds.td.nc')))
    elif 'LTIA' in f:
      pairs.append((f, os.path.join(rtma_dir,'ds.temp.nc')))

  subset_rtma_analysis(pairs)
