

from fmncast import run_data_assimilation
from fmrun import run_model
from datetime import datetime, timedelta
from swissknife import get_grid_values_from_file
import os
import sys


def time_to_input_dir(tm):
  return 'inputs/%04d%02d%02d-%02d00' % (tm.year,tm.month,tm.day,tm.hour)


if __name__ == '__main__':

  if len(sys.argv) != 4:
    print('usage: %s <obs-file> <fm-index> <prerun-length-hrs>' % sys.argv[0])
    sys.exit(0)

  obs_path = sys.argv[1]
  fm_ndx = int(sys.argv[2])
  prerun_hrs = int(sys.argv[3])

  print('INFO: **********************************************************************')
  print('INFO: run_tests.py with obs_file=%s fm-index=%d prerun-length=%d' % (obs_path,fm_ndx,prerun_hrs))
  print('INFO: **********************************************************************')

  # load the csv file
  with open(obs_path,'r') as f:
    obss = []
    days = []
    for line in f.readlines():
      dtstr,name,lats,lons,vals = line.strip().split(',')
      mon,day,yr = map(lambda x: int(x), dtstr.split('/'))
      dt = datetime(2000+yr,mon,day)
      lat,lon,val = float(lats),float(lons),float(vals)
      obss.append([name,dt,lat,lon,val])
      days.append(dt)

    udays = sorted(set(days))
    days2 = map(lambda x: x - timedelta(1), udays)
    days = udays[:]
    days.extend(days2)
    days = sorted(set(days))

    # write all dates to file
    with open('dates_needed','w') as f:
      for d in days:
        f.write('%d,%d,%d\n' % (d.year,d.month,d.day))

    #print('FATAL: terminating after reading in dates')
    #sys.exit(1)

    udays = [datetime(2013,6,11)]

    # for each uday, run the raw simulation
    for tm in udays:

      # remove pre-computed fuel moisture file
      # os.system('find state -name "fm*" | xargs rm')

      print('*************** DAY = %s ******************' % str(tm))
      # start 12 hrs back and run 11 times
      tm2 = tm.replace(hour=20)
      tms =  [tm2 - timedelta(seconds=i*3600) for i in range(prerun_hrs+1)]
      tms.reverse()

      # run the simulation (raw model)
      for i in range(len(tms)-1):
        ind0,ind1 = time_to_input_dir(tms[i]),time_to_input_dir(tms[i+1])
        # run the raw model
        run_model(ind0,ind1,'state')
        # run the simulation with data assimilation
        run_data_assimilation(ind0,ind1,'state')

      # extract the values for the particular observations for the time we just ran
      ncpathr = 'state/fmraw-%04d%02d%02d-2000.nc' % (tm.year,tm.month,tm.day)
      ncpathda = 'state/fm-%04d%02d%02d-2000.nc' % (tm.year,tm.month,tm.day)
      if os.path.exists(ncpathr) and os.path.exists(ncpathda):
        for o in obss:
          if o[1] == tm:
            dist,gi,gj,valsr = get_grid_values_from_file(ncpathr,'FMC_GC_RAW',o[2],o[3])
            dist,gi,gj,valsda = get_grid_values_from_file(ncpathda,'FMC_GC',o[2],o[3])
            print('VALUES,%s,%s,%g,%g,%g,%g,%g' % (o[0],str(o[1]),o[2],o[3],o[4],valsda[fm_ndx]*100,valsr[fm_ndx]*100))
      else:
        print('MISSING: cannot obtain values for time %s' % str(tm))
