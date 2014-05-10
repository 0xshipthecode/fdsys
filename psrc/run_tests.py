

from fmncast import run_data_assimilation
from fmrun import run_model
from datetime import datetime, timedelta
from swissknife import get_grid_values_from_file
import os


def time_to_input_dir(tm):
  return 'inputs/%04d%02d%02d-%02d00' % (tm.year,tm.month,tm.day,tm.hour)


if __name__ == '__main__':

  # load the csv file
  with open('CO_100hr2.csv','r') as f:
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

    # for each uday, run the raw simulation
    for tm in udays:
      print('*************** DAY = %s ******************' % str(tm))
      # start 12 hrs back and run 11 times
      tm2 = tm.replace(hour=20)
      tms = [tm2 - timedelta(seconds=i*3600) for i in range(6+1)]
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
            print('VALUES: %s,%s,%g,%g,%g,%g,%g' % (o[0],str(o[1]),o[2],o[3],o[4],valsda[2]*100,valsr[2]*100))
      else:
        print('MISSING: cannot obtain values for time %s' % str(tm))
