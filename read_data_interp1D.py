import numpy as np
import pandas as pd
import sys
import calendar

import numpy.ma as ma

from glob import glob
from datetime import date, datetime, timedelta

import matplotlib.pyplot as plt
import matplotlib as mpl # in python
from matplotlib.colors import BoundaryNorm
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import maskoceans
from netCDF4 import Dataset
from scipy import interpolate

from rpn.rpn import RPN
from rpn.domains.rotated_lat_lon import RotatedLatLon


'''
Read the dm files from the Samples folder of the simulations.
- variables: GZ, UU, VV, TT, HU

For each variable, read the vertical profile and transform it into levels that will be more usefull:
- height: [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]

Estimate the pressure for each level;
Estimate the density for each level;
Calculate the Inversion Strength (deltaT) at 40m

This version gets a station point and do the calculations.

Saving the results in CSV (if doing for only one point) or Netcdf (if for the whole domain)

for later:
    - use the Kahn algorithm, calculate the inversion frequency, height and strength.


'''

def main():
  datai = 1989
  dataf = 1990

    # constants
  Rd = 287  # Gas constant of dry air
  g = 9.80665 # gravity

  lats = []
  lons = []
  stnames = []

  stations = open('stations.txt', 'r')
  for line in stations:
    aa = line.replace("\n", '').split(';')
    if (aa[0] != "#"):      
      lats.append(float(aa[3]))
      lons.append(float(aa[5]))
      stnames.append(aa[1].replace(',',"_"))

  stations.close()
  # simulation
  exp = "cPanCan_011deg_675x540_SPN_ERA5_90lvl"
  #exp = "cPanCan_011deg_675x540_SPN_ERA5_80lvl"
  #exp = "cPanCan_011deg_675x540_SPN_CanESM2_histo_r1i1p1_90lvl"

  folder = "/home/cruman/projects/rrg-sushama-ab/cruman/storage_model/Output/{0}".format(exp)

  # Sounding Data
  #sounding_file = "/home/cruman/project/cruman/Scripts/soundings/inv_list_DJF.dat"

  period = ["DJF", "JJA", 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Nov', 'Dec']
  period = ["DJF", "JJA", "SON", "MAM"] #, "JFM", "JAS"]
  period = ["DJF", "JJA"]
  period = ["DJF"]

  height = [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]

  for per in period:

    #read file
    month_range = getMonths(per)

    for month in month_range:
      for year in range(datai, dataf+1):  
        print(month, year)

        file_gem = "{0}/Samples/{1}_{2}{3:02d}/dm*".format(folder, exp, year, month)

        file_list = glob(file_gem)
        file_list.sort()

        ini = True
        for f in file_list:

          with RPN(f) as r:

            # position of elements: (time, z, x, y)
            # top altitude is the 0 array index
            gz = r.variables['GZ'][:]
            # [...0.982782, 0.988501, 0.994254, 0.997123, 1]
            # [...,gz_uu, gz_tt, gz_uu, gz_tt, surface]

            tt = r.variables['TT'][:]
            # [..., 0.988501, 0.997123, 1.5]

            uu = r.variables['UU'][:]
            vv = r.variables['VV'][:]
            # [..., 0.982782, 0.994254, 10m]

            hu = r.variables['HU'][:] # specific humidity

            gz_0 = gz[:,-1,:,:]

            # removing the last level (surface)
            gz = gz[:,:-1,:,:]

            # spliting between tt/hu and uu levels
            gz_tt = gz[:,1::2,:,:]
            gz_uu = gz[:,::2,:,:]

            tt_0 = tt[:,-1,:,:]+273.15
            tt = tt[:,:-1,:,:]+273.15

            hu_0 = hu[:,-1,:,:]
            hu = hu[:,:-1,:,:]

            uv = np.sqrt(np.power(uu,2) + np.power(vv,2))
            uv_0 = uv[:,-1,:,:]
            uv = uv[:,:-1,:,:]

            uu_0 = uu[:,-1,:,:]
            uu = uu[:,:-1,:,:]

            vv_0 = vv[:,-1,:,:]
            vv = vv[:,:-1,:,:]

            mslp = r.variables['PN'][:]
            aux = r.variables['PN']
            
            if ini:
              ini = False
              data_uu = uu
              data_vv = vv
              data_uv = uv
              data_tt = tt
              data_hu = hu
              data_mslp = mslp
              data_gz_tt = gz_tt
              data_gz_uu = gz_uu

              dates = [str(d) for d in aux.sorted_dates]

              lons2d, lats2d = r.get_longitudes_and_latitudes_for_the_last_read_rec()
            else:
              data_uu = np.vstack( (data_uu, uu) )
              data_vv = np.vstack( (data_vv, vv) )
              data_uv = np.vstack( (data_uv, uv) )
              data_tt = np.vstack( (data_tt, tt) )
              data_hu = np.vstack( (data_hu, hu) )
              data_mslp = np.vstack( (data_mslp, mslp) )
              data_gz_tt = np.vstack( (data_gz_tt, gz_tt) )
              data_gz_uu = np.vstack( (data_gz_uu, gz_uu) )

              dates += [str(d) for d in aux.sorted_dates]
        
        print("Initiating calculations on pre-selected locations")
        for lat, lon, name in zip(lats, lons, stnames):
          i, j = geo_idx([lat, lon], np.array([lats2d, lons2d]))

          # Virtual Temperature
          # Tv ~ T*(1 + 0.61*w)
          # Tv ~ T*(1 + 0.61*(hu/(1-hu)))
          Tv = data_tt[:,:,i,j]*(1 + 0.61*(data_hu[:,:,i,j]/(1-data_hu[:,:,i,j])))
          #Tv_0 = tt_0*(1 + 0.61*(hu_0/(1-hu_0)))

          # Using the hypsometric equation
          # Z2 - Z1 = (Rd*Tv)*ln(p1/p2)/g
          # at sea level, z1 = 0 and p1 = mslp
          # mslp/p2 = np.exp(Z2*g/(Rd*Tv))
          # p2 = mslp/np.exp(Z2*g/(Rd*Tv))
          p = data_mslp[:,:,i,j]/(np.exp(data_gz_tt[:,:,i,j]*g/(Rd*Tv)))

            # estimating the air density
          pho = p/(Rd*Tv)
          #pho_0 = mslp/(Rd*Tv_0)

          # interpolate values to nice levels
          #for dd, i in enumerate(dates):

          tt_i = interpData(data_gz_tt[:,:,i,j], height, data_tt[:,:,i,j])

          print(tt_i[0])
          print(data_tt[0,:,i,j])

            
          sys.exit()

def geo_idx(dd, dd_array, type="lat"):
  '''
    search for nearest decimal degree in an array of decimal degrees and return the index.
    np.argmin returns the indices of minium value along an axis.
    so subtract dd from all values in dd_array, take absolute value and find index of minimum.
    
    Differentiate between 2-D and 1-D lat/lon arrays.
    for 2-D arrays, should receive values in this format: dd=[lat, lon], dd_array=[lats2d,lons2d]
  '''
  if type == "lon" and len(dd_array.shape) == 1:
    dd_array = np.where(dd_array <= 180, dd_array, dd_array - 360)

  if (len(dd_array.shape) < 2):
    geo_idx = (np.abs(dd_array - dd)).argmin()
  else:
    if (dd_array[1] < 0).any():
      dd_array[1] = np.where(dd_array[1] <= 180, dd_array[1], dd_array[1] - 360)

    a = abs( dd_array[0]-dd[0] ) + abs(  np.where(dd_array[1] <= 180, dd_array[1], dd_array[1] - 360) - dd[1] )
    i,j = np.unravel_index(a.argmin(), a.shape)
    geo_idx = [i,j]

  return geo_idx

def interpData(levels, new_levels, data, interp='linear'):
  """
    Interpolate data to custom levels
    levels: Original level
    custom_levels: new level
    data: Original variable to be interpolated to custom pressure level
    returns: new_val, the original variable interpolated.
  """
  new_val = np.zeros_like(new_levels)

  f = interpolate.interp1d(levels, data, kind=interp)

  for level in range(new_val.shape[0]):
    new_val[level] = f(new_levels[level])

  return new_val


def getMonths(period):

  # 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Nov', 'Dec'
  if period == "DJF":
    months = [12, 1, 2]
  elif period == "JJA":
    months = [6, 7, 8]
  elif period == "MAM":
    months = [3, 4, 5]
  elif period == "SON":
    months = [9, 10, 11]
  elif period == "JFM":
    months = [1, 2, 3]
  elif period == "JAS":
    months = [7, 8, 9]
  elif period == "Jan":
    months = [1]
  elif period == "Feb":
    months = [2]
  elif period == "Mar":
    months = [3]
  elif period == "Apr":
    months = [4]
  elif period == "May":
    months = [5]
  elif period == "Jun":
    months = [6]
  elif period == "Jul":
    months = [7]
  elif period == "Aug":
    months = [8]
  elif period == "Sep":
    months = [9]
  elif period == "Oct":
    months = [10]
  elif period == "Nov":
    months = [11]
  elif period == "Dec":
    months = [12]

  return months

if __name__ == "__main__":
    main()
