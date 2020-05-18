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
#from netCDF4 import Dataset

#from rpn.rpn import RPN
#from rpn.domains.rotated_lat_lon import RotatedLatLon


'''
Read the dm files from the Samples folder of the simulations.
- variables: GZ, UU, VV, TT, HU

For each variable, read the vertical profile and transform it into levels that will be more usefull:
- height: [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]

Estimate the pressure for each level;
Estimate the density for each level;
Calculate the Inversion Strength (deltaT) at 40m

Saving the results in CSV (if doing for only one point) or Netcdf (if for the whole domain)

for later:
    - use the Kahn algorithm, calculate the inversion frequency, height and strength.


'''

def main():
    datai = 1981
    dataf = 1990

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

                file_list = glob.glob(file_gem)
                print(file_list)
                file_list.sort()

                print(file_list)


        #arq_gem.close()
        #arq_gz.close()


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