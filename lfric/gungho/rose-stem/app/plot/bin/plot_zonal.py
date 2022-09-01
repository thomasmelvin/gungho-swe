#!/usr/bin/env python

# Need to set a non-interactive backend for suites
from __future__ import absolute_import
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import sys
import iris

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata

if iris.__version__ < "3.0.0":
    iris.FUTURE.netcdf_promote = True

# Size of regular grid
ny, nx = 200, 400

def make_figures(filein, plotpath, fields):

   if fields is None:
      # Set the standard default fields
      fields = ['air_potential_temperature', 'eastward_wind']

   for field in fields:
      cube = iris.load_cube(filein, field)
      # Vertical levels will be last entry in dimension coords
      levels_name = cube.dim_coords[-1].name()
      #Set some levels:
      levels=None
      if field=='air_potential_temperature':
         levels=np.arange(270,360,5.)
      if field=='eastward_wind':
         levels=np.arange(-20,36,4.)

      n_levs = len(cube.coord(levels_name).points)

      zonal_data=np.zeros((ny,n_levs))

      time = np.around(cube.coord('time').points, decimals=1)

      interp_fig = plt.figure(figsize=(20,15))

      x = np.around(cube.coord('longitude').points, decimals=5)
      y = np.around(cube.coord('latitude').points, decimals=5)

      xmin = np.amin(x)
      xmax = np.amax(x)
      ymin = np.amin(y)
      ymax = np.amax(y)

      # Generate a regular grid to interpolate the data.
      xi = np.linspace(xmin, xmax, nx)
      yi = np.linspace(ymin, ymax, ny)

      xf, yf = np.meshgrid(xi, yi)
      
      for l in range(n_levs):
         data = cube.data[:,l].mean(axis=0)

         dmin  = min(data)
         dmax = max(data)

         # Interpolate using delaunay triangularization 
         zi = griddata((x, y), data, (xf, yf), method='linear')
         zi_n = griddata((x, y), data, (xf, yf), method='nearest')
         zi[np.isnan(zi)] = zi_n[np.isnan(zi)] 

         zonal_data[:,l]=zi.mean(axis=1)

      ax = interp_fig.add_subplot(1,1,1)

      ys=np.tile(yi,(n_levs,1))

      CS=plt.contour(zonal_data.T, levels=levels)
      plt.colorbar()
      plt.clabel(CS, CS.levels[1::2], fontsize=15, inline=1,
                 fmt='%2.0f')
      
      plt.suptitle(field, fontsize=25 )

      pngfile='%s/zonal-%s.png' % (plotpath, field)
      plt.savefig(pngfile)
      plt.close()


if __name__ == "__main__":

  try:
     args=sys.argv[:]
     filein, plotpath = args[1:3]
     field_list=None
     if len(args[:])>3: field_list=args[3].split(':')
  except ValueError:
     print("Usage: {0} <filein> <plotpath> [<fields_list>]".format(sys.argv[0]))
     exit(1)
  
  make_figures(filein, plotpath, field_list)
