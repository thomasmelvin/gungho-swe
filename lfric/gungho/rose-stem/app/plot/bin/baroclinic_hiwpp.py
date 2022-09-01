#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2018 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Python script to plot xz slices along the y=0 and xy slices on a specified
level

This version takes nodal format output files and
interpolates onto a regular grid.

Filename hardcoded.

Levels are determined from the data

This version stitches together a directory of files
and extracts all levels so it can work in the serial
case where there is one file or the parallel case where
there is a file for each processor.

'''

import numpy as np
from scipy.interpolate import griddata
import sys
from magma import magma_data
from read_data import read_nodal_data
# Need to set a non-interactive backend for suites
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# Use magma colormap
from matplotlib.colors import ListedColormap


magma = ListedColormap(magma_data, name='magma')
plt.register_cmap(name='magma', cmap=magma)

levels = None
data = None


def make_figure(plotpath, field, component, timestep, levels, plotlevel):

    val_col = 'c' + str(component)

    # get min and max of x,y data for plot axes

    min_lev = min(levels)

    deltaz = data.loc[data['level'] == (min_lev+1)][val_col].min() \
        - data.loc[data['level'] == (min_lev)][val_col].min()
    xmin = data.loc[data['level'] == min_lev]['x'].min()
    xmax = data.loc[data['level'] == min_lev]['x'].max()
    ymin = data.loc[data['level'] == min_lev]['y'].min()
    ymax = data.loc[data['level'] == min_lev]['y'].max()

    zmin = min(levels)*1000.0
    zmax = max(levels)*1000.0

    r2d = 180.0/np.pi
    nx, ny, nz = 720, 360, len(levels)

    # Create 2D plot
    x2d = np.linspace(xmin, xmax, nx)
    z2d = np.linspace(zmin, zmax, nz)
    y2d = np.linspace(ymin, ymax, ny)
    zi = np.zeros([ny, nx, len(levels)])
    xi, yi = np.meshgrid(x2d, y2d)
    for p in xrange(len(levels)):
        p_data = data.loc[data['level'] == levels[p]]
        zi[:, :, p] = griddata((p_data['x'].values, p_data['y'].values),
                               p_data[val_col].values, (xi, yi),
                               method='linear')

    # Specific settings for plotting theta and p
    if field == 'theta':
        cc = np.linspace(220, 330, 12)
        plotlevel = 1
    elif field == 'exner':
        rd = 287.05
        p0 = 100000.0
        kappa = rd/1005.0
        zi = 0.01*zi**(1.0/kappa) * p0
        cc = np.linspace(916, 1020, 14)
        plotlevel = 0
    elif field == 'w3projection_u2':
        cc = np.linspace(-25, 40, 14)
    elif field == 'w3projection_u3':
        plotlevel = 1
        cc = np.linspace(-0.01, 0.01, 11)
    else:
        p_data = data.loc[data['level'] == levels[plotlevel]]
        cc = np.linspace(np.amin(p_data[val_col].values),
                         np.amax(p_data[val_col].values), 13)

    fig = plt.figure(figsize=(10, 5))
    # Extrapolate exner field to surface
    if field == 'exner':
        dz = 1.5*zi[:, :, 0] - 0.5*zi[:, :, 1]
    else:
        dz = zi[:, :, plotlevel]

    plt.clf()
    cf = plt.contourf(xi * r2d, yi * r2d, dz, cc, cmap=magma)
    plt.colorbar(cf, cmap=magma)
    cl = plt.contour(xi * r2d, yi * r2d, dz, cc, linewidths=0.5, colors='k')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    if field == 'w3projection_u3':
        plt.axis([130, 140, -47.5, -22.5])

    if (field != 'u' and field != 'xi'):
        out_file_name = plotpath + "/" "baroclinic_xy_" + field + "_" \
            + timestep + ".png"
    else:
        out_file_name = plotpath + "/" "baroclinic_xy_" + field \
            + str(component) + "_" + timestep + ".png"

    plt.savefig(out_file_name, bbox_inches='tight')


if __name__ == "__main__":

    try:
        config, datapath, fields, timesteps, plotpath, plotlevel \
            = sys.argv[1:7]
    except ValueError:
        print("Usage: {0} <config> <datapath> <fields_list> <timestep_list>"
              "<plotpath> <plotlevel>".format(sys.argv[0]))
        exit(1)

    # Split out the list of fields
    field_list = fields.split(':')

    # Split out the list of timesteps
    ts_list = timesteps.split(':')

    any_plots = False

    for field in field_list:

        if field in ['rho', 'theta', 'exner']:
            # Scalar fields
            ncomp = 1
            comp = 1
        else:
            # Vector fields
            ncomp = 3
            # W3 projected U, V, W and XI components
            if field in ['w3projection_u1', 'w3projection_u2',
                         'w3projection_u3', 'w3projection_xi1',
                         'w3projection_xi2', 'w3projection_xi3']:
                comp = 1
            elif (field == 'u' or field == 'xi'):
                comp = [1, 2, 3]

        for ts in ts_list:

            filestem = datapath + "/" + config + "_nodal_" + field + "_" \
                + ts + "*"

            if (field != 'u' and field != 'xi'):
                data = read_nodal_data(filestem, ncomp, comp)
                if (not data.empty):
                    levels = data.level.unique()
                    make_figure(plotpath, field, comp, ts, levels,
                                int(plotlevel))
                    any_plots = True

            else:
                for comp_u in comp:
                    data = read_nodal_data(filestem, ncomp, comp_u)
                    if (not data.empty):
                        levels = data.level.unique()
                        make_figure(plotpath, field, comp_u, ts, levels,
                                    int(plotlevel))
                        any_plots = True

    if not any_plots:
        print("Error: No plots made.")
        exit(2)
