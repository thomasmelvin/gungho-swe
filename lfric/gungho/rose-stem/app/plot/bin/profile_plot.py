#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (C) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Plots profiles of horizontally averaged fields.
'''

import sys
import iris

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

if iris.__version__ < "3.0.0":
    iris.FUTURE.netcdf_promote = True

#------------------------------------------------------------------------------#
# Use a variables name to obtain its data from an iris cube
#------------------------------------------------------------------------------#
def load_cube_by_varname(filename, var):

   variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == var))
   return iris.load_cube(filename, constraint=variable_constraint)

#------------------------------------------------------------------------------#
# Set heights of levels for selected extrusion
#------------------------------------------------------------------------------#
def make_extrusion(extrusion, number_of_layers, domain_top):

    if extrusion == 'linear':
        z1d = np.linspace(0.0, domain_top, nlayers+1)
    elif extrusion == 'quadratic':
        z1d = domain_top * np.array([(float(k) / float(number_of_layers)) ** 2
                                      for k in range(number_of_layers+1)])
    elif extrusion == 'geometric':
        stretching = 1.03
        z1d = np.zeros([number_of_layers+1])
        dz = 1.0
        z1d[0] = 0.0
        for k in range(number_of_layers):
            z1d[k+1] = z1d[k] + dz
            dz = dz * 1.03
        z1d = z1d * domain_top / z1d[number_of_layers]
    elif extrusion == 'umL38':
        if number_of_layers != 38:
            raise ValueError('For UM extrusion need 38 layers not %d' % number_of_layers)

        z1d = domain_top * np.array([.0000000, .0005095, .0020380,
                 .0045854, .0081519, .0127373, .0183417, .0249651, .0326074,
                 .0412688, .0509491, .0616485, .0733668, .0861040, .0998603,
                 .1146356, .1304298, .1472430, .1650752, .1839264, .2037966,
                 .2246857, .2465938, .2695209, .2934670, .3184321, .3444162,
                 .3714396, .3998142, .4298913, .4620737, .4968308, .5347160,
                 .5763897, .6230643, .6772068, .7443435, .8383348, 1.000000])

    else:
        raise NotImplementedError('Extrusion %s is not implemented' % extrusion)

    return z1d

#################################################################################

def make_figures(data_path, field_list, extrusion, number_of_layers, domain_top, plot_path):

    z1d = make_extrusion(extrusion, number_of_layers, domain_top)

    file_name = data_path+'/lfric_diag.nc'

    for field_name in field_list:

        field_cube = load_cube_by_varname(file_name, field_name)

        field_shape = field_cube.shape

        if len(field_shape) == 3:
           (ntims, nlevs, npts) = field_shape
        else:
           (nlevs, npts) = field_shape
           ntims = 1

        nz1d = len(z1d)
        if (nlevs != nz1d):
            print("profile_plot.py: only works for fields in Wtheta")
            raise ValueError("Number of levels in file does't match specified number_of_layers %d"
                             % number_of_layers)

        for itim in range(ntims):
            if len(field_shape) == 3:
                field = field_cube.data[itim, :, :]
            else:
                field = field_cube.data[:,:]

            field_prof = np.zeros(nlevs)
            for k in range(nlevs):
                field_prof[k] = np.sum(field[k, :]) / float(npts)

            plt.plot(field_prof, z1d)

        plt.xlabel('Domain mean '+field_name, style='italic')
        plt.ylabel('height (km)', style='italic')
        plt.savefig(plot_path+'/'+field_name+'_profile.png')
        plt.close()

if __name__ == "__main__":

    try:
        data_path, fields, extrusion, number_of_layers, domain_top, plot_path = sys.argv[1:7]
    except ValueError:
        print(
            "Usage: {0} <data_path> <field1:field2:...> <extrusion> <number_of_layers>"
            "<domain_top(km)> <plot_path>".format(sys.argv[0]))
        exit(1)

    field_list = fields.split(':')

    make_figures(data_path, field_list, extrusion, int(number_of_layers), float(domain_top), plot_path)
