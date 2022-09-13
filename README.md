# gungho-swe
Mixed FEM model of the shallow water equations

This is a reduced version of the full LFRic model code

###
Copyright of code existing on 3rd April 2017:
Copyright (c) 2017, Met Office, on behalf of HMSO and Queen's Printer
###
Copyright of code created subsequently to 3rd April 2017:
(C) Crown copyright 2017 Met Office. All rights reserved.
### 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

###
User Guide
###

The shallow water model is part of LFRic. The code in this repository will only run the shallow water model.

LFRic can be built on any Linux like system or virtual machine. 
LFRic has a number of dependencies on third-party software. 
These libraries and tools have to be built before LFRic can be successfully built and run.
The following set of instructions documents a build of LFRic on a clean install of Ubuntu v16.04. 
Your system may be set up slightly differently, so the installation may differ in some of the details, 
but the overall sequence of tasks should be the same.

The following notes assume you already have a fully installed (and working) 
copy of Ubuntu 16.04 running and that you have access to the command line (through a terminal window). 
The versions of the software packages used below are just those used to build LFRic on the Met Office Linux desktop system in June 2017. 
Other versions (particularly newer versions) are also likely to work.

To run the SW model you need to do the following:

* Install a recent GNU compiler (ftp://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-8.1.0/gcc-8.1.0.tar.gz)
* Download and install mpich (http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz)
* Install dependencies (sudo apt-get install m4 zlib1g-dev)
* Download and install hdf5 (http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8/hdf5-1.8.16/src/hdf5-1.8.16.tar.gz)
* Download and install netcdf (ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz)
* Install netcdf fortran binding (ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.2.tar.gz)
* Install netcdf C++ binding (ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx-4.2.tar.gz)
* Download and install YAXT (https://www.dkrz.de/redmine/attachments/download/498/yaxt-0.9.0.tar.gz)
* Download and install subversion (sudo apt-get install subversion)
* Download and install XIOS (svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk@2252 XIOS)
* Install pFUnit (https://sourceforge.net/projects/pfunit/files/latest/download/pFUnit-3.2.9.tgz)
* Install Jinja2 (sudo pip install Jinja2)
* Install PSyclone (sudo pip install PSyclone==2.3.1)
* Note: this requires Python 3. 

Once these are installed you need the following: 

* Go to directory gungho-swe/lfric/miniapps/shallow_water
* Type "make build" to compile the shallow water model
* Go to directory gungho-swe/lfric/miniapps/shallow_water/example
* Type "../bin/shallow_water configuration.nml" to run the shallow water model using the settings given in configuration.nml (found in this directory)
* The configurations for other tests can be found in gungho-swe/lfric/miniapps/shallow_water/configurations
* Plotting scripts can be found in gungho-swe/lfric/miniapps/shallow_water/rose-stem/app/plot/bin
* New meshes can be generated from gungho-swe/lfric/mesh_tools/ using "make build". In the example directory type "../bin/cubedsphere_mesh_generator mesh_generation.nml" to create the mesh specified in mesh_generation.nml. The configuration files for the shallow water meshes can be found in gungho-swe/lfric/miniapps/shallow_water/configurations
