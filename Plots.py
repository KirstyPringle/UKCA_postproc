# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""

import numpy as np
import iris
import sys
sys.path.append('/nfs/a107/eejvt/PYTHON_CODE')
 
from glob import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
from scipy.io import netcdf
import datetime
import UKCA_lib as ukl
import getpass
#test svn commits
#second test svn commit
#third test commit 27/04
#test
username=getpass.getuser()

plt.interactive(0)
files_directory='/nfs/a107/earkpr/DataVisualisation/UKCA/tebxd/'
files_directory='/nfs/a201/'+username+'/UKCA_TEST_FILES/tebxd/2008jan/'
saving_folder='/nfs/a201/'+username+'/PLOTS_UKCA_KIRSTY/'
saving_folder=files_directory+'PLOTS'
ukl.create_folder(files_directory+'PLOTS')


nc_files=glob(files_directory+'*nc')

file_name=nc_files[5]
file_name=nc_files[0]
#cube=iris.load_cube(file_name)
short_file_name=file_name[len(files_directory):-3]


variable='field38438'
variables_units_name=['t','latitude','longitude','hybrid_ht']
mb=netcdf.netcdf_file(file_name,'r')


for variable in mb.variables:
    print variable
    lat=mb.variables['latitude'].data
    lon=mb.variables['longitude'].data
    pressures=mb.variables['hybrid_ht']
    month=ukl.month_names[ukl.get_months(mb.variables['t'])[0]]
    if mb.variables[variable].data.shape == (1, 85, 145, 192) and variable not in variables_units_name:
        ukl.plot(mb.variables[variable].data[0,0,:,:],
                lat=lat,lon=lon,cblabel=mb.variables[variable].units,
                title=mb.variables[variable].long_name+' '+mb.variables[variable].date+' '+month,
                file_name=saving_folder+short_file_name+mb.variables[variable].long_name+'_'+month,
                saving_format='png',show=0)




#%%


