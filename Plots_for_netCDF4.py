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
####sys.path.append('/nfs/see-fs-01_users/eejvt/UKCA_visualization')
 
from glob import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
from scipy.io import netcdf
import datetime
import UKCA_lib as ukl
import netCDF4
import matplotlib
import getpass


username=getpass.getuser()

plt.interactive(0)
files_directory='/nfs/a201/'+username+'/UKCA_TEST_FILES/tebxd/'
#saving_folder=files_directory+'PLOTS'

year='2008'
folders=[]
for imon in ukl.months_str:
    folders.append(year+imon+'/')
if os.path.isdir(files_directory+path):
    folders.append('Annual_mean')




for month_folder in folders:
    saving_folder=files_directory+month_folder+'PLOTS/'
    ukl.create_folder(saving_folder)
    nc_files=glob(files_directory+month_folder+'*nc')
    for nc_file in nc_files:
        #mb=netCDF4.Dataset(nc_file,'r')
        cube=iris.load(nc_file)[0]
        print cube.shape
        if cube.ndim==2 and cube.shape==(145, 192):
            qplt.contourf(cube,cmap=plt.cm.RdBu_r)
            plt.gca().coastlines()
            if isinstance(cube.long_name,str):
                name=ukl.get_stash(cube)+'_'+cube.long_name
            else:
                name=ukl.get_stash(cube)
            saving_str=saving_folder+'2D_plot'+name+'.png'
            plt.savefig(saving_str)
            plt.close()
        if cube.ndim==3 and cube.shape==(85, 145, 192):
            print 1
            if isinstance(cube.long_name,str):
                name=ukl.get_stash(cube)+'_'+cube.long_name
            else:
                name=ukl.get_stash(cube)
            ukl.zonal_mean_plot(cube,saving_folder,name)
            ukl.level_plot(cube,saving_folder,name)
            ukl.level_plot(cube,saving_folder,name,level=22)#around 600hpa
#%%
#cube.aux_coords[4][22]
#%%
'''
nc_file='/nfs/a201/'+username+'/UKCA_TEST_FILES/tebxd/2008jan/2008jan_m01s34i079_mass_fraction_of_hydrogen_sulfide_in_air.nc'

a=ukl.Obtain_name(files_directory+folders[0],'mass')
nc_file=a[8]
cube=iris.load(nc_file)[0]
logscale=1
level=0
data=cube.data[level,]

'''

#%%



'''
nc_file=nc_files[5]
file_name=nc_files[14]
#cube=iris.load_cube(file_name)
short_file_name=file_name[len(files_directory):-3]


variable='field38438'
variables_units_name=['t','latitude','longitude','hybrid_ht']
mb=netCDF4.Dataset(file_name,'r')


#%%
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
'''
