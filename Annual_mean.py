
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
import os
import UKCA_lib as ukl
import getpass
iris.FUTURE.netcdf_promote = True
username=getpass.getuser()

files_directory_UKCA='/nfs/a201/'+username+'/UKCA_TEST_FILES/'
run='tebxd/'
files_directory=files_directory_UKCA+run

folder_annual_mean=files_directory+'Annual_mean/'
ukl.create_folder(folder_annual_mean)
year='2008'
sample_month_folder=year+'jan/'

list_variable_names=glob(files_directory+sample_month_folder+'*.nc')
for i in range(len(list_variable_names)):
    list_variable_names[i]=list_variable_names[i][len(files_directory+sample_month_folder):]
    print list_variable_names,'\n'
#Loop over all variables and save in annual_mean folder

#file_variable_name='2008apr_m01s00i101_mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air.nc'
for file_variable_name in list_variable_names:
    names=[]
    print file_variable_name, '\n'
    for imon in range(12):
        print files_directory
        name=files_directory+'2008'+ukl.months_str[imon]+'/'+file_variable_name[:4]+ukl.months_str[imon]+file_variable_name[7:]
        print name, '\n'
        names.append(name)
    cube_list=iris.load(names)

    cube_list_concatenated=cube_list.concatenate_cube()
    cube_annual_mean=cube_list_concatenated.collapsed(['time'],iris.analysis.MEAN)

    stash_code=ukl.get_stash(cube_annual_mean)
    if isinstance(cube_annual_mean.long_name,str):
        saving_name=folder_NETCDF+date+'_'+stash_code+'_'+cube_annual_mean.long_name+'.nc'
    else:
        saving_name=folder_annual_mean+'Annual_mean'+'_'+stash_code+'.nc'

    iris.save(cube_annual_mean,saving_name, netcdf_format="NETCDF4")



#plt.plot([cube.data.mean() for cube in cube_list])
