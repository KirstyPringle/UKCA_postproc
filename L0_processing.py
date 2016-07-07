# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""

import sys
dir_scripts='/nfs/see-fs-01_users/eejvt/UKCA_postproc'#Change this to the downloaded folder
sys.path.append(dir_scripts)
import UKCA_lib as ukl
import numpy as np
import time
import iris
from glob import glob
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import datetime
from scipy.io import netcdf
import os
import getpass
import multiprocessing
os.chdir(dir_scripts)
username=getpass.getuser()
iris.FUTURE.netcdf_promote = False
import variable_dict as vd
reload(vd)



files_directory='/nfs/a201/eejvt/CASSIM/first_runs/SO/673/'
pp_files=glob(files_directory+'umnsaa_*')
# pp_files=[pp for pp in pp_files if not 'stash' in pp and not 'xhist' in pp]


def from_pp_to_nc_single_var_single_ts(step_file):
    print step_file
    cubes=iris.load(step_file)#long and heavy bit. Time: around 15 minutes
    for cube in cubes:
        #capturing stash code from pp file
        stash_code=ukl.get_stash(cube)
        if stash_code in vd.variable_reference_stash:
            if not isinstance(cube.long_name,str):
                cube.long_name=vd.variable_reference_stash[stash_code].long_name
                print 'added long_name',cube.long_name, 'to', stash_code
                if not isinstance(cube._var_name,str):
                    cube._var_name=vd.variable_reference_stash[stash_code].short_name
                    print 'added short_name as cube._var_name',cube._var_name, 'to', stash_code
        times=cube.coord('time').points
        for it in range(len(times)):
            cube_single_t=cube.extract(iris.Constraint(time=times[it]))

            folder_NETCDF=files_directory+str(int(times[it]))+'/'
            ukl.create_folder(folder_NETCDF)

            if cube._var_name:
                saving_name=folder_NETCDF+str(int(times[it]))+'_'+stash_code+'_'+cube._var_name+'.nc'
            else:
                saving_name=folder_NETCDF+str(int(times[it]))+'_'+stash_code+'.nc'

            iris.save(cube,saving_name, netcdf_format="NETCDF4")


jobs=[]
start=time.time()
for step_file in pp_files:
    p = multiprocessing.Process(target=from_pp_to_nc_single_var_single_ts, args=(step_file,))
    jobs.append(p)
    p.start()

for job in jobs:
    job.join()

end=time.time()

print end-start
'''
#file_variable_name='2008apr_m01s00i101_mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air.nc'
def join_variables(list_variables):
    for file_variable_name in list_variables:
        names=[]
        print file_variable_name, '\n'
        for imon in range(12):
            print files_directory
            name=files_directory+year+ukl.months_str[imon]+'/'+file_variable_name[:4]+ukl.months_str[imon]+file_variable_name[7:]
            print name, '\n'
            names.append(name)
        cube_list=[]
        cube_list=iris.load(names)
        try:
            print cube_list[0].long_name
        except:
            jfskjsf=1
        if 'm01s00i033' in file_variable_name:
            print 'orography skipped'
            continue
        indx=0
        if len(cube_list)>1:
            if cube_list[0]==cube_list[1]:
                indx=-1
        cube_list_concatenated=cube_list[indx]
        stash_code=ukl.get_stash(cube_list_concatenated)
        if cube_list_concatenated.long_name:
            saving_name=folder_all_months+'All_months_'+stash_code+'_'+cube_list_concatenated.long_name+'.nc'
        elif cube_list_concatenated._var_name:
            saving_name=folder_all_months+'All_months_'+stash_code+'_'+cube_list_concatenated._var_name+'.nc'
        else:
            saving_name=folder_all_months+'All_months_'+stash_code+'.nc'

        iris.save(cube_list_concatenated,saving_name, netcdf_format="NETCDF4")

        cube_annual_mean=cube_list_concatenated.collapsed(['time'],iris.analysis.MEAN)

        stash_code=ukl.get_stash(cube_annual_mean)
        if cube_annual_mean.long_name:
            saving_name=folder_annual_mean+'Annual_mean_'+stash_code+'_'+cube_annual_mean.long_name+'.nc'
        elif cube_annual_mean._var_name:
            saving_name=folder_annual_mean+'Annual_mean_'+stash_code+'_'+cube_annual_mean._var_name+'.nc'
        else:
            saving_name=folder_annual_mean+'Annual_mean_'+stash_code+'.nc'
        iris.save(cube_annual_mean,saving_name, netcdf_format="NETCDF4")


folder_annual_mean=files_directory+'Annual_mean/'

folder_all_months=files_directory+'All_months/'
ukl.create_folder(folder_annual_mean)
ukl.create_folder(folder_all_months)

sample_month_folder=year+'jan/'
print p
list_variable_names=glob(files_directory+sample_month_folder+'*.nc')
#cut variable names
for i in range(len(list_variable_names)):
    list_variable_names[i]=list_variable_names[i][len(files_directory+sample_month_folder):]
    print list_variable_names,'\n'
#Loop over all variables and save in annual_mean folder

processes=20
list_of_chunks=np.array_split(list_variable_names,processes)
jobs=[]
start=time.time()
for chunk in list_of_chunks:
    p = multiprocessing.Process(target=join_variables, args=(chunk.tolist(),))
    jobs.append(p)
    p.start()

for job in jobs:
    job.join()
end=time.time()
print end-start
'''
