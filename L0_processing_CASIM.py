# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""

# KP Moved to UCKA_ControlFile.py dir_scripts='/nfs/see-fs-01_users/earkpr/CODE/UKCA_postproc'#Change this to the downloaded folder

#%%


#import UKCA_ControlFile
#from UKCA_ControlFile import *
#reload(UKCA_ControlFile)

from __main__ import *


# %%
## Done in "UKCA_ControlFile.py"
###sys.path.append(dir_scripts)
import sys
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
import getpass
import multiprocessing
import os
###os.chdir(dir_scripts)
username=getpass.getuser()
iris.FUTURE.netcdf_promote = False
iris.FUTURE.netcdf_no_unlimited =False
import variable_dict as vd
reload(vd)
from multiprocessing import Process, Manager

#%%

#########################################
log.info('START OF LOGGING LEVEL 0 ')
log.info('================')
log.info('Using IRIS version ='+iris.__version__)
#########################################

#KP Moved to UKCA_ControlFile orog_file = '/group_workspaces/jasmin2/gassp/jvergaratemprado/n96_hadgem1_qrparm.orog_new.pp'#jasmin
#KP Moved to UKCA_ControlFile orog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'#leeds foe-linux

#KP Moved to UKCA_ControlFile files_directory='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/'
print input_files_directory+'umnsaa_p*'
pp_files=glob(input_files_directory+'umnsaa_p*')

#########################################
log.info('pp_files = '+str(pp_files))
log.info('================')
#########################################

# pp_files=[pp for pp in pp_files if not 'stash' in pp and not 'xhist' in pp]
manager = Manager()
stashcodes=[]
step_folders=[]
stashcodes = manager.list()
rotate_cube=False
latlon0 = manager.list([0.0,0.0])
step_folders = manager.list()

#%%
def from_pp_to_nc_single_var_single_ts(step_file):
    #print step_file
    # global stashcodes
    # global step_folders
    cubes=iris.load(step_file)#long and heavy bit. Time: around 15 minutes

    #########################################
    log.info("step_file="+str(step_file))
    #########################################

    for cube in cubes:
        #capturing stash code from pp file
        stash_code=ukl.get_stash(cube)
        #print stash_code
        stashcodes.append(stash_code)
        #print stashcodes
        if stash_code in vd.variable_reference_stash:
            if not isinstance(cube.long_name,str):
                cube.long_name=vd.variable_reference_stash[stash_code].long_name
                # print 'added long_name',cube.long_name, 'to', stash_code
                if not isinstance(cube._var_name,str):
                    if not vd.variable_reference_stash[stash_code].short_name=='':
                        cube._var_name=vd.variable_reference_stash[stash_code].short_name
                        # print 'added short_name as cube._var_name',cube._var_name, 'to', stash_code
        times=cube.coord('time').points

        #########################################
        log.info("cube.long_name= "+str(cube.long_name))
        log.info("times= "+str(times))
        #########################################

        if rotate_cube:
            cube.coord('grid_latitude').points=cube.coord('grid_latitude').points+latlon0[0]
            cube.coord('grid_longitude').points=cube.coord('grid_longitude').points+latlon0[1]+180
        for it in range(len(times)):
            cube_single_t=cube.extract(iris.Constraint(time=times[it]))

            folder_NETCDF=output_files_directory+str(int(times[it]))+'/'

            #########################################
            log.info("folder_NETCDF= "+str(folder_NETCDF))
            #########################################

            ukl.create_folder(folder_NETCDF)
            step_folders.append(folder_NETCDF)

            if cube_single_t._standard_name:
                saving_name=folder_NETCDF+str(int(times[it]))+'_'+stash_code+'_'+cube_single_t._standard_name+'.nc'
            elif isinstance(cube.long_name,str):
                saving_name=folder_NETCDF+str(int(times[it]))+'_'+stash_code+'_'+cube_single_t.long_name+'.nc'
            else:
                saving_name=folder_NETCDF+str(int(times[it]))+'_'+stash_code+'.nc'

            iris.save(cube_single_t,saving_name, netcdf_format="NETCDF4")
jobs=[]

processes=20
print 'Number of pp_files for L0', len(pp_files)
list_of_chunks=np.array_split(pp_files,len(pp_files)/processes+1)
start=time.time()
for chunk in list_of_chunks:
    jobs=[]
    for step_file in chunk:
        p = multiprocessing.Process(target=from_pp_to_nc_single_var_single_ts, args=(step_file,))
        print step_file,p
        jobs.append(p)
        p.start()
    for job in jobs:
        job.join()


end=time.time()
print stashcodes
stashcodes=list(set(stashcodes))
step_folders=list(set(step_folders))
step_folders.sort()
print stashcodes

print 'time to convert from pp to single nc:',end-start
#%%
#file_variable_name='2008apr_m01s00i101_mass_fraction_of_sulfur_dioxide_expressed_as_sulfur_in_air.nc'
def join_variables(list_variables):
    for stash_code in list_variables:
        names=[]
        print stash_code, '\n'

        for step_folder in step_folders:
            # print input_files_directory
            #print step_folder
            file_name=ukl.Obtain_name(step_folder,stash_code)
            if len(file_name)>=1:names.append(file_name[0])
        cube_list=[]
        #print len(names)
        #print names
        cube_list=iris.load(names)
        if 'm01s00i033' == stash_code:
            print 'orography skipped'
            continue
        #print cube_list
        cube_list_concatenated=cube_list.concatenate()[0]

        if stash_code in vd.variable_reference_stash:
            if not isinstance(cube_list_concatenated.long_name,str):
                cube_list_concatenated.long_name=vd.variable_reference_stash[stash_code].long_name
                # print 'added long_name',cube_list_concatenated.long_name, 'to', stash_code

        # print cube_list_concatenated.standard_name
        if cube_list_concatenated.standard_name:
            saving_name=folder_all_time_steps+'All_time_steps_'+stash_code+'_'+cube_list_concatenated._standard_name+'.nc'
        elif isinstance(cube_list_concatenated.long_name,str):
            saving_name=folder_all_time_steps+'All_time_steps_'+stash_code+'_'+cube_list_concatenated.long_name+'.nc'
        else:
            saving_name=folder_all_time_steps+'All_time_steps_'+stash_code+'.nc'

        iris.save(cube_list_concatenated,saving_name, netcdf_format="NETCDF4")

        #########################################
        log.info("cube_list_concatenated="+str(cube_list_concatenated.long_name))
        #########################################

#       This part is for calculating mean values across the whole run time
#        cube_run_mean=cube_list_concatenated.collapsed(['time'],iris.analysis.MEAN)
#
#        stash_code=ukl.get_stash(cube_run_mean)
#        if cube_run_mean.long_name:
#            saving_name=folder_run_mean+'run_mean_'+stash_code+'_'+cube_run_mean.long_name+'.nc'
#        elif cube_run_mean._var_name:
#            saving_name=folder_run_mean+'run_mean_'+stash_code+'_'+cube_run_mean._var_name+'.nc'
#        else:
#            saving_name=folder_run_mean+'run_mean_'+stash_code+'.nc'
#        iris.save(cube_run_mean,saving_name, netcdf_format="NETCDF4")

#folder_run_mean=files_directory+'run_mean/'

folder_all_time_steps=input_files_directory+'All_time_steps/'
#ukl.create_folder(folder_run_mean)
ukl.create_folder(folder_all_time_steps)


#########################################
log.info("folder_all_time="+str(folder_all_time_steps))
#########################################


processes=20
print 'Number of variables', len(stashcodes)
list_of_chunks=np.array_split(stashcodes,processes)
jobs=[]
start=time.time()
for chunk in list_of_chunks:
    p = multiprocessing.Process(target=join_variables, args=(chunk.tolist(),))
    print p
    jobs.append(p)
    p.start()

for job in jobs:
    job.join()
end=time.time()
print end-start

############################################
log.info('removing single time step folders')
###########################################


for folder_name in step_folders:
    print 'removing:',folder_name
    os.system('rm -rf %s'%folder_name)







#########################################
log.info('End of L0_processing.py ')
log.info('================')
#########################################
