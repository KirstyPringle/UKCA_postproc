# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:42:27 2016

@author: earkpr
"""


import os
import sys
import getpass
from glob import glob
import logging
import datetime

#%%

########################################################

## Key user defined settings.

####IMPORTANT####
dir_scripts='/nfs/see-fs-01_users/eejvt/UKCA_final/'
# This is the directory where the code is. Change it to the path where you downloaded the file
####IMPORTANT####

model_name = "UKCA"
#model_name = "TOMCAT_GLOMAP"
jobID = "tebxe"
#jobID = "glo301"

run_L0=True
run_L1=True
run_plots=True



## Key user defined paths

#dir_scripts='/nfs/a107/earkpr/DataVisualisation/Jesus/git_area/UKCA_postproc-master/'
#orog_file = '/group_workspaces/jasmin2/gassp/jvergaratemprado/n96_hadgem1_qrparm.orog_new.pp'
#jasminorog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'#leeds foe-linux

## Automatically defined paths

username=getpass.getuser()
LOGFILES_Directory_Path = dir_scripts+'LOGFILES/'

## Location of the pp files (raw data files)
#it can be given as:
input_files_directory='/nfs/a201/'+str(username)+'/'+str(model_name)+'_TEST_FILES/'+str(jobID)+'/'
#or just as:
input_files_directory='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/'
## Location of where to write Level 0 (data files in nc format). Typically will be the same as input_files_directory.
# output_files_directory='/nfs/a201/'+str(username)+'/'+str(model_name)+'_TEST_FILES/'+str(jobID)+'/'
output_files_directory=input_files_directory

##  Set path to look in dir_scripts for Python routines
sys.path.append(dir_scripts)

########################################################


#%%   Setting that shouldn't be edited by users very often


################################################

# For UKCA runs, need to also define an orography file.
if model_name is "UKCA":
    orog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'


#  Setup the level of logging information required
#  Options are: DEBUG (most info), WARNING, INFO (least info)
#  Default is:  DEBUG

if not os.path.exists(LOGFILES_Directory_Path):
        os.makedirs(LOGFILES_Directory_Path)

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
log = logging.getLogger()
hdlr = logging.FileHandler(LOGFILES_Directory_Path+'/test_log_'+str(jobID)+'_'+sys.argv[0].split('/')[-1][:-3]+'_'+str(now).replace(' ', '_')+'.log')

##FORMAT='%(asctime)s\t%(levelname)s\t%(message)s'
FORMAT='%(levelname)s\t%(message)s'

formatter = logging.Formatter(FORMAT)
logging.basicConfig(format=FORMAT) # log sur console
hdlr.setFormatter(formatter)
log.addHandler(hdlr)

log.setLevel(logging.DEBUG) #set verbosity to show all messages of severity >= DEBUG
#log.setLevel(logging.WARNING) #set verbosity to show all messages of severity >= DEBUG
#log.setLevel(logging.INFO) #set verbosity to show all messages of severity >= DEBUG

log.info('START OF LOGGING')
log.info('================')
log.info('username = '+str(username))
log.info("jobID = "+str(jobID))
log.info("input_files_directory = "+str(input_files_directory))
log.info("output_files_directory = "+str(output_files_directory))


if run_L0:
    if model_name == "TOMCAT_GLOMAP":
        execfile("L0_processing_TOMCAT_GLOMAP.py")
    else:
        execfile("L0_processing.py")
if run_L1:
    execfile("L1_processing.py")
if run_plots:
    execfile("Plots_for_netCDF4")
        



##log.flush()
#log.removeHandler(hdlr)
#hdlr.close()

################################################
