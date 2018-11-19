# -*- coding: utf-8 -*-
"""
Control file for L0 and L1 processing
"""


import os
import sys
import getpass
from glob import glob
import logging
import datetime

#%%
try:
    sys.path.append('/nfs/a107/eejvt/PYTHON_CODE')
    import Jesuslib as jl
except:
    print 'Jesuslib not accesible'
########################################################


#Set defaults - DO NOT change values in here
execfile("UKCA_Control_defaults.py")

#Override the defaults using the following user script (can change the path to this if needed)
execfile("UKCA_Control_USER_config.py")


## Automatically defined paths
username=getpass.getuser()
LOGFILES_Directory_Path = dir_scripts+'LOGFILES/'

## Location of the pp files (raw data files)
#The location can be given as an argument when you call the script (ipython UKCA_ControlFile.py /nfs/a201/...../) or
#it can be given as:
# input_files_directory='/nfs/a201/'+str(username)+'/'+str(model_name)+'_TEST_FILES/'+str(jobID)+'/'
#or just as:
if len(sys.argv)>1:
    input_files_directory=sys.argv[1]
    if input_files_directory[-1]!='/':
        input_files_directory=input_files_directory+'/'
    try:
        run_L0=int(sys.argv[2])
        run_L1=int(sys.argv[3])
        run_plots=int(sys.argv[4])
        print 'run booleans writen from command line'
        print run_L0,run_L1,run_plots
    except:
        asdfasdf=3452

## Location of where to write Level 0 (data files in nc format). Will be the same as input_files_directory unless the flag below is True
if l_set_output_dir:
    output_files_directory=output_files_directory_set
else:
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
hdlr = logging.FileHandler(LOGFILES_Directory_Path+'/Log_'+str(jobID)+'_'+sys.argv[0].split('/')[-1][:-3]+'_'+str(now).replace(' ', '_')+'.log')

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
    elif model_name == "CASIM":
        execfile("L0_processing_CASIM.py")
    else:
        execfile("L0_processing.py")

if run_L1:
    if model_name == "CASIM":
        execfile("L1_processing_CASIM.py")
    else:
        execfile("L1_processing.py")
if run_plots:
    execfile("Plots_for_netCDF4.py")
        
#if send_mail:
#    jl.send_email()


##log.flush()
#log.removeHandler(hdlr)
#hdlr.close()

################################################
