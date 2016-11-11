# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 12:42:27 2016

@author: earkpr
"""


##import sys
import os
import sys
import getpass
from glob import glob
import logging
import datetime

#%%

########################################################

## Key user defined settings.

jobID = "tebxd"

## Key user defined paths

#dir_scripts='/nfs/a107/earkpr/DataVisualisation/Jesus/git_area/UKCA_postproc-master/'
#orog_file = '/group_workspaces/jasmin2/gassp/jvergaratemprado/n96_hadgem1_qrparm.orog_new.pp'
#jasminorog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'#leeds foe-linux

dir_scripts='/nfs/a107/earkpr/DataVisualisation/Jesus/git_area/UKCA_postproc-master/'
orog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'

## Automatically defined paths

username=getpass.getuser()
LOGFILES_Directory_Path = dir_scripts+'LOGFILES/'

## Location of the Level 0 (raw data files)
input_files_directory='/nfs/a201/'+str(username)+'/UKCA_TEST_FILES/'+str(jobID)+'/'

## Location of where to write Level 1 (processed data files)
output_files_directory='/nfs/a201/'+str(username)+'/UKCA_TEST_FILES/'+str(jobID)+'/'


##  Set path to look in dir_scripts for Python routines
sys.path.append(dir_scripts)

########################################################


#%%   Setting that shouldn't be edited by users very often


################################################

#  Setup the level of logging information required
#  Options are: DEBUG (most info), WARNING, INFO (least info)
#  Default is:  DEBUG

now = datetime.datetime.now()
log = logging.getLogger()
hdlr = logging.FileHandler(LOGFILES_Directory_Path+'/test_log_'+str(jobID)+'_'+str(now).replace(' ', '_')+'.log')

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

##log.flush()
#log.removeHandler(hdlr)
#hdlr.close() 

################################################

