# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 11:06:34 2016

@author: eejvt
"""

import os
import sys
import getpass
from glob import glob
import logging
import datetime
import numpy as np


log=np.genfromtxt('Mapping_ensable_numbers.csv',dtype=str,delimiter=',')

names=log[:,0].tolist()
numbers=map(int,log[:,1])




now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
log = logging.getLogger()
hdlr = logging.FileHandler('Ensemble_log_'+sys.argv[0].split('/')[-1][:-3]+'_'+str(now).replace(' ', '_')+'.log')

FORMAT='%(levelname)s\t%(message)s'

formatter = logging.Formatter(FORMAT)
logging.basicConfig(format=FORMAT) # log sur console
hdlr.setFormatter(formatter)
log.addHandler(hdlr)

log.setLevel(logging.DEBUG) 
username=getpass.getuser()

log.info('START OF LOGGING')
log.info('================')
log.info('username = '+str(username))


failed=[]
analysed=[]
for i in range(len(names)):
    print 'Number:   ',i    
    run_name=names[i]
    run_number=numbers[i]
    out=os.system('python UKCA_ControlFile.py %s %i'%(run_name,run_number))        
    
    if out==0:
        log.info('Run: %s Number: %i succesfully analysed'%(run_name,run_number))
        analysed.append(run_name)
    else:
        failed.append(run_name)
        log.info('#########################ERROR############################')
        log.info('Run: %s Number: %i Failed'%(run_name,run_number))
        log.info('#########################ERROR############################')
        
np.savetxt('Failed_runs',failed)
np.savetxt('Succesfull_runs',analysed)

