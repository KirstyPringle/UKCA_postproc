#  Routine to write output to a log file.
#  Needs to be run first.
#
#  Kirsty Pringle

import logging
import datetime

now = datetime.datetime.now()

LOGFILES_Directory_Path = '/nfs/see-fs-02_users/earkpr/arch2/DataVisualisation/Jesus/git_area/UKCA_postproc-master/LOGFILES/'

log = logging.getLogger()
hdlr = logging.FileHandler(LOGFILES_Directory_Path+'/log'+str(now).replace(' ', '_')+'.csv')
##FORMAT='%(asctime)s\t%(levelname)s\t%(message)s'
FORMAT='%(levelname)s\t%(message)s'
formatter = logging.Formatter(FORMAT)
logging.basicConfig(format=FORMAT) # log sur console
hdlr.setFormatter(formatter)
log.addHandler(hdlr)
log.setLevel(logging.DEBUG) #set verbosity to show all messages of severity >= DEBUG

log.info('START OF LOGGING')
log.info('================')

##log.flush()
#log.removeHandler(hdlr)
#hdlr.close() 


#%%