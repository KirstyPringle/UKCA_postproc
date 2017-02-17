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


from __main__ import *


sys.path.append(dir_scripts)



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

folders=[files_directory+'All_time_steps/',files_directory+'L1/']



for data_folder in folders:
    saving_folder=data_folder+'PLOTS/'
    ukl.create_folder(saving_folder)
    nc_files=glob(data_folder+'*nc')
    for nc_file in nc_files:
        #mb=netCDF4.Dataset(nc_file,'r')
        cube=iris.load(nc_file)[0]
        cube=cube.collapsed(['time'],iris.analysis.MEAN)
        print cube.var_name
        print cube.shape
        try:
            if cube.ndim==2:
                qplt.contourf(cube,cmap=plt.cm.RdBu_r)
                plt.gca().coastlines()
                if isinstance(cube.var_name,str) or isinstance(cube.var_name,unicode):
                    name=ukl.get_stash(cube)+'_'+cube.var_name
                else:
                    name=ukl.get_stash(cube)
                saving_str=saving_folder+'2D_plot'+name+'.png'
                plt.savefig(saving_str)
                plt.close()
            if cube.ndim==3:
                if isinstance(cube.var_name,str) or isinstance(cube.var_name,unicode):
                    name=ukl.get_stash(cube)+'_'+cube.var_name
                else:
                    name=ukl.get_stash(cube)
                ukl.zonal_mean_plot(cube,saving_folder,name)
                ukl.zonal_mean_plot(cube,saving_folder,name,logscale=True)
                ukl.level_plot(cube,saving_folder,name)
                ukl.level_plot(cube,saving_folder,name,logscale=True)
                ukl.level_plot(cube,saving_folder,name,logscale=True,level=22)#around 600hpa
        except:
            print cube
            print 'not managed to plot automatically. Check wether if there is any other dimension apart from [time, level, lat, lon]'






#
