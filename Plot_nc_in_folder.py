# -*- coding: utf-8 -*-
"""
import UKCA_lib as ukl
reload(ukl)
import netCDF4
import matplotlib
import getpass


username=getpass.getuser()

plt.interactive(0)

folder='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/L1/'


saving_folder=folder+'PLOTS/'
ukl.create_folder(saving_folder)
nc_files=glob(folder+'*nc')
for nc_file in nc_files:
    #mb=netCDF4.Dataset(nc_file,'r')
    cube=iris.load(nc_file)[0]
    cube=cube.collapsed(['time'],iris.analysis.MEAN)
    print cube._var_name
    print cube.shape
    #print cube.shape
    if cube.ndim==2 and cube.shape==(145, 192):
        qplt.contourf(cube,cmap=plt.cm.RdBu_r)
        plt.gca().coastlines()
        if isinstance(cube.long_name,str):
            name=cube.long_name
        elif isinstance(cube._var_name,str):
            name=cube._var_name
        else:
            name=nc_file[len(folder):-3]

        saving_str=saving_folder+'2D_plot'+name+'.png'
        plt.savefig(saving_str)
        plt.close()
    if cube.ndim==3 and cube.shape==(85, 145, 192):
        if isinstance(cube.long_name,str):
            name=cube.long_name
        elif isinstance(cube._var_name,str):
            name=cube._var_name
        else:
            name=nc_file[len(folder):-3]
        name_log='Logscale_'+name
        ukl.zonal_mean_plot(cube,saving_folder,name_log,logscale=1)
        ukl.zonal_mean_plot(cube,saving_folder,name)
        ukl.level_plot(cube,saving_folder,name_log,logscale=1)
        ukl.level_plot(cube,saving_folder,name)
        ukl.level_plot(cube,saving_folder,name_log,level=22,logscale=1)#around 600hpa
        ukl.level_plot(cube,saving_folder,name,level=22)#around 600hpa
