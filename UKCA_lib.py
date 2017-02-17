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
#sys.path.append('/nfs/a107/eejvt/PYTHON_CODE')
#
from glob import glob
import matplotlib as mpl
#mpl.use('Agg')
from glob import glob
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from matplotlib import colors, ticker, cm
import matplotlib
import matplotlib.pyplot as plt
import iris.plot as iplt
from scipy.io import netcdf
import iris.quickplot as qplt
import datetime
import os
import scipy as sp


test_run_path='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/'
months_str=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
months_str_upper_case=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
month_names=['January','February','March','April','May','June','July','August','September','October','November','December']


def plot(data,title=' ',projection='cyl',file_name=datetime.datetime.now().isoformat(),show=1,cblabel='$\mu g/ m^3$',cmap=plt.cm.CMRmap_r,clevs=np.zeros(1),return_fig=0,dpi=300,lon=0,lat=0,colorbar_format_sci=0,saving_format='svg',scatter_points=0,f_size=20):
    fig=plt.figure(figsize=(20, 12))
    m = fig.add_subplot(1,1,1)
    if projection=='merc':
        m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20)
    else:
        m = Basemap(projection=projection,lon_0=0)
    m.drawcoastlines()

    if isinstance(lon, int):
        lon=readsav('/nfs/a107/eejvt/IDL_CODE/glon.sav')
    if isinstance(lat, int):
        lat=readsav('/nfs/a107/eejvt/IDL_CODE/glat.sav')

        X,Y=np.meshgrid(lon.glon,lat.glat)
    else:
        if lon.ndim==1:
            X,Y=np.meshgrid(lon,lat)
        else:
            X=np.copy(lon)
            Y=np.copy(lat)
    if type(clevs) is list:
        cs=m.contourf(X,Y,data,clevs,latlon=True,cmap=cmap,norm= colors.BoundaryNorm(clevs, 256))
        if colorbar_format_sci:
            def fmt(x, pos):
                a, b = '{:.1e}'.format(x).split('e')
                b = int(b)
                return r'${} \times 10^{{{}}}$'.format(a, b)
            cb = m.colorbar(cs,"right",format=ticker.FuncFormatter(fmt),ticks=clevs)
        else:
            cb = m.colorbar(cs,format='%.2e',ticks=clevs)
    else:
        cs=m.contourf(X,Y,data,15,latlon=True,cmap=cmap)
        cb = m.colorbar(cs)

    if not isinstance(scatter_points,int):
        m.scatter(scatter_points[:,0],scatter_points[:,100])
    cb.set_label(cblabel,fontsize=f_size)
    cb.ax.tick_params(labelsize=f_size)
    plt.title(title,fontsize=f_size)
    if os.path.isdir("PLOTS/"):
        plt.savefig('PLOTS/'+file_name+'.'+saving_format,format=saving_format,dpi=dpi, bbox_inches='tight')
        plt.savefig('PLOTS/'+file_name+'.svg',format='svg', bbox_inches='tight')
    else:
        plt.savefig(file_name+'.'+saving_format,format=saving_format,dpi=dpi, bbox_inches='tight')
        plt.savefig(file_name+'.svg',format='svg', bbox_inches='tight')
    if show:
        plt.show()
    if return_fig:
        return fig


def zonal_mean_plot(cube,saving_path,name,cmap='CMRmap_r',logscale=0):
    for coord in cube.coords():
        #print coord.name()
        if coord.name()=='surface_altitude':
            cube.remove_coord('surface_altitude')
            #print 'removed'
    #coords_name=[coord.name() for coord in cube.coords()]
    #print coords_name
    cube_zonal_mean=cube.collapsed(['longitude'],iris.analysis.MEAN)
    if logscale:
        qplt.contourf(cube_zonal_mean,cmap=cmap,norm=matplotlib.colors.LogNorm())
    else:
        qplt.contourf(cube_zonal_mean,cmap=cmap)
    # data=cube_zonal_mean.data
    #
    # title=cube.var_name+' max:%1.2f mean:%1.2f min:%1.2f'%(data.max(),data.mean(),data.min())
    # plt.title(title)
    plt.yscale('log')
    plt.savefig(saving_path+'Zonal_mean_'+name+'.png',bbox_inches='tight')
    plt.close()


def print_cube_single_value(cube):
    """
    Prints the value of the first element of a cube.
    Flexible to deal with different cube sizes

    This function is required as sometimes an element of the cube needs to be
    printed out in order to force the cube to be calculated.
    """
    val=cube.shape
    ind=[0 for v in val]
    print eval("cube.data"+str(ind))




def level_plot(cube,saving_path,name='',level=0,color_levels=9,cmap=plt.cm.CMRmap_r,logscale=0,saving_format='.png'):

    '''
    This function works for 3 dimensional cubes (model_level_number, latitude, longitude)


    It plots and saves a png file (by default)
    You can use it like:
        ukl.level_plot(cube_time_mean,saving_path)

    By default, it plots the cube at level 0 (surface_level) in linear scale and saves it in the path given.

    you can change 'level' for plotting a different level
    For example

    lev=22
    ukl.level_plot(cube_time_mean,saving_path,level=lev)


    Other kargs:

    'name' sets a different name in the saved file. By default it uses cube.var_name
    'color_levels' is an integrer number for setting how many levels you want
    'logscale' if set to true, the plot will be in logarithmic scale
    'cmap' changes the mapping colors
    'saving_format' can be set to something different than png to change the format of the plot

    '''

    if cube.ndim!=3:
        raise NameError('The cube has to have 3 dimensions (model_level_number, latitude, longitude) \n \
        Currently its shape is: %s' % (cube.shape,) )

    if logscale:
        qplt.contourf(cube[level,],color_levels,cmap=cmap,norm=matplotlib.colors.LogNorm())
        log_str='_log_scale'
    else:
        qplt.contourf(cube[level,],color_levels,cmap=cmap)
        log_str=''

    plt.gca().coastlines()
    if name=='':
        name=cube.var_name
    if level==0:
        saving_str=saving_path+'Surface_level_'+name+log_str+saving_format
    else:
        saving_str=saving_path+'Level_%i_'%level+name+log_str+saving_format
    plt.savefig(saving_str,bbox_inches='tight')
    plt.close()








def lognormal_PDF(rmean,r_list,std):
   X=(1/(r_list*np.log(std)*np.sqrt(2*np.pi)))*np.exp(-(np.log(r_list)-np.log(rmean))**2/(2*np.log(std)**2))
   return X
def lognormal_cummulative(N,r,rbar,sigma):
    total=(N/2)*(1+sp.special.erf(np.log(r/rbar)/np.sqrt(2)/np.log(sigma)))
    return total
def lognormal_cummulative_forcubes(N,r,rbar,sigma):
    total=(N/2)*(1+sp.special.erf(iris.analysis.maths.log(r/rbar)/np.sqrt(2)/np.log(sigma)))
    return total

class log_steps():
    def __init__(self,start,final,points=10000):
        #note that start and final have to be given in logaritmic units
        self.step_limits=np.logspace(start,final,points)
        self.step_size=self.step_limits[1:]-self.step_limits[:-1]
        self.mid_points=self.step_limits[:-1]+self.step_size[:]/2




def Obtain_name_from_list(str_list,string,end=0,just_one=0):
    names=[name for name in str_list if string in name]
    return names

def Obtain_name(folder,string,end=0,just_one=0):
    if folder[-1]!='/':
        folder=folder+'/'
    files_list=glob(folder+'*')
    names=Obtain_name_from_list(files_list,string)
    if end:
        names=[name for name in names if end == name[:-len(end)]]
    if just_one and len(names)>1:
        raise NameError('more than one value of %s in %folder and just one needed'%(string,folder))
    if just_one and len(names)==1:
        raise NameError('not file founded with %s in %s'%(string,folder))
    if just_one:
        return names[0]#to return a single value instead of a list
    else:
        return names

#%%
class VariableAttributes:
    """
     Class = VariableAttributes
    """
    def __init__(self,stash_code,name,short_name,long_name,units,description='None'):
        self.stash_code=stash_code
        self.name=name
        self.short_name=short_name
        self.long_name=long_name
        self.units=units
        self.description=description
#%%
class SpeciesAttributes:
    def __init__(self,name,mm,rhocomp,kappa,description='None'):
        self.name=name
        self.mm=mm
        self.rhocomp=rhocomp
        self.kappa=kappa
        self.description=description


class ModalAttributes:
    def __init__(self,name,sigma,ddplim0,ddplim1,modesol,mode_choice,description='None'):
        self.name=name
        self.sigma=sigma
        self.ddplim0=ddplim0
        self.ddplim1=ddplim1
        self.modesol=modesol
        self.mode_choice=mode_choice
        self.description=description
    def plot_mode(self,r_mean,N=1,real_PDF=False,limits=False):
        '''
        The PDF returned is weigthed by the step size of the radius in order to obtain a nice ilustrative plot
        for scientific purposes, note that the function shown is not a real PDF unless real_PDF is set to True
        '''
        rs=log_steps(-10,-4,1000)
        if self.modesol==0:
            pl_ls='-.'
        else:
            pl_ls='-'#*rs.step_size
        if real_PDF:
            plt.plot(rs.mid_points,N*lognormal_PDF(r_mean, rs.mid_points,self.sigma),label=self.name,ls=pl_ls)
        else:
            plt.plot(rs.mid_points,N*lognormal_PDF(r_mean, rs.mid_points,self.sigma)*rs.step_size,label=self.name,ls=pl_ls)

        plt.axvline(r_mean,ls='-',c='k')
        if limits:
            plt.axvline(self.ddplim0,ls='--',c='k')
            plt.axvline(self.ddplim1,ls='--',c='k')
        plt.xscale('log')
        if real_PDF:
            plt.yscale('log')
            plt.ylim(1e2,1e10)
            plt.ylabel('$PDF$ (Probability Density Function)')
        else:
            plt.ylabel('Arbitrary units')
        plt.xlabel('Radius $(m^{-3})$')
        return rs

def get_months(time,date_list=[2006,12,1]):
    st_year=date_list[0]
    st_month=date_list[1]
    st_day=date_list[2]
    t=datetime.datetime.fromtimestamp(0)
    diff_sec=(t-datetime.datetime(st_year,st_month,st_day)).total_seconds()
    time_arr=time[:]-diff_sec
    if isinstance(time_arr,float):
        time_arr_struct=datetime.datetime.fromtimestamp(time_arr)
        months=time_arr_struct.month
    else:
        convert_time_vectorized=np.vectorize(datetime.datetime.fromtimestamp)
        time_arr_struct=convert_time_vectorized(time_arr[:,])
        months=[time_arr_struct[i].month for i in range(len(time_arr_struct))]
    return months

def normalize_date(time,date_list=[2006,12,1]):
    st_year=date_list[0]
    st_month=date_list[1]
    st_day=date_list[2]
    t=datetime.datetime.fromtimestamp(0)
    diff_sec=(t-datetime.datetime(st_year,st_month,st_day)).total_seconds()
    time_arr=time[:]-diff_sec
    if isinstance(time_arr,float):
        time_arr_struct=datetime.datetime.fromtimestamp(time_arr)
    else:
        convert_time_vectorized=np.vectorize(datetime.datetime.fromtimestamp)
        time_arr_struct=convert_time_vectorized(time_arr[:,])
    return time_arr_struct


def log_levels(data_map,levels_per_order=2):
    data_map=data_map[np.logical_not(np.isnan(data_map))]
    data_map=data_map[np.logical_not(data_map==0)]
    maxmap=data_map.max()
    minmap=data_map.min()
    lim_max=int(1000+np.log10(maxmap))-1000+1
    lim_min=int(1000+np.log10(minmap))-1000
    orders_of_magnitude=lim_max-lim_min
    levels=np.logspace(lim_min,lim_max,levels_per_order*orders_of_magnitude+1)
    return levels.tolist()


def create_folder(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def beauty_colorscale(data,levels=9):
    len_one_dim=data.shape[0]*data.shape[1]
    one_dim_array=np.reshape(data,len_one_dim)
    levels_list=[]
    #levels_list.append(one_dim_array[0])
    for i in range (levels+1):
        #if levels_list[i]==one_dim_array[i/levels*(len_one_dim-1)]:

        levels_list.append(one_dim_array[i/levels*(len_one_dim-1)])
    return levels_list

from cStringIO import StringIO


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout




def get_stash(cube):
   with Capturing() as output:
       print cube.attributes['STASH']
   stash_code=output[0]
   return stash_code

def get_stash_from_numbers(model,section,item):
   """
     Gets STASH code from the individual section and item numbers.

     eg converts ITEM=02, SEC=123 to STASH = 02123
   """
   with Capturing() as output:
       print iris.fileformats.pp.STASH(model,section,item)
   stash_code=output[0]
   return stash_code


def extract_stcodes(stcode_file):
   """
     Extract Stash codes from text file.
     Can handle space, line,comma seperated or a combination of these

     Returns a integer list of STASH codes
   """
   try:
     f = open(stcode_file, 'r')
   except:
     raise PncError(" ERROR: Unable to read stashcodes file ")

   comma = ','
   stcodes = [ ]
   for stline in f:
     stline = stline.rsplit('\n')[0]
     if len(stline) > 5:       # multiple codes on a line
       if comma in stline:       # comma separated
         stc = stline.split(',')
       else:                   # space separated
         stc = stline.split()

       for s in stc:
         stcodes.append( int(s) )
     else:
       stcodes.append( int(stline) )

     print stcodes
   return stcodes
# End def extract_stcodes

def unrotated_grid(cube):
    rotated_cube=isinstance(cube.coord('grid_longitude').coord_system,iris.coord_systems.RotatedGeogCS)
    if rotated_cube:
        pole_lat=cube.coord('grid_longitude').coord_system.grid_north_pole_latitude
        pole_lon=cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
        lons, lats =iris.analysis.cartography.unrotate_pole(cube.coord('grid_longitude').points,cube.coord('grid_latitude').points,pole_lon,pole_lat)
    else:
        lons=cube.coord('grid_longitude').points
        lats=cube.coord('grid_latitude').points
    return lons,lats
