# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""
#
# import UKCA_ControlFile
# from UKCA_ControlFile import *
# reload(UKCA_ControlFile)
from __main__ import *


# import sys
### KP Moved to UKCA_ControlFile  dir_scripts='/nfs/see-fs-01_users/eejvt/CODE/UKCA_postproc/'#Change this to the downloaded folder
# sys.path.append(dir_scripts)
import numpy as np
import iris
import UKCA_lib as ukl
reload(ukl)
import multiprocessing
import variable_dict as vd
reload(vd)
import I_MODE_SETUP_Variables as ims
reload(ims)

print 'Starting L1 processing'
def save_cube(cube):
    """
    Saves cube as a netCDF file.
    """
    saving_name=saving_folder_l1+'L1_'+cube._var_name+'_'+cube.long_name+'.nc'
    iris.save(cube,saving_name, netcdf_format="NETCDF4")
    print 'saved:',cube.long_name




####files_directory='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/'
folder=output_files_directory+'All_time_steps/'
saving_folder_l1=output_files_directory+'L1/'
ukl.create_folder(saving_folder_l1)
print folder

#Reading necesary cubes
potential_temperature=iris.load(ukl.Obtain_name(folder,'m01s00i004'))[0]
air_pressure=iris.load(ukl.Obtain_name(folder,'m01s00i408'))[0]

p0 = iris.coords.AuxCoord(1000.0,
                          long_name='reference_pressure',
                          units='hPa')
p0.convert_units(air_pressure.units)

Rd=287.05 # J/kg/K
cp=1005.46 # J/kg/K
Rd_cp=Rd/cp

temperature=potential_temperature*(air_pressure/p0)**(Rd_cp)


print temperature.data[0,0,0,0]
temperature.long_name='Temperature'
temperature._var_name='temperature'

save_cube(temperature)
R_specific=iris.coords.AuxCoord(287.058,
                          long_name='R_specific',
                          units='J-kilogram^-1-kelvin^-1')#J/(kg·K)

air_density=(air_pressure/(temperature*R_specific))
print air_density.data[0,0,0,0]
molar_mass_air=iris.coords.AuxCoord(28.991e-3,
                          long_name='Molar mass of air',
                          units='kilogram-mole^-1')#J/(kg·K)
avogadro_number=iris.coords.AuxCoord(6.022e23,
                          long_name='Avogadros number - particles per mol',
                          units='mole^-1')#J/(kg·K)

particle_density_of_air=air_density/molar_mass_air*avogadro_number

print particle_density_of_air.data[0,0,0,0]
air_density._var_name='air_density'
air_density.long_name='Density of air'
save_cube(air_density)

try:
    cloud_number=iris.load(ukl.Obtain_name(folder,'m01s00i075_CLOUD_NUMBER_AFTER_TIMESTEP'))[0]

    CDNC=cloud_number*air_density
    CDNC.long_name='Cloud droplet number concentration'
    CDNC._var_name='CDNC'
    CDNC.units='meter^-3'

    save_cube(CDNC)

    cloud_ice_number=iris.load(ukl.Obtain_name(folder,'m01s00i271_CLOUD_ICE_(CRYSTALS)_AFTER_TIMESTEP'))[0]

    CINC=cloud_ice_number*air_density
    CINC.long_name='Cloud ice number concentration'
    CINC._var_name='CINC'
    CINC.units='meter^-3'
    save_cube(CINC)
except:
    print 'could not calculate ice and water number'
    
import copy
#m01s15i101_height_above_reference_ellipsoid.nc
height=iris.load(ukl.Obtain_name(folder,'m01s15i101'))[0]
base=np.zeros(height.data.shape[1:])
length_gridbox_cube=height.copy()#copy.deepcopy(height)
length_gridbox=np.zeros(height.data.shape)
#length_gridbox.data=np.zeros(length_gridbox.data.shape)
#%%
for i in range(height.data.shape[0]):
    if i==0:
        length_gridbox[0,].data=height[0,].data
    else:
        length_gridbox[i,]=height[i,].data-height[i-1,].data
        print height[i,0,0].data
        print height[i-1,0,0].data
        print '---'
        print height[i,0,0].data-height[i-1,0,0].data
        print '---'

#%%
length_gridbox_cube.data=length_gridbox
length_gridbox_cube.remove_coord('forecast_reference_time')
length_gridbox_cube.remove_coord('forecast_period')
#length_gridbox_cube.remove_coord('time')

#length_gridbox_cube=iris.coords.AuxCoord(length_gridbox,
#                          long_name='length_gridbox',
#                          units='meter^1')#J/(kg·K)
#%%

stash_code='m01s00i254'#mass_fraction_of_cloud_liquid_water_in_air
liquid_water_mmr=iris.load(ukl.Obtain_name(folder,stash_code))[0]
liquid_water_mc=air_density*liquid_water_mmr
liquid_water_mc.long_name='Mass_concentration_of_cloud_liquid_water'
liquid_water_mc._var_name='mcon_lw'
save_cube(liquid_water_mc)
LWP_column=np.empty(liquid_water_mc.data.shape[0]).tolist()

for i in range(liquid_water_mc.data.shape[0]):
    LWP_column[i]=(liquid_water_mc[i,]*length_gridbox_cube)

LWP_cube_list=iris.cube.CubeList(LWP_column)
LWP=LWP_cube_list.merge()[0]
LWP=LWP.collapsed(['model_level_number'],iris.analysis.SUM)
LWP._var_name='LWP'
LWP.long_name='Liquid water path'
save_cube(LWP)


stash_code='m01s00i012'#mass_fraction_of_cloud_ice_in_air
ice_water_mmr=iris.load(ukl.Obtain_name(folder,stash_code))[0]
ice_water_mc=air_density*ice_water_mmr
ice_water_mc.long_name='Mass_concentration_of_cloud_ice'
ice_water_mc._var_name='mcon_iw'
save_cube(ice_water_mc)
IWP_column=np.empty(ice_water_mc.data.shape[0]).tolist()

for i in range(ice_water_mc.data.shape[0]):
    IWP_column[i]=(ice_water_mc[i,]*length_gridbox_cube)

IWP_cube_list=iris.cube.CubeList(IWP_column)
IWP=IWP_cube_list.merge()[0]
IWP=IWP.collapsed(['model_level_number'],iris.analysis.SUM)
IWP._var_name='IWP'
IWP.long_name='Ice water path'
save_cube(IWP)


#CLOUD TOP TEMPERATURE   
cube_l = iris.load(ukl.Obtain_name(folder,'m01s00i254'))[0]
cube_i = iris.load(ukl.Obtain_name(folder,'m01s00i012'))[0]
cloud_mass=cube_l.data[:,:,:,:]+cube_i.data[:,:,:,:]
cloud_mass[cloud_mass<1e-6]=0
temp_cloud=temperature.copy()
temp_cloud_data=temp_cloud.data
temp_cloud_data[cloud_mass==0]=999
temp_cloud.data=temp_cloud_data
temp_cloud_top=temp_cloud.collapsed(['model_level_number'],iris.analysis.MIN)
temp_cloud_bottom=temp_cloud.collapsed(['model_level_number'],iris.analysis.MAX)

temp_cloud_top_data=temp_cloud_top.data
temp_cloud_top_data[temp_cloud_top_data==999]=np.nan
temp_cloud_top.data=temp_cloud_top_data
temp_cloud_top._var_name='CTT'
temp_cloud_top.long_name='Cloud_top_temperature'


temp_cloud_bottom_data=temp_cloud_bottom.data
temp_cloud_bottom_data[temp_cloud_bottom_data==999]=np.nan
temp_cloud_bottom.data=temp_cloud_bottom_data
temp_cloud_bottom._var_name='CTT'
temp_cloud_bottom.long_name='Cloud_bottom_temperature'

save_cube(temp_cloud_top)
save_cube(temp_cloud_bottom)

# CDNC max cloud water   
CDNC_max_cloud_water=CDNC.copy()

CDNC_max_cloud_water=CDNC_max_cloud_water.collapsed(['model_level_number'],iris.analysis.MAX)

cdnc=CDNC.data
lw=cube_l.data
args=np.argmax(lw,axis=1)
data=np.zeros(args.shape)
for it in range(data.shape[0]):
    # print it
    for ilat in range(data.shape[1]):
        for ilon in range(data.shape[2]):
            data[it,ilat,ilon]=cdnc[it,args[it,ilat,ilon],ilat,ilon]

CDNC_max_cloud_water.data=data
CDNC_max_cloud_water._var_name='CDNC_max_cloud_water'
CDNC_max_cloud_water.long_name='Cloud_droplet_concentratio_at_maximum_cloud_water_content'

save_cube(CDNC_max_cloud_water)
