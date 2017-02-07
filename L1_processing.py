# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""
import sys
dir_scripts='/nfs/see-fs-01_users/eejvt/UKCA_postproc'
sys.path.append(dir_scripts)#Change this to the downloaded folder
import numpy as np
import iris
sys.path.append('/nfs/a107/eejvt/PYTHON_CODE')
import Jesuslib as jl
import UKCA_lib as ukl
reload(ukl)
import multiprocessing
import variable_dict as vd
reload(vd)
import I_MODE_SETUP_Variables as ims
reload(ims)
from matplotlib.colors import LogNorm

def save_cube(cube):
    saving_name=saving_folder_l1+'L1_'+cube._var_name+'_'+cube.long_name+'.nc'
    iris.save(cube,saving_name, netcdf_format="NETCDF4")
    print 'saved:',cube.long_name

#directory='/nfs/a201/eejvt/UKCA_TEST_FILES/tebxd/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/NO_INITIAL_ICE/BASE_RUN/'
directory='/nfs/a201/eejvt/CASIM/SO_KALLI/TRY2/ALL_ICE_PROC/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/TRY2/3_ORD_LESS_762/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/SECOND_DOMAIN/NOICE/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/TRY2/BASE_CONTACT_242/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/TRY2/BASE_CONTACT_242/'
#folder=directory+'All_time_steps/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/TRY2/2_ORD_MORE/'
#directory='/nfs/a201/eejvt/CASIM/SO_KALLI/CLOUD_SQUEME/BASE/'
directory='/nfs/a201/eejvt/CASIM/SO_KALLI/NO_CLOUD_SQUEME/MEYERS/'
folder=directory+'All_time_steps/'
saving_folder_l1=directory+'L1/'
ukl.create_folder(saving_folder_l1)

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
liquid_water_mc.long_name='Mass concentration of cloud liquid water'
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
ice_water_mc.long_name='Mass concentration of cloud ice'
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

temp_cloud_top_data=temp_cloud_top.data
temp_cloud_top_data[temp_cloud_top_data==999]=np.nan
temp_cloud_top.data=temp_cloud_top_data
temp_cloud_top._var_name='CTT'
temp_cloud_top.long_name='Cloud_top_temperature'

save_cube(temp_cloud_top)



jl.send_email()

#%%
#plt.figure()
#plt.contourf(LWP[13,].data,norm = LogNorm(),cmap=plt.cm.Blues)
#plt.colorbar()
#plt.show()
'''
#%%
#Obtain mass mixing ratios from netcdf
particles_mixing_ratio={}
particles_concentration={}
for name in vd.variable_reference_name.keys():
    if 'nbr_' in name:
        #print name
        nc_name=ukl.Obtain_name(folder,vd.variable_reference_name[name].stash_code)
        if len(nc_name)==0:
            print name, 'not founded'
            continue
        particles_mixing_ratio[name]=iris.load(nc_name)[0]
        particles_concentration['n'+name[3:]]=particles_mixing_ratio[name]*particle_density_of_air
        particles_concentration['n'+name[3:]].long_name=particles_mixing_ratio[name].long_name
        particles_concentration['n'+name[3:]]._var_name='n'+name[3:]
        save_cube(particles_concentration['n'+name[3:]])

        #print particles_concentration.keys()
        print particles_concentration['n'+name[3:]].data[0,0,0,0]


mass_mixing_ratio={}
for name in vd.variable_reference_name.keys():
    if 'mmr' in name:
        if 'sol' in name or 'ins' in name:
            #print name
            nc_name=ukl.Obtain_name(folder,vd.variable_reference_name[name].stash_code)
            if len(nc_name)==0:
                print name, 'not founded'
                continue
            mass_mixing_ratio[name]=iris.load(nc_name)[0]

#Calculate mass concentrations from air density and mass mixing ratios
mass_concentration={}
for name in mass_mixing_ratio.keys():
    new_short_name='mcon_'+name[3:]
    #print new_short_name
    mass_concentration[new_short_name]=air_density*mass_mixing_ratio[name]
    mass_concentration[new_short_name].long_name='mass_concentration'+mass_mixing_ratio[name].long_name[13:-7]
    mass_concentration[new_short_name]._var_name=new_short_name
    save_cube(mass_concentration[new_short_name])
'''
'''
mass_concentration_per_mode={}
for mode_name in ims.mode_names:
    cubes_to_add=[mass_concentration[key] for key in mass_concentration.keys() if mode_name in key]
    mass_concentration_per_mode['mcon_'+mode_name]=np.sum(cubes_to_add)
    mass_concentration_per_mode['mcon_'+mode_name].long_name='total_mass_concentration_in_%s'%mode_name
    mass_concentration_per_mode['mcon_'+mode_name]._var_name='mcon_'+mode_name
    save_cube(mass_concentration_per_mode['mcon_'+mode_name])

mass_concentration_per_component={}
for comp_name in ims.component_names:
    cubes_to_add=[mass_concentration[key] for key in mass_concentration.keys() if comp_name in key]
    mass_concentration_per_mode['mcon_'+comp_name]=np.sum(cubes_to_add)
    mass_concentration_per_mode['mcon_'+comp_name].long_name='total_mass_concentration_of_%s'%comp_name
    mass_concentration_per_mode['mcon_'+comp_name]._var_name='mcon_'+comp_name
    save_cube(mass_concentration_per_mode['mcon_'+comp_name])


#Calculate volume of modes per component
mode_volume_per_component={}
list_volumes=[]
for mode_name in ims.mode_names:
    keys=[key for key in mass_concentration.keys() if mode_name in key]
    #print keys
    for comp_name in ims.component_names:
        keys_comp=[key for key in keys if comp_name in key]
        if any(keys_comp):
            mode_volume_per_component['vol_'+keys_comp[0][5:]]= \
            mass_concentration[keys_comp[0]]/ims.species_attributes[comp_name].rhocomp
            mode_volume_per_component['vol_'+keys_comp[0][5:]].long_name='Volume_fraction'+mass_concentration[keys_comp[0]].long_name[18:]
            mode_volume_per_component['vol_'+keys_comp[0][5:]]._var_name='vol_'+keys_comp[0][5:]
            save_cube(mode_volume_per_component['vol_'+keys_comp[0][5:]])

#Calculate volume and radius of modes
mode_radius={}
mode_volume={}
for mode_name in ims.mode_names:
    keys=[key for key in mode_volume_per_component.keys() if mode_name in key]
    #print mode_name

    cubes_to_add=[mode_volume_per_component[key] for key in keys]
    mode_volume['vol_'+mode_name]=np.sum(cubes_to_add)
    mode_volume['vol_'+mode_name]._var_name='vol_'+mode_name
    mode_volume['vol_'+mode_name].long_name='Total_volume_fraction_'+mode_name
    save_cube(mode_volume['vol_'+mode_name])

    mode_radius['rad_'+mode_name]=0.5*iris.analysis.maths.exponentiate((6*mode_volume['vol_'+mode_name]/particles_concentration['n_'+mode_name])/(np.pi*np.exp(4.5*(np.log(ims.modal_attributes[mode_name].sigma)))),1./3.)
    mode_radius['rad_'+mode_name]._var_name='rad_'+mode_name
    mode_radius['rad_'+mode_name].long_name='Radius_of_mode_'+mode_name
    save_cube(mode_radius['rad_'+mode_name])


# factor to convert from mean Radious in the number PDF to the mean Radious in the volume/mass PDF conv_V = np.exp(3.0*np.log(sigma[imode])**2) #check this
cubes_to_add_N2p5=[]
r_125=1.25e-6#meters long_name='Radious for calculating PM2.5',
for mode_name in ims.mode_names:
    #print mode_name
    N=particles_concentration['n_'+mode_name]
    rbar=mode_radius['rad_'+mode_name].data
    sigma=ims.modal_attributes[mode_name].sigma
    contribution_to_N2p5=ukl.lognormal_cummulative(N,r_125,rbar,sigma)
    cubes_to_add_N2p5.append(contribution_to_N2p5)
N2p5=np.sum(cubes_to_add_N2p5)
N2p5._var_name='N2p5'
N2p5.long_name='Particules_smaller_than_a_diameter_of_2.5_um'
save_cube(N2p5)

cubes_to_add_N10=[]
r_5=5e-6#meters long_name='Radious for calculating N10',
for mode_name in ims.mode_names:
    #print mode_name
    N=particles_concentration['n_'+mode_name]
    rbar=mode_radius['rad_'+mode_name].data
    sigma=ims.modal_attributes[mode_name].sigma
    contribution_to_N10=ukl.lognormal_cummulative(N,r_5,rbar,sigma)
    cubes_to_add_N10.append(contribution_to_N10)
N10=np.sum(cubes_to_add_N10)
N10._var_name='N10'
N10.long_name='Particules_smaller_than_a_diameter_of_10_um'
save_cube(N10)

cubes_to_add_PM25=[]
r_125=1.25e-6#meters long_name='Radious for calculating PM2.5',
for mode_name in ims.mode_names:
    #print mode_name
    Mode_mass=mass_concentration_per_mode['mcon_'+mode_name]
    rbar=mode_radius['rad_'+mode_name].data
    sigma=ims.modal_attributes[mode_name].sigma
    rbar_volume=rbar*np.exp(3.0*np.log(sigma)**2)
    contribution_to_PM25=ukl.lognormal_cummulative(Mode_mass,r_125,rbar_volume,sigma)
    cubes_to_add_PM25.append(contribution_to_PM25)
PM25=np.sum(cubes_to_add_PM25)
PM25._var_name='PM25'
PM25.long_name='Mass_of_particules_smaller_than_a_diameter_of_2.5_um'
save_cube(PM25)

cubes_to_add_PM10=[]
r_5=5e-6#meters long_name='Radious for calculating PM10',
for mode_name in ims.mode_names:
    #print mode_name
    Mode_mass=mass_concentration_per_mode['mcon_'+mode_name]
    rbar=mode_radius['rad_'+mode_name].data
    sigma=ims.modal_attributes[mode_name].sigma
    rbar_volume=rbar*np.exp(3.0*np.log(sigma)**2)
    contribution_to_PM10=ukl.lognormal_cummulative(Mode_mass,r_5,rbar_volume,sigma)
    cubes_to_add_PM10.append(contribution_to_PM10)
PM10=np.sum(cubes_to_add_PM10)
PM10._var_name='PM10'
PM10.long_name='Mass_of_particules_smaller_than_a_diameter_of_10_um'
save_cube(PM10)


#CCN

#This method works for k>0.2, sc<1%
#Petters, M. D., & Kreidenweis, S. M. (2007).
#A single parameter representation of hygroscopic growth and cloud condensation nucleus activity.
#Atmospheric Chemistry and Physics, 7(8), 1961–1971. doi:10.5194/acp-7-1961-2007

surface_tension_water = 0.0761 - 1.55e-4 * (temperature.data - 273.)#J m^-2 Nenes/Seinfeld
Mw=0.018015#molecular weight of water kg/mol
R=8.3144#Universal gas constant J mol-1 K-1
pw=1000#Density of water  kg m-3
A_kohler = 4. * Mw * surface_tension_water / ( R * temperature.data * pw )
kappa_weighted_per_mode={}
for mode_name in ims.mode_names:
    cubes_to_add=[]
    if ims.modal_attributes[mode_name].modesol:
        for comp_name in ims.component_names:
            if 'vol_'+comp_name+'_'+mode_name in mode_volume_per_component.keys():
                print 'kappa_weighted_per_mode',mode_name,comp_name
                cubes_to_add.append(mode_volume_per_component['vol_'+comp_name+'_'+mode_name]*ims.species_attributes[comp_name].kappa/mode_volume['vol_'+mode_name])
                kappa_weighted_per_mode[mode_name]=np.sum(cubes_to_add)

for supersaturation in [0.1,0.2,0.5,1]:
    #print supersaturation
    critical_diameter_per_mode={}
    CCN_per_mode={}
    CCN=0
    for mode_name in ims.mode_names:
        if ims.modal_attributes[mode_name].modesol:
            critical_diameter_per_mode[mode_name]=(4.*A_kohler**3/(27.*kappa_weighted_per_mode[mode_name].data*np.log(1+supersaturation/100.)))**(1./3.)
            N=particles_concentration['n_'+mode_name]
            rbar=mode_radius['rad_'+mode_name].data
            sigma=ims.modal_attributes[mode_name].sigma
            Rcrit=critical_diameter_per_mode[mode_name]/2.
            CCN_per_mode[mode_name]=N-ukl.lognormal_cummulative(N,Rcrit,rbar,sigma)
            CCN=CCN+CCN_per_mode[mode_name]
    CCN.long_name='Cloud_condensation_nuclei_at_a_supersaturation_of_%1.4f'%supersaturation
    CCN._var_name='ccn'+str(supersaturation)
    save_cube(CCN)
'''
