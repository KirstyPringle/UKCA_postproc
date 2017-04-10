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
# air_pressure=iris.load(ukl.Obtain_name(folder,'m01s00i408'))[0]
try:
    air_pressure=iris.load(ukl.Obtain_name(folder,'m01s00i408'))[0]
    if air_pressure.shape!=potential_temperature.shape:
        raise NameError('air pressure with different shape as potential_temperature')
except:
    air_pressure_nodim=iris.load(ukl.Obtain_name(folder,'m01s00i255'))[0]
    p_convert = iris.coords.AuxCoord(100000.0,
                              long_name='convert_units',
                              units='Pa')
    air_pressure=air_pressure_nodim*p_convert



print 'potential temperature and air_pressure loaded'
p0 = iris.coords.AuxCoord(1000.0,
                          long_name='reference_pressure',
                          units='hPa')
p0.convert_units(air_pressure.units)

Rd=287.05 # J/kg/K
cp=1005.46 # J/kg/K
Rd_cp=Rd/cp
print 'constants defined'
temperature=potential_temperature*(air_pressure/p0)**(Rd_cp)
# temperature=potential_temperature.data*(air_pressure.data/1000)**(Rd_cp)

print 'temperature calculated'

temperature._var_name='temperature'
temperature.long_name='Temperature'
# Need to print a value to force the code to calculate the cube
# Do not remove call to ukl.print_cube_single_value.

ukl.print_cube_single_value(temperature)
save_cube(temperature)



temperature._var_name='temperature'
R_specific=iris.coords.AuxCoord(287.058,
                          long_name='R_specific',
                          units='J-kilogram^-1-kelvin^-1')#J/(kg·K)

air_density=(air_pressure/(temperature*R_specific))

# Need to print a value to force the code to calculate the cube
ukl.print_cube_single_value(air_density)  # Do not remove


#%%
molar_mass_air=iris.coords.AuxCoord(28.991e-3,
                          long_name='Molar mass of air',
                          units='kilogram-mole^-1')#J/(kg·K)
avogadro_number=iris.coords.AuxCoord(6.022e23,
                          long_name='Avogadros number - particles per mol',
                          units='mole^-1')#J/(kg·K)

particle_density_of_air=air_density/molar_mass_air*avogadro_number

ukl.print_cube_single_value(particle_density_of_air) # Do not remove

air_density._var_name='air_density'
air_density.long_name='Density of air'
save_cube(air_density)

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

        print particles_concentration.keys()
        ukl.print_cube_single_value(particles_concentration['n'+name[3:]]) # Do not remove


mass_mixing_ratio={}
for name in vd.variable_reference_name.keys():
    if 'mmr' in name:
        if '_' in name:
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
        if len(keys_comp):
            mode_volume_per_component['vol_'+keys_comp[0][5:]]= \
            mass_concentration[keys_comp[0]]/ims.species_attributes[comp_name].rhocomp
            mode_volume_per_component['vol_'+keys_comp[0][5:]].long_name='Volume_fraction'+mass_concentration[keys_comp[0]].long_name[18:]
            mode_volume_per_component['vol_'+keys_comp[0][5:]]._var_name='vol_'+keys_comp[0][5:]
            save_cube(mode_volume_per_component['vol_'+keys_comp[0][5:]])

#%%
#Calculate volume and radius of modes
mode_radius={}
mode_volume={}

for mode_name in ims.mode_names:
    keys=[key for key in mode_volume_per_component.keys() if mode_name in key]
    print mode_name

    cubes_to_add=[mode_volume_per_component[key] for key in keys]
    mode_volume['vol_'+mode_name]=np.sum(cubes_to_add)
    mode_volume['vol_'+mode_name]._var_name='vol_'+mode_name
    mode_volume['vol_'+mode_name].long_name='Total_volume_fraction_'+mode_name
    save_cube(mode_volume['vol_'+mode_name])

    mode_radius['rad_'+mode_name]=0.5*iris.analysis.maths.exponentiate((
        6*mode_volume['vol_'+mode_name]/particles_concentration['n_'+mode_name]
        )/(np.pi*np.exp(4.5*(np.log(ims.modal_attributes[mode_name].sigma))))
        ,1./3.)


    mode_radius['rad_'+mode_name]._var_name='rad_'+mode_name
    mode_radius['rad_'+mode_name].long_name='Radius_of_mode_'+mode_name
    save_cube(mode_radius['rad_'+mode_name])

#%%

# factor to convert from mean Radius in the number PDF to the mean Radius in the volume/mass PDF conv_V = np.exp(3.0*np.log(sigma[imode])**2) #check this
cubes_to_add_N2p5=[]
r_125=1.25e-6#meters long_name='Radius for calculating PM2.5',
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
r_5=5e-6#meters long_name='Radius for calculating N10',
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

#%%

cubes_to_add_PM25=[]
r_125=1.25e-6#meters long_name='Radius for calculating PM2.5',
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
r_5=5e-6#meters long_name='Radius for calculating PM10',
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
