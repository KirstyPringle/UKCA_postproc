"""

Routine to store model variables not given in the ouput files.

Variables are set within the modules file, depending on the I_MODE_SETUP used


"""

import UKCA_lib as ukl
reload(ukl)
import matplotlib.pyplot as plt
import iris

I_MODE_SETUP=8#is this the number for this setup? I set it to 8 as example
#we can add later on the I_MODE_SETUP to be read as an argument when you run the script

if I_MODE_SETUP==8:
    # Modal based paramters that are dependent on the I_MODE_SETUP of the model

    mode_names=["nucsol","aitsol","accsol","corsol","aitins","accins","corins"]

    # Mode switches (1=on, 0=0ff)
    mode_choice=[1,1,1,1,1,1,1]

    # Specify which modes are soluble
    # Mode names
    modesol=[1,1,1,1,0,0,0]

    # Mode width
    sigma=[1.59,1.59,1.40,2.0,1.59,1.59,2.0]

    # Specify size limits of geometric mean diameter for each mode (ddp - dry diameter of particle, m)
    # Lower Limit of Mode: ddplim0
    ddplim0=[1.0e-9,1.0e-8,1.0e-7,0.5e-6,1.0e-8,1.0e-7,1.0e-6]
    # Upper Limit of Mode: ddplim1
    ddplim1=[1.0e-8,1.0e-7,0.5e-6,1.0e-5,1.0e-7,1.0e-6,1.0e-5]


    # Species (Component) based paramters that are dependent on the I_MODE_SETUP of the model

    # Component names
    component_names= ["so4","bc","oc","ss","dust"]#,"sec_org"]#do we have secondary organics??
    kappa_component=[   0.61,    0.00,    0.10,    1.28,    0.00]#, 0.1]#check seconday organics I setted it to the same as organics but need to check. Jesus
    # Molar Mass kg / mol
    mm= [0.098,0.012,0.0168,0.05844,0.100]#,0.0168]

    # Density
    rhocomp=[1769.0,1500.0,1500.0,1600.0,2650.0]#,1500.0]#kg/m3

# Put species (component) and modal parameters into dictionaries
species_attributes={}
for i in range (len (component_names)):
    species_attributes[component_names[i]]=ukl.SpeciesAttributes(component_names[i],mm[i],
    iris.coords.AuxCoord(rhocomp[i],long_name='density',units='kilograms-meter^-3'),
    kappa_component[i])

modal_attributes={}
for i in range (len (mode_names)):
    modal_attributes[mode_names[i]]=ukl.ModalAttributes(mode_names[i],sigma[i],ddplim0[i],ddplim1[i],modesol[i],mode_choice[i])

#  To get the values out use, e.g.:
#print species_attributes['so4'].name
#print species_attributes['so4'].mm
#print species_attributes['so4'].rhocomp

#print modal_attributes['nuc_sol'].name
#print modal_attributes['nuc_sol'].sigmag
#print modal_attributes['acc_sol'].ddplim1

#Section to have a overview of how the modes looks like
plot=0#to plot, set to non 0
if plot:
    plt.figure(figsize=(15,10))
    ax=plt.subplot(211)
    for mode in modal_attributes.itervalues():
        print mode
        mode.plot_mode((mode.ddplim1-mode.ddplim0)/2.,limits=1)
    plt.legend()
    #plt.show()
    bx=plt.subplot(212)
    for mode in modal_attributes.itervalues():
        print mode
        mode.plot_mode((mode.ddplim1-mode.ddplim0)/2.,real_PDF=1)
    plt.legend()
    plt.show()

#modal_attributes['acc_sol'].plot_mode(0.5e-6)
