#DO NOT RUN this script on the gateway ;-)
#When the model produces monthly mean output, this will be much faster, one needs
# then to concatenate and to save, but not to aggregate, which is the slow part
import numpy as np
import iris
import iris.coord_categorisation
import sys
import variable_dict_TOMCAT_UKCA as vdTOMCAT
from glob import glob
import time
import multiprocessing

# Hamish 28/04/17
# The nc_diag_files is not currently fully compatible with the main output file, because monthly means are 
#only written for the main output. To allow L1_processing to work with just the main output, I have added 
#processing for the pressure and temperature variables from the main output.
# Also variable_dict_TOMCAT_UKCA.py has a new array which allows conversion of the TOMCAT output to mmr for
# full compatibility with UKCA. 
# It is now optional to produce netcdf files of budget diagnostics.

#replace by output_files_directory
files_directory='/nfs/a201/earhg/CLOUD/GLOMAP-nc/Hamish-nitrate-8p1/'
path_out='/nfs/a201/earhg/CLOUD/GLOMAP-nc/Hamish-nitrate-8p1/'
prefix_out='nitrate_testsFeb1_'
include_diagnostics=1

nc_diag_files = glob(files_directory+'gld305_diag_t042*2007*nc')
nc_files = glob(files_directory+'gld305*nc')
bigcubelist = iris.load(nc_files[0])
print bigcubelist

# order of components is SU, BC, OC, SS, DU

#in my data, Aerosol104, Aerosol110, and Aerosol114-117 are all zeroes
# So Aerosol092 is nucleation N, 093 is nucleation SU, 094 is nucleation OC
# 095 is Aitken N, 096 Aitken SU, 97 Aitken BC, 98 Aitken OC
# 99 Acc N, 100-103 SU,BC,OC,SS 
# 104 Acc DU which is zero
# 105 cor nd, 106-109 su,ss,bc,oc; 
# 110 cor-sol dust, zero
# 111 Ait-ins ND, 112-113 Ait-ins BC/OC
# 114/115 and 116/117 Acc-ins and cor-ins dust nd/md
#check with SOA_TO_NDMD in prog.f

mass_of_air=29.0

cubes_to_load = vdTOMCAT.TOMCAT_TO_STASH.keys()
print cubes_to_load

# for one variable from the model, concatenate cubes by time and calculate monthly average. This code for monthly
# averages is slow but should soon not be needed (Richard Pope's output doesn't seem to have the averages, but the 
# model does produce them.
def aggregate_cubes(cubelistToAggregate, variable_dict, mass_for_mmr,out_base,mmflag):
    # give it the TOMCAT name if cannot match to a name in the UKCA dictionary
    if mmflag==1:
        print str(cubelistToAggregate[0].long_name)[:12]
        if str(cubelistToAggregate[0].long_name)[:12]=='Monthly mean':
            getname = str(cubelistToAggregate[0].long_name)
        else:
            getname = (str(cubelistToAggregate[0].long_name))[:-13]
    else:
        getname = (str(cubelistToAggregate[0].long_name))
    ukca_name = variable_dict.get(getname, getname)
    print str(ukca_name)
    newlist = cubelistToAggregate.concatenate_cube()
    iris.coord_categorisation.add_month(newlist,'time',name='month')
    monthmean1 = newlist.aggregated_by(['month'],iris.analysis.MEAN)
    sf=1
    if monthmean1.units=='ppbv':
        sf=1e-9
    if monthmean1.units=='ppmv':
        sf=1e-6
    if monthmean1.units=='pptv':
        sf=1e-12
    mass_of_molecule = mass_for_mmr.get(getname, 1)
    print str(monthmean1.long_name), mass_of_molecule
    if mass_of_molecule==1:
        monthmean=monthmean1
    else:
        monthmean=monthmean1*(sf*mass_of_molecule/mass_of_air)
        monthmean.units='kg/kg'
        
    monthmean.rename(str(ukca_name))
    iris.save(monthmean,path_out+out_base+str(ukca_name)+'.nc')

# Read in cubes (call once for each set of netcdf files)
def read_cubes(nc_files, variable_dict,mass_for_mmr,diag_list,path_out,out_base, mmflag, include_diagnostics):
    listtoload=variable_dict.keys()
    if include_diagnostics==1:
        listtoload =listtoload+diag_list.keys()
    listcopy=[]
    if mmflag==1:
        for item in listtoload:
            if item[:5]=='Month':
                listcopy.append(item)
            else:
                item = item + " monthly mean"
                listcopy.append(item)
        listtoload =listcopy
    cubes = iris.load(nc_files[0],listtoload)
    modeloutputfields=[]
    doconcatenation=[]
    for cube in cubes:
        # for each model variable, if it's time dependent add it to a list of time-dependent
        # cubes which is then aggregated into a monthly mean. If it's not time-dependent, save it
        concatenateflag=1
       
        try:
            cube.coord_dims('time')
        except iris.exceptions.CoordinateNotFoundError:
            concatenateflag=0
            ukca_name = variable_dict.get(str(cube.long_name), str(cube.long_name))
            cube.rename(str(ukca_name))
            iris.save(cube,path_out+out_base+str(ukca_name)+'.nc')
        if concatenateflag==1:
            modeloutputfields.append(iris.cube.CubeList([cube]))
        doconcatenation.append(concatenateflag)

    # Load cubes at subsequent times and add each variable to the list of time-dependent cubes
    # The following code relies on python loading cubes in the same order for each file
    # I think ths is safe due to the sort method in _generate_cubes here https://github.com/SciTools/iris/blob/master/lib/iris/__init__.py
    # but check anyway
    cubeit=0
    for cubefile in nc_files[1:]:
        cubelist = iris.load(cubefile, listtoload)
        cubeit=0
        for cube in cubelist:
            if doconcatenation[cubeit]==1 and cube.name()==((modeloutputfields[cubeit])[0]).name():
                (modeloutputfields[cubeit]).append(cube)
            else:
                print 'skip cube'
                print cube.name()
                #print ((modeloutputfields[cubeit])[0]).name()
            cubeit=cubeit+1

    # Concatenate cubes by time and perform monthly averages
    cubeit=0
    jobs=[]
    start=time.time()
    for cubelistToAgg in modeloutputfields:
        print cubeit
        print str(cubelistToAgg[0].long_name)
        cubeit=cubeit+1
        jobs=[]
        p = multiprocessing.Process(target=aggregate_cubes, args=(cubelistToAgg,variable_dict,mass_for_mmr,out_base,mmflag,))
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()
        
    end=time.time()
    print end-start

# The last argument specifies whether or not the input netcdfs are monthly means. Usually set to 1.

read_cubes(nc_files, vdTOMCAT.TOMCAT_TO_STASH,vdTOMCAT.mass_for_mmr,vdTOMCAT.TOMCAT_TO_STASH_DIAGS,
           path_out,prefix_out,1, include_diagnostics)
#read_cubes(nc_diag_files, vdTOMCAT.TOMCAT_TO_STASH,vdTOMCAT.mass_for_mmr,path_out,'nitrate_GLO305_7_diags_',0)

