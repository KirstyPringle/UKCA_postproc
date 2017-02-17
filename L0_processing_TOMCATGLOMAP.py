#DO NOT RUN this script on the gateway ;-)
#When the model produces monthly mean output, this will be much faster, one needs
# then to concatenate and to save, but not to aggregate, which is the slow part
import numpy as np
import iris
import iris.coord_categorisation
import sys
import variable_dict_TOMCAT as vdTOMCAT
from glob import glob
import time
import multiprocessing
from __main import *

#replace by input_files_directory
# files_directory='/nfs/a105/eerjp/TOMCAT/data_2012_2015/'
# files_directory='/nfs/a201/earhg/CLOUD/GLOMAP-nc/tests-Jan2/'

files_directory=input_files_directory
#replace by output_files_directory
# path_out='/nfs/a201/earhg/CLOUD/GLOMAP-nc/Richard-2013/'
path_out=files_directory
nc_files=glob(files_directory+'glo301_t042_2006*nc')
nc_diag_files = glob(files_directory+'glo301_diag_t042_2006*nc')

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

cubes_to_load = vdTOMCAT.TOMCAT_TO_STASH.keys()
print cubes_to_load

# for one variable from the model, concatenate cubes by time and calculate monthly average. This code for monthly
# averages is slow but should soon not be needed (Richard Pope's output doesn't seem to have the averages, but the
# model does produce them.
def aggregate_cubes(cubelistToAggregate, variable_dict, out_base,mmflag):
    # give it the TOMCAT name if cannot match to a name in the UKCA dictionary
    if mmflag==1:
        getname = (str(cubelistToAggregate[0].long_name))[:-13]
    else:
        getname = (str(cubelistToAggregate[0].long_name))
    ukca_name = variable_dict.get(getname, getname)
    print str(ukca_name)
    print cubelistToAggregate
    newlist = cubelistToAggregate.concatenate_cube()
    iris.coord_categorisation.add_month(newlist,'time',name='month')
    # The GLOMAP code will ultimately do this for us; this code should just work faster when you give it monthly
    # averages as input (not tested)
    monthmean = newlist.aggregated_by(['month'],iris.analysis.MEAN)
    print monthmean
    print 'concatenated'
    #if mmflag==1:
    #    writename = (str(ukca_name))[:-13]
    #else:
    #    writename = str(ukca_name)
    iris.save(monthmean,path_out+out_base+str(ukca_name)+'.nc')

# Read in cubes (call once for each set of netcdf files)
def read_cubes(nc_files, variable_dict,path_out,out_base, mmflag):
    listtoload=variable_dict.keys()
    listcopy=[]
    if mmflag==1:
        for item in listtoload:
            item = item + " monthly mean"
            listcopy.append(item)
        listtoload =listcopy
    cubes = iris.load(nc_files[0],listtoload)
    print cubes
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
        p = multiprocessing.Process(target=aggregate_cubes, args=(cubelistToAgg,variable_dict,out_base,mmflag,))
        jobs.append(p)
        p.start()

    for job in jobs:
        job.join()

    end=time.time()
    print end-start

# The last argument specifies whether or not the input netcdfs are monthly means
read_cubes(nc_files, vdTOMCAT.TOMCAT_TO_STASH,path_out,'monthmean_2013_testsJan2_',1)
#read_cubes(nc_diag_files, vdTOMCAT.TOMCAT_TO_STASH,path_out,'monthmean_2013_testsJan2diag_',0)
