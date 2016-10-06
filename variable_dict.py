# -*- coding: utf-8 -*-
"""

Code developed by Jesus Vergara Temprado and Kirsty Pringle

eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""

import UKCA_lib as ukl
import iris

PRESM_A_Summary_File='PRESM_A_Summary_File.txt'

#section=[]
#item=[]
#stash_name=[]
#stash_code=[]


#crs = open(PRESM_A_Summary_File, "r")
#for columns in ( raw.strip().split() for raw in crs ):
#    print columns[0], columns[1], columns[2]
#    section.append(columns[0])
#    item.append(columns[1])
#    stash_name.append(columns[2])
#
#    with ukl.Capturing() as output:
#        print iris.fileformats.pp.STASH(1, columns[0],columns[1])
#    stash_code.append(output[0])
variable_reference_stash={}
variable_reference_name={}

def add_variable(stash_code,name,short_name,long_name,units,description='None'):
    variable_reference_stash[stash_code]=ukl.VariableAttributes(stash_code,name,short_name,long_name,units,description)
    variable_reference_name[name]=ukl.VariableAttributes(stash_code,name,short_name,long_name,units,description)

#%%
#%% Use data from pp file and mapping file to fill the VariableAttributes dictionary.

# VariableAttributes class is defined in UKCA_lib.py

stash_code='empty'
long_name_mapping_file='long_name'
name_mapping_file='name'
units_mapping_file='units'


name = 'empty'
short_name='empty'
long_name='empty'
units='empty'

##KP_Comment:   From Mohit, need to update
mapping_file = 'ukca_stdname_vn92_v2'
STASH_File_From_UMUI='teafy.A.diags_short'

#print mapping_file
#print STASH_File_From_UMUI

f = open('STASHmaster_A', 'r')
for _ in range(12):
    header1 = f.readline()
for line in f:
    #header1 = f.readline()
    line=line.split('|')
    if line[0]=='1':
        st1=str(int(line[1]))
        st2=str(int(line[2]))
        st3=str(int(line[3]))
        name=str(line[4]).replace (" ", "_").replace('/','div')
        #if len(st1)
        if len(st1)==1:
            st1='0'+st1
        if len(st2)==1:
            st2='0'+st2
        if len(st3)==1:
            st3='0'+st3
        if len(st3)==2:
            st3='0'+st3
        st_code='m'+st1+'s'+st2+'i'+st3
        add_variable(st_code,'','',name,'')





#%%

f = open(STASH_File_From_UMUI, 'r')
# Read first 8 lines as header.  There must be a better way to do this!
for _ in range(8):
    header = f.readline()

s = '_'

#umui_table={}
umui_table={'stash':[],'name':[]}

for i, line in enumerate(f):
    #print i, line
    line = line.strip()
    columns = line.split()
    stash = ukl.get_stash_from_numbers('1',columns[0],columns[1])
    umui_table['stash'].append(stash)
    umui_table['name'].append(s.join(columns[2:]))


#%%%

mapping_file_dict={'stash':[],'short_name':[],'long_name':[],'units':[]}

f = open(mapping_file, 'r')
header1 = f.readline()

data = []
for line in f:
    line = line.strip()
    columns = line.split()
    #print columns[0], columns[1], columns[2], columns[3]
    stash=columns[0]
    short_name=columns[1]
    long_name=columns[2]
    units=columns[3]
    add_variable(stash,short_name,short_name,long_name,units)
    mapping_file_dict['stash'].append(columns[0])
    mapping_file_dict['short_name'].append(columns[1])
    mapping_file_dict['long_name'].append(columns[2])
    mapping_file_dict['units'].append(columns[3])








# Add the variables names / stash numbers from the stash file to the following dictionaries:
# variable_reference_name and variable_reference_stash

#for i in range (len (section)):
#    print stash_code, stash_name
#    add_variable(stash_code[i],stash_name[i])

#%%
#add_variable(ukl.get_stash_from_numbers('1','34','101'),'nd_nuc_sol_vmr','')
#add_variable('34103','nd_ait_sol_vmr','')
#add_variable('34107','nd_acc_sol_vmr','')
#add_variable('34113','nd_cor_sol_vmr','')
#add_variable('34119','nd_ait_ins_vmr','')
#add_variable('34122','nd_acc_ins_vmr','')
#add_variable('34124','nd_cor_ins_vmr','')
#
#
#print variable_reference_stash
#print ""
#print variable_reference_name
