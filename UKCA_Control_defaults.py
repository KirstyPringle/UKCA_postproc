####IMPORTANT - DO NOT EDIT values in here ###
#Default values for flags, etc.


# This is the directory where the code is. Change it to the path where you downloaded the file
#dir_scripts='/nfs/see-fs-01_users/eejvt/CASIM_postproc/'
dir_scripts='/home/users/kpringle/DataVisualisation/AerVis/UKCA_postproc/'
#dir_scripts='/nfs/a107/earkpr/DataVisualisation/Jesus/git_area/UKCA_postproc-master/'


#  Which model would you like to read data from?  Options are:  UKCA, TOMCAT_GLOMAP or CASIM
model_name = "CASIM"
model_name = "UKCA"
#model_name = "TOMCAT_GLOMAP"

jobID = "whatever"
jobID = "tebiz"
#jobID = "glo301"

run_L0=True                 # if true will convert raw output to Level1 nc files.
#run_L0=0
run_L1=True                 # if true will calculate additional output as a postprocessing step e.g. CCN 
run_plots=False             # if true will produce a series of standard plots
send_mail=True              # ?
## Key user defined paths


#orog_file = '/group_workspaces/jasmin2/gassp/jvergaratemprado/n96_hadgem1_qrparm.orog_new.pp'
#jasminorog_file = '/nfs/a107/earkpr/ACID-PRUFF/Masaru/OAT5/teafw/ppfiles/n96_hadgem1_qrparm.orog_new.pp'#leeds foe-linux


# Location of the UM output files to proces
input_files_directory='/nfs/a201/eejvt/CASIM/SO_KALLI/NO_CLOUD_SQUEME/GP_HIGH_CSED/'
input_files_directory='/nfs/a201/eejvt/CASIM/SECOND_CLOUD/GP_HAM_DMDUST/'    


## Location of where to write Level 0 (data files in nc format). Will be the same as input_files_directory unless the flag below is set to True
l_set_output_dir=False #output dir will be set to output_files_directory_set if True.
#output_files_directory_set='/nfs/a201/'+str(username)+'/'+str(model_name)+'_TEST_FILES/'+str(jobID)+'/'

