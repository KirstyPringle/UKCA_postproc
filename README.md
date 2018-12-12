# UKCA_postproc

Python library for post processing and plotting of .pp files outputed by the UM.

### Introduction.

The Python library is designed to convert output data from the Unified Model (pp format) into netCDF files.  The code converts all diagnosed variables to netCDF, and also calculates additional useful diagnostics (e.g. CCN, N10).  The code works 


The working squeme is the following:
-<img height='700' src='http://www.imageurl.ir/images/04941817947862852344.png'/>

#### Version

The current version is version v0.2, and has code status Amber.  This means the code runs, but is still in the process of being validated and has not yet been published. Further developments will be added to the code.

Documentation, test, and examples will be added soon to the repository.

Please note, although we provide the code, users are responsible for checking the code output.  

Please report any bugs / faults to the code owners.


### Authors and Contributors
The main developers, and current code owners, are Jesus Vergara Temprado (@Numlet) and Kirsty Pringle (@KirstyPringle) and Hamish Gordon (@hgordo). 


### Support or Contact
For any questions contact Jesus Vergara Temprado (eejvt@leeds.ac.uk) or Kirsty Pringle (K.pringle@leeds.ac.uk)


## List of files

UKCA_Control_File.py = 

Settings:

dir_scripts = Location of the python code files
jobID = ID of the job e.g. tefxf for UKCA runs 
input_files_directory= Location of the model files, e.g. netCDF, PP or PDG files
output_files_directory= Location to where processed files will be written


UKCA_Control_USER_config_NAME.py = This script controls the input / output directories.  Users can create thier own (replace NAME with e.g. Kirsty), add to the repository and edit. 

## Logical Options

In UKCA_Control_File.py choose whether to do:




## Instructions for running



On JASMIN run with:  python2.7 UKCA_ControlFile.py



