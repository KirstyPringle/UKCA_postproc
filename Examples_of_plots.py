# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:53:19 2016

@author: eejvt
"""

import iris
import sys

##CHANGE THE PATH TO THE DIRECTORY WHERE YOUR SCRIPTS ARE
directory_scripts='/nfs/see-fs-01_users/eejvt/tests/UKCA_postproc/'
sys.path.append(directory_scripts)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import os
cwd = os.getcwd()+'/'
if not os.path.isfile(directory_scripts+'UKCA_lib.py'):
    raise NameError(directory_scripts+'UKCA_lib.py does not exist \n Change directory_scripts to the folder where you downloaded UKCA_postproc ')

'''
I recomend using Spyder for going through this tutorial

In Spyder, you can run single lines with F9 

or whole blocks separated by #%% with Control+Enter (or Control+Shift+Enter)


'''



'''
Examples of ploting routines

'''


#Loading a cube is as simple as:

#set the path to the file:
path_to_nc_file=directory_scripts+'test_cube_PM25.nc'#Run this line in Spyder with F9 and the following.
cube=iris.load(path_to_nc_file)[0]#and loading it. We add [0] at the end, as iris.load load the cubes in a list.
#Alternativelly you can use iris.load_cube(path_to_nc_file) for just a single cube

#Let's see how the cube looks like

print cube

#We can check things such as the shape, number of dimensions etc...
print cube.shape
print cube.ndim

#we can subset our cube, for example, for a single month, let's say, january (0 as in python lists indexes start in 0 not 1)
cube_january=cube[0,:,:,:]

#we can also calculate the mean through all the time period
cube_time_mean=cube.collapsed(['time'],iris.analysis.MEAN)

#see how the dimensions have changed
print 'cube  shape ',cube.shape
print 'cube_time_mean shape ',cube_time_mean.shape

#if we want to get the surface level:
cube_surface=cube_time_mean[0,:,:]
#This is specifically of the sample cube included in the repo.
#For other formats of nc files, the surface level might be other than 0 

#%%
#lets do a quick plot of the surface
#In spyder you can run this whole block with Control+Enter (Comand+Enter in mac)

import iris.quickplot as qplt

plt.figure()
qplt.contourf(cube_surface)

# add coastlines
plt.gca().coastlines()
plt.show()
#%%
#it doesn't look very good on linear scale, let's do logarithmic scale
plt.figure()

qplt.contourf(cube_surface,norm=matplotlib.colors.LogNorm())
plt.gca().coastlines()
plt.show()


#%%

#The colors do not look very good, let's choose another colorscale

cmap=plt.cm.CMRmap_r#more on : http://matplotlib.org/examples/color/colormaps_reference.html

plt.figure()
qplt.contourf(cube_surface,cmap=cmap,norm=matplotlib.colors.LogNorm())
plt.gca().coastlines()
plt.show()

#%%


#Now lets use the functions of UKCA_lib
import UKCA_lib as ukl
reload(ukl)# This bit is necessary if you have modified UKCA_lib

#ukl.level_plot

print 'Current working directory:', os.getcwd()+'/'
print 'The files will be saved there'
saving_path=os.getcwd()+'/'


#Plotting something is as easy as:
ukl.level_plot(cube_time_mean,saving_path)
#ukl.level_plot(cube,saving_path)
#Now the file will be in saving_path
print saving_path
#go to the saving_path there and check it

#If you are using Spyder, you can look at the documentation any function
#by setting the cursor in the name of the function 
# and pressin Control+I
#you should see somethig like:

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


#let's change some of the optional parameters (kwargs) of the function


lev=22#this variable sets the level 
name='personalized_plot'# this will be added to the file name
color_levels=15#more levels in the colors
cmap=plt.cm.RdBu_r#color map
saving_format='.jpg'#saving format
ukl.level_plot(cube_time_mean,saving_path,name=name,level=lev,color_levels=color_levels,cmap=cmap,saving_format=saving_format)

#Now in logarithmic scale. We need to set logscale=True (or anything different from 0 (False))

ukl.level_plot(cube_time_mean,saving_path,name=name,level=lev,cmap=cmap,logscale=True)










