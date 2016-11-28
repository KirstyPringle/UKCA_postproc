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

import UKCA_lib as ukl
import matplotlib.pyplot as plt
import matplotlib

'''
Examples of ploting routines

'''


#Load a cube is as simple as:
path_to_nc_file=directory_scripts+'test_cube_PM25.nc'

cube=iris.load(path_to_nc_file)[0]#we add [0] at the end, as iris.load load the cubes in a list

#Let's see how the cube looks like

print cube

#We can see thins such as the shape, number of dimensions etc...
print cube.shape
print cube.ndim

#we can subset our cube, for example, for a single month, let's say, january
cube_january=cube[0,:,:,:]

#we can also calculate the mean through all the time period
cube_time_mean=cube.collapsed(['time'],iris.analysis.MEAN)
#see how the dimensions have changed
print 'cube  shape ',cube.shape
print 'cube_time_mean shape ',cube_time_mean.shape

#if we want to get the surface level:
cube_surface=cube_time_mean[0,:,:]

#lets do a quick plot of the surface
#%%
import iris.quickplot as qplt

qplt.contourf(cube_surface)
#coastlines
plt.gca().coastlines()
#%%
#it doesn't look very good on linear scale, let's do logarithmic scale
plt.figure()

qplt.contourf(cube_surface,norm=matplotlib.colors.LogNorm())
plt.gca().coastlines()


#%%

#The colors do not look very good, let's choose another colorscale

cmap=plt.cm.CMRmap_r#more on : http://matplotlib.org/examples/color/colormaps_reference.html

plt.figure()
qplt.contourf(cube_surface,cmap=cmap,norm=matplotlib.colors.LogNorm())
plt.gca().coastlines()





