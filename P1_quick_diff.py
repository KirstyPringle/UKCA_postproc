# -*- coding: utf-8 -*-
"""

Code developed by Robin Stevens

rbnstevens@gmail.com

2018

Make quick difference plots between two different sets of L1 output:
Annual means of differences (absolute and relative)
Zonal means of differences (absolute and relative)

This should help you answer the questions "Are these two simulations
significantly different? In what variables? Where? When? By how much?"

"""

from glob import glob
import iris
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
import numpy as np
import sys

# --- PARAMETERS ---

# These are more or less in the order of most-frequently-changed to least-frequently-changed

saving_folder = '/nfs/a201/earrste/Plots/GLOMAP/'

# Note that differences will be plotted as b minus a
files_dir_a  = '/nfs/a201/earrste/GLOMAP_DATA/Nov2018-clouds13-SimoneHOMNO/L1/'
files_dir_b  = '/nfs/a201/earrste/GLOMAP_DATA/Nov2018-clouds12-SimoneHOMNO-P6/L1/'

vert_level = -2 # vertical level to use for plotting 3D data.


clim_rel = 55. # limit for colourbars for relative difference plots [%]

rel_diff_thresh = 5. # threshold in relative difference for making plots [%]
corr_coef_thresh = 0. # threshold in correlation between relative and absolute differences for making plots
eps = 1e-10 # a small number to prevent divide-by-zeroes. Should be smaller than any variable value of interest.

month_chars = ['J','F','M','A','M','J','J','A','S','O','N','D']

# locations to place lat ticks on zonal mean plots. Must correspond to lat_tick_strs
lat_tick_locs = [-90, -60, -30, 0, 30, 60, 90] 

# strings to identify latitudes on zonal mean plots. Must correspond to lat_tick_locs
lat_tick_strs = ['90$\degree$ S', '60$\degree$ S', '30$\degree$ S', '0$\degree$', '30$\degree$ N', '60$\degree$ N', '90$\degree$ N']

cnticks = 12 # number of ticks on colourbar. Number of colours will be this minus one.

posnegclr = [ # colour scheme as rgb values
   (0.0, 0.0, 0.5),
   (0.3, 0.4, 0.9),
   (0.5, 0.9, 1.0),
   (1.0, 1.0, 1.0),
   (1.0, 0.8, 0.5),
   (0.9, 0.4, 0.3),
   (0.5, 0.0, 0.0),
   ]


# --- INITIALIZE ---

nc_files_a = glob(files_dir_a+'*nc')
nc_files_b = glob(files_dir_b+'*nc')

cmap=matplotlib.colors.LinearSegmentedColormap.from_list('tmp', posnegclr, cnticks-1)

map = Basemap(resolution='c', llcrnrlon=0,llcrnrlat=-90, urcrnrlon=360, urcrnrlat=90)

def drawmap():
   # draw coastlines, country boundaries, fill continents.
   map.drawcoastlines(linewidth=0.25)
   
   # draw the edge of the map projection region (the projection limb)
   map.drawmapboundary()

   # draw lat/lon grid lines every 30 degrees.
   map.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1],fontsize=10,linewidth=0.25) # draw meridians
   map.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0],fontsize=10,linewidth=0.25) #draw parallels

# --- BEGIN LOOPING ---

for nc_file_a in nc_files_a:

    # get the filename without the path
    filename_a = nc_file_a.split('/')[-1]

    for nc_file_b in nc_files_b:

        # get the filename without the path
        filename_b = nc_file_b.split('/')[-1]

        # I'll assume the same naming convention for both files, but I won't 
        # assume that the list sizes are equal or in the same order.
        if filename_a == filename_b:
            cube_a = iris.load(nc_file_a)[0]
            cube_b = iris.load(nc_file_b)[0]

            # Check some reasons why we would want to skip this file
            same_name = (cube_a.var_name == cube_b.var_name)
            same_shape = (cube_a.shape == cube_b.shape)
            valid_dim = ( (cube_a.ndim == 3) or (cube_a.ndim == 4) )

            donotskip = same_name and same_shape and valid_dim

            varname = cube_a.var_name

            if donotskip:
                print(varname)

                if cube_a.ndim == 4:
                    # Cubes have data at multiple vertical levels. Just extract one level.
                    cube_a = cube_a[:, vert_level, :, :]
                    cube_b = cube_b[:, vert_level, :, :]

                abs_diff = cube_b - cube_a
                # add a small number to denominator to prevent divide-by-zeros
                rel_diff = abs_diff / (cube_a + eps) * 100.
                
                # At this point we're done with cube_a and cube_b. Remove them to free up memory.
                cube_a = None
                cube_b = None

                # Check some more reasons why we would skip this variable

                # Are the difference values all NaNs?
                all_NaNs = np.isnan(abs_diff.data).all()

                # Is the largest relative difference still small?
                tiny_rel_diff =  abs(rel_diff.data).max() < rel_diff_thresh

                # Is there no correlation between absolute and relative differences?
                # This would indicate that the large relative differences are probably
                # just small differences between smaller values - potentially noise.
                abs_rel_corr_coef = np.corrcoef(abs_diff.data.flatten(), rel_diff.data.flatten())[0,1]
                no_corr = abs_rel_corr_coef < corr_coef_thresh
                    
                donotskip = not tiny_rel_diff and not no_corr and not all_NaNs

                if donotskip:

                 ann_mean_rel = rel_diff.collapsed(['time'],iris.analysis.MEAN)
                 zonal_mean_rel = rel_diff.collapsed(['longitude'],iris.analysis.MEAN)

                 # Is the largest relative difference after averaging small?
                 max_rel_diff = max( abs(ann_mean_rel.data).max(), abs(zonal_mean_rel.data).max() )
                 tiny_rel_diff = max_rel_diff  < rel_diff_thresh

                 if not tiny_rel_diff:

                    ann_mean_abs = abs_diff.collapsed(['time'],iris.analysis.MEAN)
                    zonal_mean_abs = abs_diff.collapsed(['longitude'],iris.analysis.MEAN)

                    # I'm not using Jesus' plotting functions because I don't want
                    # each plot as a separate figure.

                    lats = ann_mean_abs.coord('latitude').points
                    lons = ann_mean_abs.coord('longitude').points

                    plt.figure(figsize=(7,7))

                    #try:
                    #   figtitle = filename_a[(len(varname)+4):-3]
                    #except:
                    figtitle = varname

                    # calculate an appropriate colourbar limit for absolute difference plots
                    clim_abs = max( -np.percentile(  ann_mean_abs.data, 5),
                                    -np.percentile(zonal_mean_abs.data, 5),
                                    np.percentile(   ann_mean_abs.data, 95),
                                    np.percentile( zonal_mean_abs.data, 95),
                                 )

                    plt.subplot(221)
                    plt.title(figtitle)
                    plt.pcolormesh(lons, lats, ann_mean_abs.data, vmin=-clim_abs, vmax=clim_abs, cmap=cmap)
                    drawmap()
                    cb = plt.colorbar(orientation='horizontal')
                    #print('   abs max ann_mean_abs ' + str(abs(ann_mean_abs.data).max()) )

                    plt.subplot(222)
                    plt.pcolormesh(lons, lats, ann_mean_rel.data, vmin=-clim_rel, vmax=clim_rel, cmap=cmap)
                    drawmap()
                    cb = plt.colorbar(orientation='horizontal')
                    #print('   abs max ann_mean_rel ' + str(abs(ann_mean_rel.data).max()) )

                    plt.subplot(223)
                    plt.pcolormesh(range(1,13), lats, zonal_mean_abs.data.transpose(), vmin=-clim_abs, vmax=clim_abs, cmap=cmap)
                    plt.xticks(range(1,13), month_chars)
                    plt.yticks(lat_tick_locs, lat_tick_strs)
                    cb = plt.colorbar(orientation='horizontal')
                    #print('   abs max zonal_mean_abs ' + str(abs(zonal_mean_abs.data).max()) )

                    plt.subplot(224)
                    plt.pcolormesh(range(1,13), lats, zonal_mean_rel.data.transpose(), vmin=-clim_rel, vmax=clim_rel, cmap=cmap)
                    plt.xticks(range(1,13), month_chars)
                    plt.yticks(lat_tick_locs, lat_tick_strs)
                    cb = plt.colorbar(orientation='horizontal')
                    #print('   abs max zonal_mean_rel ' + str(abs(zonal_mean_rel.data).max()) )

                    saving_str=saving_folder+'diff_plot_'+varname+'.png'
                    plt.savefig(saving_str)

                 else:
                    print(' Skipping '+filename_a)
                    print('  Maximum meaned relative difference is ' + format(max_rel_diff, '0.2f')
                          + '%, < threshold ' + format(rel_diff_thresh, '0.2f') + '%')

                else:
                    print(' Skipping '+filename_a)
                    if tiny_rel_diff:
                       print('  Maximum relative difference is ' + format(abs(rel_diff.data).max(), '0.2f')
                              + '%, < threshold ' + format(rel_diff_thresh, '0.2f') + '%')
                    if no_corr:
                       print('  Correlation between abs and rel diffs is ' + format(abs_rel_corr_coef, '0.2f')
                              + ' < threshold ' + format(corr_coef_thresh, '0.2f') )

                    if all_NaNs:
                       print('  All difference values are NaNs.')

            else:
                print(' Skipping '+filename_a)
                if not same_name:
                    print('  Cubes have different variable names.')
                if not same_shape:
                    print('  Cubes have different shapes.')
                if not valid_dim:
                    print('  Cubes have an unexpected number of dimensions.')

#plt.show()



#
