# -*- coding: utf-8 -*-
"""
Created on Fri Feb 06 09:50:48 2015

@author: Windows User
"""
from __future__ import print_function
from mpl_toolkits.basemap import Basemap
import scipy.interpolate as il
#from images2gif import writeGif
from PIL import Image
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.dates as dt
import numpy as np
import os
import glob
import datetime
import time
import sys
import matplotlib.mlab as ml


def xyz2grd(x,y,z):
#    x,y,z = np.loadtxt(finp, unpack = True)
    xmin,xmax,ymin,ymax = [min(x),max(x),min(y),max(y)]
    
    #size of 1 m grid
    nx = (int(xmax - xmin + 1))
    ny = (int(ymax - ymin + 1))

    #generate a regular grid
    xi = np.linspace(xmin,xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    xi, yi = np.meshgrid(xi, yi)



    #interpolate 1
    #zi = ml.griddata(x,y,z,xi,yi,interp = 'linear')

    #interpolate 2
    zi = il.griddata((x, y), z, (xi, yi),method='nearest')
    return xi,yi,zi


def surgemap(datafolder, mapfolder, tyname, start_datetime, chunk):
    surgefilename_list = sorted(glob.glob(os.path.join(datafolder, 'zz*.xyz')))
    uwindfilename_list = sorted(glob.glob(os.path.join(datafolder, 'wu*.xyz')))
    vwindfilename_list = sorted(glob.glob(os.path.join(datafolder, 'wv*.xyz')))
    timestep = start_datetime

    count = 0
    start_calc = time.time()
    file_count = len(surgefilename_list[::chunk])
    for (surgefilename, uwindfilename, vwindfilename) in \
            zip(surgefilename_list[::chunk],
                uwindfilename_list[::chunk],
                vwindfilename_list[::chunk]
               ):
        longitude = []
        latitude = []
        surge = []
        uwind = []
        vwind = []
        figpath = os.path.join(mapfolder, datetime.datetime.strftime(timestep,
            format='%Y-%m-%d_%H:%M')+'.png')
        #~ if os.path.exists(figpath) and timestep != start_datetime:
            #~ timestep = timestep+datetime.timedelta(hours=1)
            #~ continue
        with open(surgefilename, 'r') as surgefile:
            surgelines = surgefile.readlines()
            for i in range(len(surgelines)):
                surgelines[i] = surgelines[i].strip('\n').split()
                longitude.append(float(surgelines[i][0]))
                latitude.append(float(surgelines[i][1]))
                surge.append(float(surgelines[i][2]))

        with open(uwindfilename, 'r') as uwindfile:
            uwindlines = uwindfile.readlines()
            for i in range(len(uwindlines)):
                uwindlines[i] = uwindlines[i].strip('\n').split()
                uwind.append(float(uwindlines[i][2]))

        with open(vwindfilename, 'r') as vwindfile:
            vwindlines = vwindfile.readlines()
            for i in range(len(vwindlines)):
                vwindlines[i] = vwindlines[i].strip('\n').split()
                vwind.append(float(vwindlines[i][2]))

        #~ print ('Done reading', os.path.basename(surgefilename), \
              #~ os.path.basename(uwindfilename), \
              #~ os.path.basename(vwindfilename), '...')

        if timestep == start_datetime:
            print ('Plotting basemap...')
            minX = min(longitude)
            minY = min(latitude)
            maxX = max(longitude)
            maxY = max(latitude)
            surgemap = Basemap(llcrnrlon=minX, llcrnrlat=minY,
                               urcrnrlon=maxX, urcrnrlat=maxY,
                               projection='merc', resolution='f'
                              )
            surgemap.drawcoastlines()
            surgemap.drawmapboundary()
            meridians = range(int(minX), int(maxX) + 2, 2)
            parallels = range(int(minY), int(maxY) + 2, 2)
            surgemap.drawmeridians(meridians, labels=[0, 0, 0, 1])
            surgemap.drawparallels(parallels, labels=[1, 0, 0, 0])
            surgemap.fillcontinents(color='g')
            zoomIn = raw_input('Zoom  in to a specific region? (y/n): ')
            if zoomIn.lower() == 'y':
                crop_minX = float(raw_input('Left longitude (decimal-degree): '))
                crop_maxX = float(raw_input('Right longitude (decimal-degree): '))
                crop_minY = float(raw_input('Bottom latitude (decimal-degree): '))
                crop_maxY = float(raw_input('Top latitude (decimal-degree): '))
                x1 = ((crop_minX - minX) / (maxX - minX)) * plt.xlim()[1]
                x2 = ((crop_maxX - minX) / (maxX - minX)) * plt.xlim()[1]
                y1 = ((crop_minY - minY) / (maxY - minY)) * plt.ylim()[1]
                y2 = ((crop_maxY - minY) / (maxY - minY)) * plt.ylim()[1]
            print('Done plotting basemap...')

        _, _, Z = xyz2grd(longitude, latitude, surge, quality=100)
        _, _, Uwind = xyz2grd(longitude, latitude, uwind, quality=20)
        _, _, Vwind = xyz2grd(longitude, latitude, vwind, quality=20)
        ny = Z.shape[0]
        nx = Z.shape[1]
        windny = Uwind.shape[0]
        windnx = Uwind.shape[1]
        # get lat/lons of ny by nx evenly space grid.
        lons, lats = surgemap.makegrid(nx, ny)
        # get lat/lons of ny by nx evenly space grid.
        windlons, windlats = surgemap.makegrid(windnx, windny)
        # compute map proj coordinates.
        x_map, y_map = surgemap(lons, lats)
        # compute map proj coordinates.
        windx_map, windy_map = surgemap(windlons, windlats)
        minlevel = -2.00
        maxlevel = 2.00
        levelstep = 0.10
        left_clevs = np.arange(minlevel - levelstep, 0.0, levelstep)
        right_clevs = np.arange(0.0, maxlevel + levelstep, levelstep)
        cs_left = surgemap.contourf(x_map, y_map, Z, left_clevs,
                      cmap=plt.cm.Blues_r)
        cs_right = surgemap.contourf(x_map, y_map, Z, right_clevs,
                      extend='max', cmap=plt.cm.hot_r)
        if zoomIn.lower() == 'y':
            plt.xlim((x1, x2))
            plt.ylim((y1, y2))
        vectors = surgemap.quiver(windx_map, windy_map, Uwind, Vwind,
                      cmap=plt.cm.Blues_r, pivot='middle',
                      scale=500.0, zorder=100, headwidth=5)
        title = plt.title(tyname.title() + ' ' +
                    datetime.datetime.strftime(timestep,
                        format='%b %d, %Y %I:%M %p'))
        if timestep == start_datetime:
            cbar_right = surgemap.colorbar(cs_right,
                             location='right', pad='10%')
            cbar_right.set_label('m', rotation=0)
            plt.quiverkey(vectors, -0.25, 0.50, 50, '50 knots',
                            labelpos='S', coordinates='axes')
            plt.tight_layout()
        plt.savefig(figpath)
        vectors.remove()
        title.set_text('')
        timestep = timestep + datetime.timedelta(hours=1)
        count += 1
        percent = count * 100.0 / file_count
        now_calc = time.time()
        remaining_time = (100.0 - percent) * (now_calc - start_calc) / percent
        sys.stdout.write('***** %.2f %% finished. Estimated time left: %.2f seconds *****\r'
            % (percent, remaining_time))
        sys.stdout.flush()
    end_calc = time.time()
    duration = end_calc - start_calc
    print('***** 100 %% finished. Total duration: %.2f seconds *****\r'
        % duration)
    plt.clf()
    
    
def timeseries(datafolder, stationfolder, timeseriesfolder, start_datetime, timestep=600):
    surgefilename_list = sorted(glob.glob(os.path.join(datafolder, 'zz*.xyz')))
    file_count = len(surgefilename_list)
    stationfilepath = os.path.join(stationfolder, 'station.csv')
    with open(stationfilepath, 'r') as stationfile:
        sta_lines = stationfile.readlines()
        sta_lines = [line.strip('\n').split('\t') for line in sta_lines]
        sta_lines = np.array(sta_lines)
        sta_lon = sta_lines[:, 4].astype(np.float)
        sta_lat = sta_lines[:, 5].astype(np.float)
        sta_namelist = sta_lines[:, 0:4]
        xi = np.array([sta_lon, sta_lat]).transpose()
    first = True
    present_time = start_datetime
    count = 0
    start_calc = time.time()
    for surgefilename in surgefilename_list:
        with open(surgefilename, 'r') as surgefile:
            surge_lines = surgefile.readlines()
            surge_lines = [line.strip('\n').split() for line in surge_lines]
            surge_lines = np.array(surge_lines)
            surge_lon = surge_lines[:, 0]
            surge_lat = surge_lines[:, 1]
            surge_val = surge_lines[:, 2]
        points = np.array([surge_lon, surge_lat]).transpose()
        values = np.array([surge_val]).reshape(len(surge_val), 1)
        interpolated_values = interpolate.griddata(points, values, xi,
            method='cubic')
        for (sta_name, value) in zip(sta_namelist, interpolated_values):
            txtfile = os.path.join(timeseriesfolder, ('_').join(
                          sta_name) + '.txt')
            if first:
                with open(txtfile, 'w') as outfile:
                    header_line = (' ').join(sta_name)
                    outfile.write(header_line + '\n')
                    line = datetime.datetime.strftime(present_time,
                        format='%Y/%m/%d %H:%M') + '\t' + str(round(value[0], 3))
                    outfile.write(line + '\n')
            else:
                with open(txtfile, 'a') as outfile:
                    line = datetime.datetime.strftime(present_time,
                        format='%Y/%m/%d %H:%M') + '\t' + str(round(value[0], 3))
                    outfile.write(line + '\n')
        count += 1
        percent = count * 100.0 / file_count
        now_calc = time.time()
        remaining_time = (100.0 - percent) * (now_calc - start_calc) / percent
        sys.stdout.write('***** %.2f %% finished. Estimated time left: %.2f seconds *****\r'
            % (percent, remaining_time))
        sys.stdout.flush()
        present_time = present_time + datetime.timedelta(seconds=timestep)
        first = False
    end_calc = time.time()
    duration = end_calc - start_calc
    print('***** 100 %% finished. Total duration: %.2f seconds *****\r'
        % duration)
        
        
def timeseries_plot(timeseriesfolder):
    file_count = len(glob.glob(os.path.join(timeseriesfolder, '*.txt')))
    count = 0
    start_calc = time.time()
    for txtfile in sorted(glob.iglob(os.path.join(timeseriesfolder, '*.txt'))):
        with open(txtfile, 'r') as timeseriesfile:
            data_lines = timeseriesfile.readlines()
            data_lines = [data.strip('\n').split() for data in data_lines]
            y = [row[2] for row in data_lines[1:]]
            if 'nan' in y:
                continue
            x = [datetime.datetime.strptime(row[0] + row[1],
                    '%Y/%m/%d%H:%M') for row in data_lines[1:]]
            x = dt.date2num(x)
            hfmt = dt.DateFormatter('%m/%d %H:%M')
            plt.plot_date(x, y, 'b-', linewidth='1.5')
            plt.gca().xaxis.set_major_formatter(hfmt)
            plt.xticks(rotation=70)
            plt.ylabel('Storm Surge (m)')
            filename = os.path.basename(txtfile)
            title = filename[:-4].replace('_', ' ')
            plt.title(title)
            plt.tight_layout()
            plt.savefig(os.path.join(timeseriesfolder, filename[:-4] + '.png'))
            count += 1
            percent = count * 100.0 / file_count
            now_calc = time.time()
            remaining_time = (100.0 - percent) * (now_calc - start_calc) / percent
            sys.stdout.write('***** %.2f %% finished. Estimated time left: %.2f seconds *****\r'
                % (percent, remaining_time))
            sys.stdout.flush()
            plt.clf()
    end_calc = time.time()
    duration = end_calc - start_calc
    print('***** 100 %% finished. Total duration: %.2f seconds *****\r'
        % duration)