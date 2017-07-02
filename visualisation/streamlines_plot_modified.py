# This script allows plotting of contour levels for chosen macroscopic
# variable.

import os # To get access for folders
import re # To use regular expressions
import numpy as np

# To plot haetmap and bodys boundary
from numpy import *
import matplotlib.pyplot as plt

# Obtain path to cur directory
path = os.getcwd()

# Get of all file names from fluid folder ends with *.txt
files_x = []
files_y = []
for (dirpath, dirnames, filenames) in os.walk(path):
    for cur_file in filenames:
        if cur_file.endswith('.txt'):
            if cur_file.startswith('vx'):
                files_x.append(dirpath + "\\" + cur_file)
            elif cur_file.startswith('vy'):
                files_y.append(dirpath + "\\" + cur_file)


file_count = len(files_x)

# Loop throw all *.txt files in folder
for i in range(0, file_count):
    cur_file_x = files_x[i]
    cur_file_y = files_y[i]
    
    # Extract data from file
    data_x = np.array(np.loadtxt(cur_file_x))
    data_y = np.array(np.loadtxt(cur_file_y))
    
    # Get size of data
    nx = len(data_x[0])
    ny = len(data_x)
    
    # Arrange x,y axis values in appropriate intervals [0,nx] [0, ny]
    x = np.array([n for n in range(0, nx)])
    y = np.array([n for n in range(0, ny)])

    
    # Plot heatmap using extracted data
    (X,Y) = meshgrid(x,y)
    
    
    #c = plt.contourf(x,y,heatmap_data)
    # plt.contour(x, y, heatmap_data, 100, colors='k')
    strm = plt.streamplot(x, y, data_x, data_y, color= data_x, linewidth=2, cmap=plt.cm.autumn)
    plt.colorbar(strm.lines)
    plt.title('Stremline plot for velocity')
    
    lx = plt.xlabel("X axis")
    ly = plt.ylabel("Y axis")
    ax = plt.axis([0,nx,0,ny])

    # Obtain output file name and save it
    out_file_name = re.split(r'fluid_txt\\', cur_file_x)[1]
    plt.savefig('streamlines_' + out_file_name[:-4] + '.png', format = 'png')
    
    plt.clf()   
    print('streamlines_' + out_file_name + '.png complete.')
