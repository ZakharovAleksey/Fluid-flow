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
files = []
for (dirpath, dirnames, filenames) in os.walk(path):
    for cur_file in filenames:
        if cur_file.endswith('.txt'):
            files.append(dirpath + "\\" + cur_file)

# Loop throw all *.txt files in folder
for cur_file in files:
    # Extract data from file
    heatmap_data = np.loadtxt(cur_file)
    # Get size of data
    nx = len(heatmap_data[0])
    ny = len(heatmap_data)
    # Arrange x,y axis values in appropriate intervals [0,nx] [0, ny]
    x = [n for n in range(0, nx)]
    y = [n for n in range(0, ny)]

    
    # Plot heatmap using extracted data
    (X,Y) = meshgrid(x,y)
    
    
    #c = plt.contourf(x,y,heatmap_data)
    plt.contour(x, y, heatmap_data, 25, colors='k')
    
    lx = plt.xlabel("x")
    ly = plt.ylabel("y")
    ax = plt.axis([0,nx,0,ny])

    # Obtain output file name and save it
    out_file_name = re.split(r'fluid_txt\\', cur_file)[1]
    plt.savefig('contour_' + out_file_name[:-4] + '.png', format = 'png')
    
    plt.clf()   
    print('contour_' + out_file_name + '.png complete.')
