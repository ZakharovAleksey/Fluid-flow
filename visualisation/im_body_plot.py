# This script allows plotting of heatmap profile for chosen macroscopic
# variable with figure of immersed in fluid boundary.

import os # To get access for folders
import re # To use regular expressions
import numpy as np

# To plot haetmap and bodys boundary
from numpy import *
import matplotlib.pyplot as plt

# Obtain path to cur directory (contains fluid and body data)
path = os.getcwd()
# Set path to fluid folder
path_body = path + "\\fluid_txt\\"


# Get of all file names from fluid folder
fluid_files = []
for (dirpath, dirnames, filenames) in os.walk(path_body):
    for file in filenames:
        fluid_files.append(dirpath + file)

# Number of files in fluid folder (the same in body folder)
file_num = len(fluid_files)

for i in range (0, file_num):
    
    # Prepare data to heatmap plot
    heatmap_data = np.loadtxt(fluid_files[i])
    nx = len(heatmap_data[0])
    ny = len(heatmap_data)
    x = [n for n in range(0, nx)]
    y = [n for n in range(0, ny)]
    
    # Prepare data for body form plot
    body_file_name = re.split(r'fluid_txt', str(fluid_files[i]))[0]
    
    temp = re.split(r'fluid_txt', str(fluid_files[i]))[1]
    current_time = re.findall(r'_t\d{1,4}.txt', temp)[0]
    # Obtain name of appropriate body boundary

    body_x = []
    body_y = []
    
    for i in range(0,1):
        bfn = body_file_name
        bfn += "body_form_txt\\body_form" + str(i) + str(current_time)
        # Construct body graph
        body_data = np.loadtxt(bfn)
        
        for node in body_data:
            body_x.append(node[0])
            body_y.append(node[1])

    # Plot heatmap with immersed boundary
    (X,Y) = meshgrid(x,y)

    c = plt.contourf(x,y,heatmap_data)#,linspace(min(input),max(input),11))
    plt.plot(body_x, body_y, 'kx')
    plt.plot(body_x, body_y, 'k--')
    
    b = plt.colorbar(c, orientation='vertical')
    lx = plt.xlabel("x")
    ly = plt.ylabel("y")
    ax = plt.axis([0,nx,0,ny])
    #plt.show()
    plt.savefig('profile' + str(current_time[:-4])+ '.png', format = 'png')
    plt.clf()   # clear plot area
    print('profile' + str(current_time[:-4])+ '.png')
