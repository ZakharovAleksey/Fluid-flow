# This script allows plotting of contour levels for chosen macroscopic
# variable.

import os # To get access for folders
import re # To use regular expressions
import numpy as np

# To plot haetmap and bodys boundary
from numpy import *
import matplotlib.pyplot as plt

from matplotlib.ticker import NullFormatter  # useful for `logit` scale

def get_measure_points_number(file_name):
    with open(file_name, 'r') as f:
        first_line = f.readline()
        point_ids = re.split(r'\t', first_line)
        # remove '\n' and there are 2 columns on each point : 't' 'physVal'
        return (len(point_ids) - 1) / 2


# Obtain path to cur directory
path = os.getcwd()

# Get of all file names from fluid folder ends with *.txt
files = []
for (dirpath, dirnames, filenames) in os.walk(path):
    for cur_file in filenames:
        if cur_file.endswith('measure.txt'):
            files.append(dirpath + "\\" + cur_file)

for cur_file in files:
    # Obtain physical Value name
    value_name = re.split(r'\\spectrums\\', cur_file)[1]
    value_name = re.split(r'_measure', value_name)[0]
        
    # Calculate number of points in which we make measurements
    measue_point_number = get_measure_points_number(cur_file)
    # Extract data from file
    data = np.genfromtxt(cur_file, names = True)

    # Find maximum value of time, including tag time and simulation time
    max_pos_time = 0
    for point_id in range(0, measue_point_number):
        time = data['t' + str(point_id)]
        if time[-1] > max_pos_time:
            max_pos_time = int(time[-1] + 1)
            
    # Array for result spectrum 
    dep = np.zeros(max_pos_time)

    for point_id in range(0, measue_point_number):
        # Extract data
        time = data['t' + str(point_id)]
        signal = data['val' + str(point_id)]

        for t in range(0, len(time)):
            dep[time[t]] += signal[t]
        
        # Fourier FFT for real values
        f = np.fft.rfft(signal)
        # First value is main frequency - remove it from the spectrum because we want to find out bifurcation frequencies
        if value_name == 'rho':
            f = f[1:]
        # Find appropriate frequencies
        time_step = 1.0
        freqs = np.fft.fftfreq(f.size, time_step)
	# Shift to frequency
        xf = np.fft.fftshift(freqs)
        y_plot = np.fft.fftshift(f)

        # Plot single point measurements
        fig, ax = plt.subplots(2, 1)
        plt.gca().yaxis.set_minor_formatter(NullFormatter())
        plt.subplots_adjust(hspace=0.7)

        legend = ''
        if value_name == 'rho':
            ax[0].set_title('Dependency of density from time.')
            ax[0].legend(shadow=True)
            ax[1].set_title('Density spectrum.')
            ax[1].legend(shadow=True)
            legend = 'Density'
        else:
            ax[0].set_title('Dependency of velocity from time.')
            ax[0].legend(shadow=True)
            ax[1].set_title('Velocity spectrum.')
            ax[1].legend(shadow=True)
            legend = 'Velocity'
            
        # Plotting the signal
        ax[0].plot(time,signal, 'b', label = legend )
        ax[0].set_xlabel('Time')
        ax[0].grid()
        if value_name == 'rho':
            ax[0].set_ylabel('Density amplitude')
        else:
            ax[0].set_ylabel('Velocity amplitude')
            
        # Plotting the spectrum
        ax[1].plot(xf,y_plot.real,'r', label = legend)
        ax[1].set_xlabel('Freq (Hz)')
        ax[1].grid()
        if value_name == 'rho':
            ax[1].set_ylabel('Density amplitude')
        else:
            ax[1].set_ylabel('Velocity amplitude')
        ax[1].set_xlim([-0.005, 0.01])

        out_file_name = 'point' + str(point_id) + str(value_name) + '_spectrum.png'
        plt.savefig(out_file_name)
        print(out_file_name + ' written!')
        
    # Fourier FFT for real values
    f = np.fft.rfft(dep)
    # First value is main frequency - remove it from the spectrum because we want to find out bifurcation frequencies
    # Remove main components
    leng = len(f) / 2
    f = f[leng:]
    # Find appropriate frequencies
    time_step = 1.0
    freqs = np.fft.fftfreq(f.size, time_step)
    # Shift to frequency
    xf = np.fft.fftshift(freqs)
    y_plot = np.fft.fftshift(f)

    # Plot total spectrum
    fig, ax = plt.subplots(2, 1)
    # Format the minor tick labels of the y-axis into empty strings with
    # `NullFormatter`, to avoid cumbering the axis with too many labels.
    plt.gca().yaxis.set_minor_formatter(NullFormatter())
    plt.subplots_adjust(hspace=0.7)
    if value_name == 'rho':
        ax[0].set_title('Dependency of density from time.')
        ax[1].set_title('Density spectrum.')
    else:
        ax[0].set_title('Dependency of velocity from time.')
        ax[1].set_title('Velocity spectrum.')
    # Plotting the signal
    dep = dep[5000:35000]
    ax[0].plot([i + 5000 for i in range(0, len(dep))], dep, 'b', label = legend )
    ax[0].set_xlabel('Time')
    ax[0].grid()
    if value_name == 'rho':
        ax[0].set_ylabel('Total density amplitude')
    else:
        ax[0].set_ylabel('Velocity amplitude')
            
    # Plotting the spectrum
    ax[1].plot(xf,y_plot.real,'r', label = legend )
    ax[1].set_xlabel('Freq (Hz)')
    ax[1].grid()
    if value_name == 'rho':
        ax[1].set_ylabel('Total ensity amplitude')
    else:
        ax[1].set_ylabel('Velocity amplitude')
    out_file_name = 'total_' + str(value_name) + '_spectrum.png'
    plt.savefig(out_file_name)
    print(out_file_name + ' written!')
    plt.clf()
