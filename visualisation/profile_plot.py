# This script allows plotting of profile for chosen macroscopic
# variable

import numpy as np
from os import listdir
import matplotlib.pyplot as plt

# create list of files in current directory
fileList = listdir(".")
# choose from all of them only *.txt files
txtFileList = filter(lambda x: x.endswith('.txt'), fileList)

for currentFile in txtFileList:
    print(currentFile)
    fileName = currentFile

    value = str(fileName[:-16])

    plt.plot(*np.loadtxt(fileName,unpack=True), linewidth=2.0)
    plt.title(value + " dependency")
    plt.xlabel(u'String Length')
    plt.ylabel(value)
    plt.grid()
    pngfileName = str(fileName[:-4])
    plt.savefig(pngfileName + '.png', format = 'png')
    plt.clf()
