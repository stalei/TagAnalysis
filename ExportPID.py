#  © Shahram Talei @ 2022 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import division
import h5py as h5
import numpy as np
from os import environ
import os
environ['CFLAGS'] = "-I"+np.get_include()
import argparse
import glob

from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams["font.size"] =12


PList='/media/shahram/SD/Sample100Mpc/m12i/z0/halos_m12iC_z0_264.particles'
plist=np.genfromtxt(PList, skip_header=18,comments='#')
if plist is None:
    print ("Error, sorry, I couldn't read the halo binary file.!")
    sys.exit(1)
pids=np.array(plist[:,9])
#print(pids,end=',')
