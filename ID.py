#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#this code is a simple code to check halo and subhalo IDs for given snapshot tags



import h5py as h5
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import math

#How to use: $python TagAnalysisGal.py HDf_tag_file file_flag
#example: python ID.py StellarHalo.h5 1
#This works for a single snapshot

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    parser.add_argument("Flag", type=int)
    args = parser.parse_args()
    f=h5.File(args.TagFile,"r")
    if args.Flag:
        halo=f['FinalTag']
    else:
        halo=f['FullTag']
    x0=halo['X']
    y0=halo['Y']
    z0=halo['Z']
    Mv0=halo['Mvir']
    Hindex0=halo['HaloIndex']
    SubHIndex=halo['SubhaloIndex']
    GIndex=halo['GalIndex']
    BE0=halo['BindingEnergy']
    StellarMass=halo['StellarMass']
    print("Halo Index:")
    print(Hindex0[Hindex0!=0])
    print("Sub:")
    print(SubHIndex[SubHIndex!=0])
    print("StellarMass:")
    print(StellarMass)
