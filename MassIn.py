#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
#import yt
import numpy as np
from pygadgetreader import *
from yt.analysis_modules.halo_finding.api import *
from yt.analysis_modules.halo_analysis.api import *
from os import environ
environ['CFLAGS'] = "-I"+np.get_include()

import pyximport; pyximport.install()
#import particle_ops
import argparse


import tempfile
import shutil
import os
import sys

from scipy.spatial.transform import Rotation as R
from numpy import linalg as LA
from operator import mul
from functools import reduce
import matplotlib.pyplot as plt

import csv

plt.rcParams["font.size"] =12

    #how to run: python DEnsityProfile.py snapshot_file
    #example: $python DensityProfile.py /media/shahram/SD/Sample100Mpc/snap_500
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("snap", type=str)
    args = parser.parse_args()
    #snap = yt.load(args.snap)
    #ad = snap.all_data()
    #print(snap.field_list)
    #coordinatesDM = ad[("All","Coordinates")]
    p=readsnap(args.snap,"pos","dm")
    print(p[:,0].shape)
    #p=np.array(coordinatesDM)
    #print(p[:,1].shape)
    #print(p[1:])
    #m12i: 29.3575,31.0276,32.4926,6.37801e+11,0.139977,264,5.1374e+10,0,0,0,3.99904e+
    #m12b: 27.5708,29.1913,27.5166,8.01988e+11,0.15109,41,3.11153e+08,0,0,0,470366
    #m12f: 27.1756,33.4577,32.8438,9.1826e+11,0.158065
    #m12i
    #x0=29.3575
    #y0=31.0276
    #z0=32.4926
    #Rv=0.139977
    #m12b
    x0=27.5708
    y0=29.1913
    z0=27.5166
    Rv=0.15109
    #m12f
    #x0=27.1756
    #y0=33.4577
    #z0=32.8438
    #Rv=0.158065
    #x0=0
    #y0=0
    #z0=0
    #Rv=1000 #200kpc
    nbins=40
    px=p[:,0]
    py=p[:,1]
    pz=p[:,2]
    dx=px-x0
    dy=py-y0
    dz=pz-z0
    r2=dx*dx+dy*dy+dz*dz
    r=np.sqrt(r2)
    Rbins=np.linspace(0,Rv,nbins)
    rho=[0.0]*nbins
    NFW=[0.0]*nbins
    massIn=[0.0]*nbins
    Rs=[0.0]*nbins
    mass=0
    for i in range(0,nbins-1):
        Rin=Rbins[i]
        Rout=Rbins[i+1]
        R=(Rin+Rout)/2.
        Rs[i]=R
        Vin=(4./3.)*3.14*(Rin**3.)
        Vout=(4./3.)*3.14*(Rout**3.)
        r_in=r[(r>Rin) & (r<Rout)]
        rho[i]=(len(r_in)/(Vout-Vin))*1000
        mass+=len(r_in)*2.3726664349995558*100000.
        massIn[i]=mass
        NFW[i]=(1./((R/10.)*(1+R/10.)**2.))*25000
        #2000000 particles make a halo, say 2*10^12 M_sun  so each particle has a mass of 10^6 M_sun, 2000 is 200Kpc so r is 0.1 kpc
        #so rho*10^3 is M_sun/kpc^3
        #x=px[(r>Rin) & (r<Rout)]
        #y=py[r<Rv]
        #z=pz[r<Rv]
    # print(Rbins)
    # com_x=np.sum(x)/len(x)
    # com_y=np.sum(y)/len(y)
    # com_z=np.sum(z)/len(z)
    print("mass:")
    print(massIn)
    print("Rs:")
    print(Rs)
    fig1 = plt.figure(1)
    fig1.suptitle('Density Profile')
    ax1 = fig1.add_subplot(111)
    #ax1.title.set_text('b/a')
    ax1.set_xlabel('R [kpc]')
    ax1.set_ylabel('$Log \\rho [M_\\odot kpc^{-3}] $')
    ax1.plot(Rbins/4.,np.log10(rho),color='black',linestyle='-',label='Final Halo')
    ax1.plot(Rbins/4.,np.log10(NFW),color='black',linestyle=':',label='Initial NFW')
    ax1.legend(loc=1)
    plt.show()
    # print("COM (x,y,z):%g,%g,%g"%(com_x,com_y,com_z))
