#  © Shahram Talei @ 2021 The University of Alabama - All rights reserved.
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

if __name__ == "__main__":
    #address='*.h5'
    #address='/media/shahram/SD/Sample100Mpc/m12i/tags/rem/AllTags_161.h5'
    #address='/media/shahram/JB3/2021/AllTags/AllTags_170.h5'
    #AllTags_264.h5
    address='/media/shahram/JB3/2021/AllTagsPosFixed/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    gx=29.3575
    gy=31.0276
    gz=32.4926
    Rv=0.139977
    fig1= plt.figure(1)
    SMT=[]
    xT=[]
    zT=[]
    for h5name in glob.glob(address):
        print(h5name)
        with h5.File(h5name, "r") as f:
            # List all groups
            print("Keys: %s" % f.keys())
            f_key=list(f.keys())
            print("Read keys")
            a_group_key = f_key[0]
            age0 =np.array(f[a_group_key])
            print("0")
            a_group_key = f_key[1]
            GID0 = np.array(f[a_group_key])
            print("1")
            a_group_key = f_key[2]
            HID0 = np.array(f[a_group_key])
            print("2")
            a_group_key = f_key[3]
            ID0 =np.array(f[a_group_key])
            print("3")
            a_group_key = f_key[4]
            Metallicity0=np.array(f[a_group_key])
            print("4")
            a_group_key = f_key[5]
            StellarMass0 =np.array(f[a_group_key])
            print("5")
            a_group_key = f_key[6]
            Vx0 =np.array(f[a_group_key])
            print("6")
            a_group_key = f_key[7]
            Vy0 =np.array(f[a_group_key])
            print("7")
            a_group_key = f_key[8]
            Vz0=np.array(f[a_group_key])
            print("8")
            a_group_key = f_key[9]
            x0 =np.array(f[a_group_key])
            print("9")
            a_group_key = f_key[10]
            y0 =np.array(f[a_group_key])
            print("10")
            a_group_key = f_key[11]
            z0 =np.array(f[a_group_key])
            print("11")
            print("finished reading")
            dx=x0-gx
            dy=y0-gy
            dz=z0-gz
            r=(dx*dx+dy*dy+dz*dz)**0.5
            print("separation is done")
            print("len age=%d, len r=%d"%(len(age0),len(r)))
            print(r<Rv)
            age=age0[r<Rv]
            Metallicity=Metallicity0[r<Rv]
            StellarMass=StellarMass0[r<Rv]
            x=x0[r<Rv]
            y=y0[r<Rv]
            z=z0[r<Rv]
            f.close()
            print("finished finding particles in Rv, %d"%len(x))
            xT.extend(x)
            zT.extend(z)
            SMT.extend(StellarMass)
            #print(ID0)
            #print(float(Vx0[ID0==313488]))
            #if(len(ID0[ID0==313488])>0):
            #    print("yay!")
    plt.scatter(xT,zT , c=np.log10(SMT),cmap = 'YlGn', s =2, alpha =0.8)
    #np.log10(StellarMass)
    cbar = plt.colorbar()
    plt.clim(0.2, 7.5)
    cbar.set_label('Stellar Mass')
    plt.scatter(gx,gz,c='r',marker='+',alpha=0.4)
    plt.title("Log10(SM) $M_\odot$")
    plt.xlabel('x (Mpc)')
    plt.ylabel('z (Mpc)')
    #plt.savefig('StellarMass.png')
    plt.show()
