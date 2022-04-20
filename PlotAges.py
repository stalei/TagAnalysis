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
    #address='/media/shahram/JB3/2021/AllTags/AllTags_262.h5'
    #AllTags_264.h5
    #address='/media/shahram/JB3/2021/AllTagsPosFixed/rem/*.h5'
    #address='/media/shahram/ShahramWD1/AllTagsPosFixed/insitu/*.h5'
    address='/media/shahram/JB3/m12i_SH2/*.h5'
    #address='/media/shahram/JB3/m12i_SH/accretedsmooth/*.h5'
    #address='/media/shahram/SD/Sample100Mpc/m12b/m12b_SH/*.h5'
    #address='/media/shahram/JB3/m12b_SH/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    #m12i: 29.3575,31.0276,32.4926,6.37801e+11,0.139977,264,5.1374e+10,0,0,0,3.99904e+
    #m12b: 27.5708,29.1913,27.5166,8.01988e+11,0.15109,41,3.11153e+08,0,0,0,470366
    #m12f: 27.1756,33.4577,32.8438,9.1826e+11,0.158065
    #m12i
    gx=29.3575
    gy=31.0276
    gz=32.4926
    Rv=0.139977
    #m12b
    #gx=27.5708
    #gy=29.1913
    #gz=27.5166
    #Rv=0.15109
    #m12f
    #gx=27.1756
    #gy=33.4577
    #gz=32.8438
    #Rv=0.158065
    fig1= plt.figure(1)
    ageT=[]
    xT=[]
    zT=[]
    Rss=[]
    SM=[]
    for h5name in glob.glob(address):
        print(h5name)
        with h5.File(h5name, "r") as f:
            # List all groups
            print("Keys: %s" % f.keys())
            f_key=list(f.keys())
            print("Read keys")
            a_group_key = f_key[0]
            age0 =np.array(f[a_group_key])
            print("Age")
            a_group_key = f_key[1]
            GID0 = np.array(f[a_group_key])
            print("GID")
            a_group_key = f_key[2]
            HID0 = np.array(f[a_group_key])
            print("HID")
            a_group_key = f_key[3]
            ID0 =np.array(f[a_group_key])
            print("ID")
            a_group_key = f_key[4]
            Metallicity0=np.array(f[a_group_key])
            print("Metallicity")
            a_group_key = f_key[5]
            StellarMass0 =np.array(f[a_group_key])
            print("Stellar Mass")
            a_group_key = f_key[6]
            Vx0 =np.array(f[a_group_key])
            print("Vx")
            a_group_key = f_key[7]
            Vy0 =np.array(f[a_group_key])
            print("Vy")
            a_group_key = f_key[8]
            Vz0=np.array(f[a_group_key])
            print("Vz")
            a_group_key = f_key[9]
            x0 =np.array(f[a_group_key])
            print("X")
            a_group_key = f_key[10]
            y0 =np.array(f[a_group_key])
            print("Y")
            a_group_key = f_key[11]
            z0 =np.array(f[a_group_key])
            print("Z")
            print("finished reading")
            dx=x0-gx
            dy=y0-gy
            dz=z0-gz
            r=(dx*dx+dy*dy+dz*dz)**0.5
            print("separation is done")
            #print("len age=%d, len r=%d"%(len(age0),len(r)))
            print(r<Rv)
            age=age0[r<Rv]
            Metallicity=Metallicity0[r<Rv]
            StellarMass=StellarMass0[r<Rv]
            x=x0[r<Rv]
            y=y0[r<Rv]
            z=z0[r<Rv]
            rBin=r[r<Rv]
            f.close()
            print("finished finding particles in Rv")
            xT.extend(x)
            zT.extend(z)
            ageT.extend(age)
            Rss.extend(rBin)
            SM.extend(StellarMass)
            #print(ID0)
            #print(float(Vx0[ID0==313488]))
            #if(len(ID0[ID0==313488])>0):
            #    print("yay!")
    plt.scatter(xT,zT , c=ageT,cmap = 'gist_earth', s =1, alpha =0.3) # gist_earth YlGn
    cbar = plt.colorbar()
    cbar.set_label('Age (Gyr)')
    plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("Age")
    plt.xlabel('x (Mpc)')
    plt.ylabel('z (Mpc)')
    #plt.savefig('Age.png')
    fig2= plt.figure(2)
    #print(Rss)
    Rss=np.array(Rss)
    ageT=np.array(ageT)
    Rss2=Rss#[(ageT > 3.5) & (ageT < 4)]
    ageT2=ageT#[(ageT > 3.5) & (ageT < 4)]
    ageHeat,xedge, yedge=np.histogram2d(Rss2,ageT2,bins=[40,40])#, weights=SM) DOUBLE BINS
    #ext=[xedge[0],xedge[-1],-5,yedge[-1]]
    Xmesh, Ymesh = np.meshgrid(xedge, yedge)
    plt.title("$Age (weighted)-R$")
    plt.xlabel('d [$Mpc h^{-1}$]')
    plt.ylabel('$Age[Gyr]$')
    #plt.imshow(VrHeat.T, origin='lower')
    plt.pcolormesh(Xmesh, Ymesh, ageHeat.T,cmap='gist_earth')#,extent=ext)# cmap='gist_earth')
    cbar = plt.colorbar()
    plt.show()
