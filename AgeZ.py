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
    #address='/media/shahram/ShahramWD1/AllTagsPosFixed/rem/*.h5'
    address='/media/shahram/ShahramWD1/m12i_SH/accretedsmooth/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    gx=29.3575
    gy=31.0276
    gz=32.4926
    Rv=0.139977
    fig1= plt.figure(1)
    SMT=[]
    MetallicityT=[]
    ageT=[]
    xT=[]
    zT=[]
    Rss=[]
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
            rBin=r[r<Rv]
            f.close()
            print("finished finding particles in Rv")
            xT.extend(x)
            zT.extend(z)
            ageT.extend(age)
            MetallicityT.extend(Metallicity)
            SMT.extend(StellarMass)
            #Rss.extend(rBin)
            #print(ID0)
            #print(float(Vx0[ID0==313488]))
            #if(len(ID0[ID0==313488])>0):
            #    print("yay!")
    #plt.scatter(ageT,MetallicityT , c=SMT,cmap = 'YlGn', s =1, alpha =0.2)
    MetallicityT=np.array(MetallicityT)
    logzz=np.log10(MetallicityT/0.019)
    plt.scatter(ageT,logzz , s =1, alpha =0.9)
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.4)
    plt.title("Age Metallicity")
    plt.xlabel('Age')
    plt.ylabel('Metallicity')
    #plt.savefig('Age.png')
    fig2= plt.figure(2)
    #print(Rss)
    ageT=np.array(ageT)
    logzz=np.array(logzz)
    logzz[np.isnan(logzz)]=0.
    logzz[np.isinf(logzz)]=0.
    ageT2=ageT[logzz > -4]
    logzz=logzz[logzz>-4]
    zHeat,xedge, yedge=np.histogram2d(ageT2,logzz,bins=[50,50])#, weights=SM)
    ext=[xedge[0],xedge[-1],-5,yedge[-1]]#[0,12,-5,0]#[xedge[0],xedge[-1],yedge[0],yedge[-1]]
    Xmesh, Ymesh = np.meshgrid(xedge, yedge)
    plt.title("$Age-Metallicity$")
    plt.xlabel('$Age[Gyr]$')
    plt.ylabel('$Z$')
    #plt.imshow(VrHeat.T, origin='lower')
    plt.pcolormesh(Xmesh, Ymesh, zHeat.T,cmap='gist_earth')#,extent=ext)# cmap='gist_earth')
    cbar = plt.colorbar()
    #now just a cross section.
    plt.show()
