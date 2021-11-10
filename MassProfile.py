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
import statistics
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
    address='/media/shahram/JB3/2021/AllTagsPosFixed/test/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    gx=29.3575
    gy=31.0276
    gz=32.4926
    Rv=0.139977
    #m12i -52.883121 72.168541 100.636299
    VxH=-52.883121
    VyH=72.168541
    VzH=100.636299
    fig1= plt.figure(1)
    VxT=[]
    VyT=[]
    VzT=[]
    rT=[]
    xT=[]
    zT=[]
    SMT=[]
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
            Vx=Vx0[r<Rv]
            Vy=Vy0[r<Rv]
            Vz=Vz0[r<Rv]
            #find velocities with respect to the halo
            Vx-=VxH
            Vy-=VyH
            Vz-=VzH
            rr=r[r<Rv]
            f.close()
            print("finished finding particles in Rv")
            xT.extend(x)
            zT.extend(z)
            rT.extend(rr)
            VxT.extend(Vx)
            VyT.extend(Vy)
            VzT.extend(Vz)
            SMT.extend(StellarMass)
            #ageT.extend(age)
            #print(ID0)
            #print(float(Vx0[ID0==313488]))
            #if(len(ID0[ID0==313488])>0):
            #    print("yay!")
    plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("Vz-R")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('Vz $(Km s^{-1})$')
    #now binning
    fig2= plt.figure(2)
    NBins=80
    rBins=np.linspace(0,Rv,NBins+1)
    Rs=[0.]*NBins
    SigmaV=[0.]*NBins
    Density=[0.]*NBins
    Vc=[0.]*NBins
    #for i in range(0,NBins):
    #    Rs[i]=(Rbins[i]+Rbins[i+1])/2.
    rT=np.array(rT)
    VxT=np.array(VxT)
    VyT=np.array(VyT)
    VzT=np.array(VzT)
    SMT=np.array(SMT)
    for i in range(1,NBins+1):
        Rin=rBins[i-1]
        Rout=rBins[i]
        Rs[i]=(Rout+Rin)/2.#(rBins[i]+rBins[i+1])/2.
        VxBin=VxT[(rT>Rin) & (rT<Rout)]
        VyBin=VyT[(rT>Rin) & (rT<Rout)]
        VzBin=VzT[(rT>Rin) & (rT<Rout)]
        SMBin=SMT[(rT>Rin) & (rT<Rout)]
        dV=(4./3.)*3.1415*((Rout*1000.)**3.-(Rin*1000.)**3.)
        rho=SMBin/dV
        print("dV:%g"%dV)
        print(SMBin)
        print("SM:%g"%np.nansum(SMBin))
        Vbin=(VxBin**2.+VyBin**2.+VzBin**2.)**0.5
        #SigmaV[i]=np.std(Vbin)
        tethaBin=np.arctan2(VyBin,VxBin)
        fiBin=np.arccos(VzBin/Vbin)
        VrBin=VxBin*np.sin(fiBin)*np.cos(tethaBin)+VyBin*np.sin(fiBin)*np.sin(tethaBin)+VzBin*np.cos(fiBin)
        VtethaBin=VxBin*(-np.sin(tethaBin))+VyBin*np.cos(tethaBin)
        VfiBin=VxBin*np.cos(fiBin)*np.cos(tethaBin)+VyBin*np.cos(fiBin)*np.sin(tethaBin)-VzBin*np.sin(fiBin)
        Vcirc=np.sqrt(VtethaBin**2.+VfiBin**2.)
        SigmaV[i]=np.std(VrBin)
        #massBin=np.nansum(SMBin)
        Density[i]=np.sum(rho)#massBin/dV
        Vc[i]=statistics.mean(Vcirc)
    #plt.scatter(xT,zT , c=ageT,cmap = 'gist_earth', s =1, alpha =0.3) # gist_earth YlGn
    plt.plot(Rs,SigmaV, c='black')
    plt.title("$\\sigma_v $")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('$\\sigma_v (Km s^{-1})$')
    print(Density)
    #plt.savefig('Age.png')
    fig3= plt.figure(3)
    plt.plot(Rs,Density, c='black')
    #plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("$ \\rho $-R")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$ \\rho $")
    #Vc
    fig4= plt.figure(4)
    plt.plot(Rs,Vc, c='black')
    #plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("$ V_{circ} $-R")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$ V_c $")
    plt.show()
