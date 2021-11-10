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
from scipy.constants import G
from scipy.optimize import curve_fit
plt.rcParams["font.size"] =13
from scipy.interpolate import UnivariateSpline

#m12i
massSnap=[5543750000.0, 20971750000.0, 40641500000.0, 61843000000.0, 83163500000.0, 104176750000.0, 125208500000.0, 146756500000.0, 168737000000.0, 190738750000.0, 212401250000.0, 233703750000.0, 254553750000.0, 275254000000.0, 295184250000.0, 314619250000.0, 334300000000.0, 353782250000.0, 372899250000.0, 391263000000.0, 409546250000.0, 427567750000.0, 445733000000.0, 463717000000.0, 481010250000.0, 497740750000.0, 513780250000.0, 529309250000.0, 544975000000.0, 560099500000.0, 574898750000.0, 588946250000.0, 602522000000.0, 615285250000.0, 627393500000.0, 639123750000.0, 650461250000.0, 661348250000.0, 671958250000.0]
RSnap=[0.001794576923076923, 0.005383730769230769, 0.008972884615384615, 0.012562038461538461, 0.016151192307692307, 0.019740346153846153, 0.0233295, 0.026918653846153846, 0.03050780769230769, 0.034096961538461534, 0.03768611538461539, 0.041275269230769226, 0.04486442307692308, 0.04845357692307692, 0.05204273076923077, 0.05563188461538461, 0.059221038461538464, 0.0628101923076923, 0.06639934615384616, 0.0699885, 0.07357765384615385, 0.07716680769230769, 0.08075596153846154, 0.08434511538461538, 0.08793426923076923, 0.09152342307692307, 0.09511257692307692, 0.09870173076923076, 0.10229088461538462, 0.10588003846153846, 0.10946919230769231, 0.11305834615384615, 0.1166475, 0.12023665384615384, 0.1238258076923077, 0.12741496153846155, 0.13100411538461537, 0.13459326923076922, 0.13818242307692308]




if __name__ == "__main__":
    #address='*.h5'
    #address='/media/shahram/SD/Sample100Mpc/m12i/tags/rem/AllTags_161.h5'
    #address='/media/shahram/JB3/2021/AllTags/AllTags_262.h5'
    #AllTags_264.h5
    address='/media/shahram/JB3/2021/AllTagsPosFixed/rem/*.h5'
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
    NBins=40
    rBins=np.linspace(0,Rv,NBins+1)
    Rs=[0.]*NBins
    SigmaV=[0.]*NBins
    Density=[0.]*NBins
    Vc=[0.]*NBins
    MIn=[0.]*NBins
    #for i in range(0,NBins):
    #    Rs[i]=(Rbins[i]+Rbins[i+1])/2.
    rT=np.array(rT)
    VxT=np.array(VxT)
    VyT=np.array(VyT)
    VzT=np.array(VzT)
    SMT=np.array(SMT)
    mm=0.
    for i in range(0,NBins):
        Rin=rBins[i]
        Rout=rBins[i+1]
        Rs[i]=(rBins[i]+rBins[i+1])/2.
        VxBin=VxT[(rT>Rin) & (rT<Rout)]
        VyBin=VyT[(rT>Rin) & (rT<Rout)]
        VzBin=VzT[(rT>Rin) & (rT<Rout)]
        SMBin=SMT[(rT>Rin) & (rT<Rout)]
        SMBin[np.isnan(SMBin)]=0.
        SMBin[np.isinf(SMBin)]=0.
        #print("SMBin:")
        #print(SMBin[np.isinf(SMBin)])
        mm+=np.sum(SMBin,dtype=np.float64)
        #for j in range(0,len(SMBin)):
        #    mm+=SMBin[j]
        MIn[i]=mm
        dV=(4./3.)*3.1415*((Rout*1000.)**3.-(Rin*1000.)**3.)
        rho=SMBin/dV
        Density[i]=rho
        print("m:%g"%mm)
        #print(SMBin)
        #print("SM:%g"%np.nansum(SMBin,dtype=np.float64))
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
    #Now fitting V_sigma
    #
    SigmaVFit = UnivariateSpline(Rs, SigmaV)
    #
    #Now Jeans modeling
    M_enclosed=[0.]*NBins
    r=[0.]*NBins
    for i in range(0,NBins-1):
        #Dsigma=SigmaV[i+1]-SigmaV[i]
        Dsigma=Density[i+1]*SigmaVFit(Rs[i+1])-Density[i]*SigmaVFit(Rs[i])
        DR=Rs[i+1]-Rs[i]
        DR/=Density[i+1]
        r[i]=(Rs[i+1]+Rs[i])/2.
        M_enclosed[i]=np.abs((Dsigma/DR)*(r[i]**2.)/G)
    #plt.scatter(xT,zT , c=ageT,cmap = 'gist_earth', s =1, alpha =0.3) # gist_earth YlGn
    plt.plot(Rs,SigmaV, c='grey')
    plt.plot(Rs,SigmaVFit(Rs),c='black')
    plt.title("$\\sigma_v $")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('$\\sigma_v (Km s^{-1})$')
    print(Density)
    #plt.savefig('Age.png')
    fig3= plt.figure(3)
    plt.plot(Rs,np.log10(Density), c='black')
    #plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("$ Log(\\rho) -R$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$ \\rho $")
    #Vc
    fig4= plt.figure(4)
    plt.plot(Rs,Vc, c='black')
    #plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("$ V_{circ}-R $")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$ V_c $")
    #mass profile
    fig5= plt.figure(5)
    plt.plot(r,np.log10(M_enclosed), c='black', label='Jeans')
    plt.plot(RSnap,np.log10(massSnap), c='grey',linestyle=':',label='DM halo')
    plt.plot(Rs,np.log10(MIn), c='grey',linestyle='-.',label='SM')
    plt.title("$M_{enclosed}$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$Log(M_{enclosed})$")
    plt.legend()
    #
    plt.show()
