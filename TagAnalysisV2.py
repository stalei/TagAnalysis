#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


import h5py as h5
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import math
#import csv
#How to use: $python TagAnalysisV2.py HDf_tag_file halo_catalog galaxy_file
#example: python TagAnalysis.py StellarHalo.h5 halos_0.0.ascii gal.csv
#This code works for multiple halos.

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    parser.add_argument("HaloFile", type=str)
    parser.add_argument("GalFile", type=str)
    args = parser.parse_args()
    #f=h5.File("StellarHalo.h5","r")
    f=h5.File(args.TagFile,"r")
    halos=np.genfromtxt(args.HaloFile, skip_header=18)
    pnumh=np.array(halos[:,1])
    MvH=np.array(halos[:,2])
    RvH=np.array(halos[:,4])# in kpc
    xH=np.array(halos[:,8])
    yH=np.array(halos[:,9])
    zH=np.array(halos[:,10])
    IdH=np.array(halos[:,0])
    #extract halox in a specific mass range, MWish for instance
    LowerMass=1.0e12
    UpperMass=1.3e12
    ph=pnumh[(MvH>LowerMass) & (MvH<UpperMass)]
    Idh=IdH[(MvH>LowerMass) & (MvH<UpperMass)]
    Mvh=MvH[(MvH>LowerMass) & (MvH<UpperMass)]
    xh=xH[(MvH>LowerMass) & (MvH<UpperMass)]
    yh=yH[(MvH>LowerMass) & (MvH<UpperMass)]
    zh=zH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh=RvH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh/=1000 # convert from kpc to Mpc
    #
    datasetNames = [n for n in f.keys()]
    for n in datasetNames:
        print(n)
    #halo=f['FinalTag'] # for full tag
    #halo=f['FirstTagged']
    halo=f['FullTag'] # for individual tags
    age0=halo['Age']
    print(age0)
    StellarMass0=halo['StellarMass']
    metallicity0=halo['ZZ']
    #print(halo.shape)
    x0=halo['X']
    y0=halo['Y']
    z0=halo['Z']
    Mv0=halo['Mvir']
    Hindex0=halo['HaloIndex']
    BE0=halo['BindingEnergy']
    GalI0=halo['GalIndex']
    GalNo0=halo['GalNo']
    #Gs=GalI0[(x0<48.8) & (x0>48.7) & (z0<49.7) & (z0>49.65)]
    #Gs=GalI0[(x0<48.8) & (x0>48.7)]
    #print(Gs)
    ########
    age=age0#[BE0!=0]
    StellarMass=StellarMass0#*(1.0e10)#[BE0!=0]*(1.0e10)
    metallicity=metallicity0#[BE0!=0]
    x=x0#[BE0!=0]
    y=y0#[BE0!=0]
    z=z0#[BE0!=0]
    Mv=Mv0#[BE0!=0]
    Hindex=Hindex0#[BE0!=0]
    ###########
    #
    # get totals for different halos
    #min=np.min(Hindex)
    #max=np.max(Hindex)
    #print(min,max)
    print(len(x))
    print(len(StellarMass[StellarMass !=0]))
    #min max didn't work so let's find another way to get the total properities
    #checking power law distribution
    #
    NBins=6
    fig0=plt.figure(0)
    ax01=fig0.add_subplot(221)
    ax02=fig0.add_subplot(222)
    for i in range(0,len(Idh)):
        dx2=(xh[i]-x)**2.
        dy2=(yh[i]-y)**2.
        dz2=(zh[i]-z)**2.
        r=np.sqrt(dx2+dy2+dz2)
        #now extract tagged particles within this halos virial radius
        px=x[r<Rvh[i]]
        py=y[r<Rvh[i]]
        pz=z[r<Rvh[i]]
        pAge=age[r<Rvh[i]]
        pStellarMass=StellarMass[r<Rvh[i]]
        pMetallicity=metallicity[r<Rvh[i]]
        Rbins=np.linspace(0,Rvh[i],NBins+1)
        print(Rbins)
        Rs=[0]*NBins
        Rho=[0]*NBins
        for i in range(0,NBins):
            Rin=Rbins[i]
            Rout=Rbins[i+1]
            Rs[i]=(Rbins[i]+Rbins[i+1])/2.
            rbin=r[(r>Rin) & (r<Rout)]
            v=(4./3.)*np.pi*(Rout**3.-Rin**3.)
            print(v)
            #Rho[i]=len(rbin)/v # all p have the same mass but don't forget to convert the units
            Rho[i]=np.sum(StellarMass[(r>Rin) & (r<Rout)])/v
            print(Rho[i])
        #fig0=plt.figure(0)
        #ax01=fig0.add_subplot(221)
        ax01.plot(np.log10(Rho),Rs)
        # metalicity at 30 kpc
        metalBin=metallicity[(r>0.01) & (r<0.05)]
        met=np.sum(metalBin)/len(metalBin)
        #ax02=fig0.add_subplot(222)
        #if not(math.isnan(met)):
            #ax02.plot(np.log10(Mvh[i]),np.log10(met))
        print(met)
    #
    #metalicity-halo mass dependence
    #metalicity of the halo is the average metalicity
    #
    # Galaxy plots
    Gals=np.genfromtxt(args.GalFile, delimiter = ',')
    Gx0=np.array(Gals[:,0])
    Gy0=np.array(Gals[:,1])
    Gz0=np.array(Gals[:,2])
    GMv0=np.array(Gals[:,3])
    GRv0=np.array(Gals[:,4])
    GRd0=np.array(Gals[:,5])
    Gx=Gx0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    Gy=Gy0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    Gz=Gz0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GMv=GMv0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GRv=GRv0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GRd=GRd0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    print(GMv,GRv,GRd)
    #    #
    #
    fig = plt.figure(figsize=plt.figaspect(1))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,c='black',alpha=0.8,marker='.',s=1)
    #
    ax.set_xlabel('X (Mpc)')
    ax.set_ylabel('Y (Mpc)')
    ax.set_zlabel('Z (Mpc)')
    #
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)#, projection='3d')
    ax2.plot(x,z,'k.', markersize=1)
    #fig.set_size_inches(14,8)
    ax2.set_xlabel('X (Mpc)')
    ax2.set_ylabel('Z (Mpc)')
    #ax.set_zlabel('Z (kpc)')
    # just pick a small area to test color-map
    #x_2=x[x>17]
    #z_2=z[x>17]
    #mass_2=mass[x>17]
    #x_new=x_2[x_2<19]
    #z_new=z_2[x_2<19]
    #mass_new=mass_2[x_2<19]
    #mass_new=mass_new.reshape(len(x_new),len(z_new))
    #X, Z =np.meshgrid(x_new,z_new)
    #fig3, ax3=plt.subplots()
    fig3= plt.figure(3)
    #ax3=fig3.add_subplot(111)
    #ax3.contour(X,Z,mass_new)
    #viridis = cm.get_cmap('viridis', 256)#np.max(mass_new))
    #psm=ax3.pcolormesh([x_new,z_new],cmap=viridis, rasterized=True)
    #fig3.colorbar(psm,ax=ax3)
    #plot cmap = 'RdPu'
    plt.scatter(x,z , c=np.log10(metallicity),cmap = 'gist_earth', s =2, alpha =0.9)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("metallicity")
    fig4=plt.figure(4)
    plt.scatter(x,z , c=StellarMass,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("StellarMass ($M_{\odot}$)")
    fig5=plt.figure(5)
    plt.scatter(x,y , c=age,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gy,c='r',marker='+',alpha=0.4)
    plt.title("age (Gyr)")
    #Halo plots
    #
    print("Galaxy pos(Gal/Rockstar):")
    print(Gx,Gy,Gz)
    print(xh,yh,zh)


    plt.show()
