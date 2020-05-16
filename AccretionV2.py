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
#How to use: $python AccretionV2.py HDf_tag_file FirstTagged galaxies_file
#example: python TagAnalysis.py StellarHalo.h5 FirstTagged.h5 gals.ascii
#This works for a single halo/galaxy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    parser.add_argument("FirstTag", type=str)
    parser.add_argument("GalDirectory", type=str)
    args = parser.parse_args()
    #f=h5.File("StellarHalo.h5","r")
    TagFF=h5.File(args.TagFile,"r")
    TagIF=h5.File(args.FirstTag,"r")
    #extract halox in a specific mass range, MWish for instance
    LowerMass=1.0e12
    UpperMass=1.3e12
    NBins=6
    #
    # Galaxies
    Gals=np.genfromtxt(args.GalFile, delimiter = ',')
    Gx0=np.array(Gals[:,0])
    Gy0=np.array(Gals[:,1])
    Gz0=np.array(Gals[:,2])
    GMv0=np.array(Gals[:,3])
    GRv0=np.array(Gals[:,4])
    GRd0=np.array(Gals[:,5])
    GSnap=np.array(Gals[:,6])
    datasetNames = [n for n in TagFF.keys()]
    for n in datasetNames:
        print(n)
    #Tagged Particles
    ####################################
    #Final Tag
    TagsF=TagFF['FinalTag'] # for full tag
    #halo=f['FullTag'] # for individual tags
    IDF0=TagsF['PID']
    ageF0=TagsF['Age']
    StellarMassF0=TagsF['StellarMass']
    metallicityF0=TagsF['ZZ']
    print(TagsF.shape)
    xF0=TagsF['X']
    yF0=TagsF['Y']
    zF0=TagsF['Z']
    MvF0=TagsF['Mvir']
    HindexF0=TagsF['HaloIndex']
    BEF0=TagsF['BindingEnergy']
    TreeIndexF0=TagsF['TreeIndex']
    infallMvirF0=TagsF['LastMajorMerger']
    SnapF0=TagsF['Snap']
    #
    IDF=IDF0[BEF0!=0]
    ageF=ageF0[BEF0!=0]
    StellarMassF=StellarMassF0[BEF0!=0]*(1.0e10)
    metallicityF=metallicityF0[BEF0!=0]/0.0134
    xF=xF0[BEF0!=0]
    yF=yF0[BEF0!=0]
    zF=zF0[BEF0!=0]
    MvF=MvF0[BEF0!=0]
    HindexF=HindexF0[BEF0!=0]
    TreeIndexF=TreeIndexF0[BEF0!=0]
    UTree = set(TreeIndexF)
    infallMvirF=infallMvirF0[BEF0!=0]
    #print("TreeIndex:%d out of %d"%(len(UTree),len(TreeIndex0)))
    #print(UTree)
    ##Extract particles for this specific halo/galaxy
    #halo
    #dx2=(xh-x)**2.
    #dy2=(yh-y)**2.
    #dz2=(zh-z)**2.
    #the main galaxy
    #################
    #First tagged
    TagsI=TagIF['FirstTagged'] # for full tag
    #halo=f['FullTag'] # for individual tags
    IDI0=TagsI['PID']
    ageI0=TagsI['Age']
    StellarMassI0=TagsI['StellarMass']
    metallicityI0=TagsI['ZZ']
    print(TagsI.shape)
    xI0=TagsI['X']
    yI0=TagsI['Y']
    zI0=TagsI['Z']
    MvI0=TagsI['Mvir']
    HindexI0=TagsI['HaloIndex']
    BEI0=TagsI['BindingEnergy']
    TreeIndexI0=TagsI['TreeIndex']
    infallMvirI0=TagsI['LastMajorMerger']
    SnapI0=TagsI['Snap']
    #
    IDI=IDI0[BEI0!=0]
    ageI=ageI0[BEI0!=0]
    StellarMassI=StellarMassI0[BEI0!=0]#*(1.0e10) it is already converted in PtagPP
    metallicityI=metallicityI0[BEI0!=0]/0.0134
    xI=xI0[BEI0!=0]
    yI=yI0[BEI0!=0]
    zI=zI0[BEI0!=0]
    MvF=MvF0[BEI0!=0]
    HindexI=HindexI0[BEI0!=0]
    TreeIndexI=TreeIndexI0[BEI0!=0]
    UTree3 = set(TreeIndexI)
    infallMvirI=infallMvirI0[BEI0!=0]
    SnapI=SnapI0[BEI0!=0]
    ##########################
    # The main galaxy
    S=SnapF0[0]
    GMx0=Gx0[GSnap==S]
    GMx0=Gy0[GSnap==S]
    GMz0=Gz0[GSnap==S]
    GMMv0=GMv0[GSnap==S]
    GMRv0=GRv0[GSnap==S]
    GMRd0=GRd0[GSnap==S]
    Gx=GMx0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    Gy=GMy0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    Gz=GMz0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GMv=GMMv0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GRv=GMRv0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GRd=GMRd0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    dx2=(Gx-xF)**2.
    dy2=(Gy-yF)**2.
    dz2=(Gz-zF)**2.
    r=np.sqrt(dx2+dy2+dz2)
    #now extract tagged particles within this halos virial radius
    ###########
    # Extract particles in this galaxy
    pID=IDF[r<GRv]
    px=xF[r<GRv]
    py=yF[r<GRv]
    pz=zF[r<GRv]
    pAge=ageF[r<GRv]
    pStellarMass=StellarMassF[r<GRv]
    pMetallicity=metallicityF[r<GRv]
    pTreeIndex=TreeIndexF[r<GRv]
    #
    UTree2 = set(pTreeIndexF)
    pinfallMvir=infallMvirF[r<GRv]
    #print("TreeIndex2:%d out of %d"%(len(UTree2),len(pTreeIndex)))
    #print(UTree2)
    #print("infallMvir:")
    #print(pinfallMvir[pTreeIndex==0])
    #Rbins=np.linspace(0,Rvh,NBins+1)
    #
    # Now we can find these particles in their first tagged positions
    ppx=[0.0]*len(pID)
    ppy=[0.0]*len(pID)
    ppz=[0.0]*len(pID)
    ppSnap=[0.0]*len(pID)
    for i in range(0,len(pID)):
        ppx[i]=xI[IDI==pID[i]]
        ppy[i]=yI[IDI==pID[i]]
        ppz[i]=zI[IDI==pID[i]]
        ppSnap=SnapI[IDI==pID[i]]
        Gx=Gx0[GSnap==ppSnap]
        Gy=Gy0[GSnap==ppSnap]
        Gz=Gz0[GSnap==ppSnap]
    #Any kind of analysis after we extract these points
    Rbins=np.linspace(0,GRv,NBins+1)
    print("Bins:")
    print(Rbins)
    Rs=[0]*NBins
    Rho=[0]*NBins
    Z=[0]*NBins
    for i in range(0,NBins):
        Rin=Rbins[i]
        Rout=Rbins[i+1]
        Rs[i]=(Rbins[i]+Rbins[i+1])*500. #(Rbins[i]+Rbins[i+1])/2. x 1000
        rbin=r[(r>Rin) & (r<Rout)]
        v=(4./3.)*np.pi*(Rout**3.-Rin**3.)*1.0e9 #r_Mpc to r_kpc (r^3)
        print(v)
        #Rho[i]=len(rbin)/v # all p have the same mass but don't forget to convert the units
        Rho[i]=np.sum(StellarMass[(r>Rin) & (r<Rout)])/v
        print(Rho[i])
    #fig0=plt.figure(0)
    #ax01=fig0.add_subplot(221)

    # metalicity at 30 kpc
        metalBin=pMetallicity#[(r>0.0) & (r<0.1)]
        Z[i]=np.sum(metalBin)#/len(metalBin)
    #ax02=fig0.add_subplot(222)
    print(pMetallicity)
    print("z30:")
    print(Z)
    #if not(math.isnan(met)):
        #ax02.plot(np.log10(Mvh[i]),np.log10(met))
    #print(met)
    ########
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
    #
    fig0=plt.figure(0)
    ax01=fig0.add_subplot(221)
    ax01.plot(np.log10(Rs),np.log10(Rho))
    ax01.set_xlabel("$log(R(kpc))$")
    ax01.set_ylabel("$log(\\rho) [M_\\odot /kpc^{-3}]$")
    ax02=fig0.add_subplot(222)
    ax02.hist(pAge,linewidth=2, bins=10, log=False,cumulative=False, histtype='step', alpha=0.9,color='blue',label='age')
    ax02.set_xlabel("Age")
    ax03=fig0.add_subplot(223)
    ax03.hist(pMetallicity,linewidth=2, bins=10, log=False,cumulative=False, histtype='step', alpha=0.9,color='blue',label='metallicity')
    ax03.set_xlabel("Metallicity$(Z/Z_{\\odot})$")
    ax04=fig0.add_subplot(224)
    ax04.plot(Rs,Z)
    #for i in range(0,len(Idh)):

    #
    #metalicity-halo mass dependence
    #metalicity of the halo is the average metalicity
    #

    print(GMv,GRv,GRd)
    #    #
    #
    fig1 = plt.figure(figsize=plt.figaspect(1))
    ax = fig1.add_subplot(111, projection='3d')
    ax.scatter(px,py,pz,c='black',alpha=0.8,marker='.',s=1)
    #
    ax.set_xlabel('X (Mpc)')
    ax.set_ylabel('Y (Mpc)')
    ax.set_zlabel('Z (Mpc)')
    #
    fig2 = plt.figure(2,figsize=plt.figaspect(1))
    ax2 = fig2.add_subplot(111)#, projection='3d')
    ax2.plot(px[pTreeIndex!=0],pz[pTreeIndex!=0],'k.', markersize=1)
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
    fig3= plt.figure(3,figsize=plt.figaspect(1))
    #ax3=fig3.add_subplot(111)
    #ax3.contour(X,Z,mass_new)
    #viridis = cm.get_cmap('viridis', 256)#np.max(mass_new))
    #psm=ax3.pcolormesh([x_new,z_new],cmap=viridis, rasterized=True)
    #fig3.colorbar(psm,ax=ax3)
    #plot cmap = 'RdPu'
    plt.scatter(px,pz , c=np.log10(pMetallicity),cmap = 'gist_earth', s =2, alpha =0.9)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("metallicity $log(Z/Z_{\\odot})$")
    fig4=plt.figure(4,figsize=plt.figaspect(1))
    plt.scatter(px,pz , c=pStellarMass,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("StellarMass ($M_{\odot}$)")
    fig5=plt.figure(5,figsize=plt.figaspect(1))
    plt.scatter(px,pz , c=pAge,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("age (Gyr)")
    # TreeIndex
    fig6=plt.figure(6,figsize=plt.figaspect(1))
    plt.scatter(px,pz , c=pTreeIndex,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.scatter(Gx,Gz,c='r',marker='+',alpha=0.4)
    plt.title("Tree Index")
    #Halo plots
    #


    plt.show()
