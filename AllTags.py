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
#How to use: $python Accretion.py HDf_tag_file halo_catalog galaxy_file start_snap end_sap
#example: python AllTags.py ~/Desktop/Research/3_Tagging/TagAnalysis/StellarHaloSubselection.h5 ~/Desktop/Research/3_Tagging/TagAnalysis/FinalTag717Stars /media/shahram/SD/Sample100Mpc/717/halos_c3z717.ascii galsSM.ascii 262 264
#This works for a single halo/galaxy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    parser.add_argument("TagFiles", type=str)
    parser.add_argument("HaloFile", type=str)
    parser.add_argument("GalFile", type=str)
    parser.add_argument("Si", type=int)
    parser.add_argument("Sf", type=int)
    args = parser.parse_args()
    #f=h5.File("StellarHalo.h5","r")
    #Halo
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
    NBins=6
    #
    ph=pnumh[(MvH>LowerMass) & (MvH<UpperMass)]
    Idh=IdH[(MvH>LowerMass) & (MvH<UpperMass)]
    Mvh=MvH[(MvH>LowerMass) & (MvH<UpperMass)]
    xh=xH[(MvH>LowerMass) & (MvH<UpperMass)]
    yh=yH[(MvH>LowerMass) & (MvH<UpperMass)]
    zh=zH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh=RvH[(MvH>LowerMass) & (MvH<UpperMass)]
    Rvh/=1000 # convert from kpc to Mpc
    #
    # Galaxy
    Gals=np.genfromtxt(args.GalFile, delimiter = ',')
    Gx0=np.array(Gals[:,0])
    Gy0=np.array(Gals[:,1])
    Gz0=np.array(Gals[:,2])
    GMv0=np.array(Gals[:,3])
    GRv0=np.array(Gals[:,4])
    GRd0=np.array(Gals[:,5])
    GSM0=np.array(Gals[:,6])
    #GSM0=np.array(Gals[:,6])
    Gx=Gx0[(GMv0>LowerMass) & (GMv0<UpperMass)]#48.81
    Gy=Gy0[(GMv0>LowerMass) & (GMv0<UpperMass)]#44.62#
    Gz=Gz0[(GMv0>LowerMass) & (GMv0<UpperMass)]#49.6777#
    GMv=GMv0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GRv=GRv0[(GMv0>LowerMass) & (GMv0<UpperMass)]#*2#*10000
    GRd=GRd0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    GSM=GSM0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    print("Galaxy:%g-%g-%g"%(Gx,Gy,Gz))
    #GSM=GSM0[(GMv0>LowerMass) & (GMv0<UpperMass)]
    if len(GMv)>1 or len(Mvh)>1:
        print("I got more than one halo/galaxy. I'd better stop")
        exit(1)
    if len(GMv)==0:
        print("I got no halo/galaxy. I'd better stop")
        exit(1)
    #
    fT=h5.File(args.TagFile,"r")
    haloT=fT['FinalTag'] # for full tag
    #halo=f['FullTag'] # for individual tags
    ID0T=haloT['PID']
    x0T=haloT['X']
    y0T=haloT['Y']
    z0T=haloT['Z']
    BE0T=haloT['BindingEnergy']
    IDT=ID0T[BE0T!=0]
    xT=x0T[BE0T!=0]
    yT=y0T[BE0T!=0]
    zT=z0T[BE0T!=0]
    dx2T=(Gx-xT)**2.
    dy2T=(Gy-yT)**2.
    dz2T=(Gz-zT)**2.
    rT=np.sqrt(dx2T+dy2T+dz2T)
    targetID=IDT[rT<GRv]
    print("targetID:")
    print(targetID)
    L=len(targetID)
    print(L)
    S=args.Sf-args.Si
    pxAll=[0.0]*L*S
    pyAll=[0.0]*L*S
    pzAll=[0.0]*L*S
    pAgeAll=[0.0]*L*S
    pStellarMassAll=[0.0]*L*S
    pMetallicityAll=[0.0]*L*S
    rAll=[0.0]*L*S
    t=0
    for i in range(args.Sf,args.Si,-1):
        num="/tag_%03d.h5"%i
        tFile=args.TagFiles+str(num)
        print(tFile)
        f=h5.File(tFile,"r")
        datasetNames = [n for n in f.keys()]
        for n in datasetNames:
            print(n)
        halo=f['FullTag']
        #halo=f['FinalTag'] # for full tag
        #halo=f['FullTag'] # for individual tags
        age0=halo['Age']
        StellarMass0=halo['StellarMass']
        metallicity0=halo['ZZ']
        #print(halo.shape)
        x0=halo['X']
        y0=halo['Y']
        z0=halo['Z']
        Mv0=halo['Mvir']
        Hindex0=halo['HaloIndex']
        BE0=halo['BindingEnergy']
        TreeIndex0=halo['TreeIndex']
        infallMvir0=halo['infallMvir']
        ID0=halo['PID']
        age=age0[BE0!=0]
        StellarMass=StellarMass0[BE0!=0]#*(1.0e10)
        metallicity=metallicity0[BE0!=0]#/0.0134
        x=x0[BE0!=0]
        y=y0[BE0!=0]
        z=z0[BE0!=0]
        Mv=Mv0[BE0!=0]
        Hindex=Hindex0[BE0!=0]
        TreeIndex=TreeIndex0[BE0!=0]
        #UTree = set(TreeIndex)
        infallMvir=infallMvir0[BE0!=0]
        ID=ID0[BE0!=0]
        #print("TreeIndex:%d out of %d"%(len(UTree),len(TreeIndex0)))
        print("ID in snap:")
        print(ID)
        print(len(ID))
        print("Snap:%d-Rv=%g"%(i,GRv))
        #pxtemp=[0.0]*len(targetID)
        #pytemp=[0.0]*len(targetID)
        #pztemp=[0.0]*len(targetID)
        #pAgetemp=[0.0]*len(targetID)
        #pStellarMasstemp=[0.0]*len(targetID)
        #pMetallicitytemp=[0.0]*len(targetID)
        #pTreeIndextemp=[0.0]*len(targetID)
        #rtemp=[0.0]*len(targetID)
        for i in range(0,len(targetID)):
            id=targetID[i]
            #for j in range(0,len(ID)):
            #id=targetID[i]
                #if ID[j]==id:
            if len(ID[ID==id])>0:
                print("got the id:%d"%id)
                #pxtemp[i]=x[ID==id]
                #pytemp[i]=y[ID==id]
                #pztemp[i]=z[ID==id]
                #pAgetemp[i]=age[ID==id]
                #pStellarMasstemp[i]=StellarMass[ID==id]
                #pMetallicitytemp[i]=metallicity[ID==id]
                #pTreeIndextemp[i]=TreeIndex[ID==id]
                pxAll[t]=x[ID==id]
                pyAll[t]=y[ID==id]
                pzAll[t]=z[ID==id]
                pAgeAll[t]=age[ID==id]
                pStellarMassAll[t]=StellarMass[ID==id]
                pMetallicityAll[t]=metallicity[ID==id]
                t+=1
        #now we have those particles in either cases
        #px=pxtemp[pxtemp!=0]
        #py=pytemp[pxtemp!=0]
        #pz=pztemp[pxtemp!=0]
        #pAge=pAgetemp[pxtemp!=0]
        #pStellarMass=pStellarMasstemp[pxtemp!=0]
        #pMetallicity=pMetallicitytemp[pxtemp!=0]
        #pTreeIndex=pTreeIndextemp[pxtemp!=0]
        #pxAll.append(px)
        #pyAll.append(py)
        #pzAll.append(pz)
        #pAgeAll.append(pAge)
        #pStellarMassAll.append(pStellarMass)
        #pMetallicityAll.append(pMetallicity)
        #rAll.append(r)
        #now prepare to plot
    # now this is a (sf-si)*N array and we have to flatten them to 1D arrays.
    px=pxAll[pxAll !=0]#np.array(pxAll).ravel()
    py=pyAll[pxAll !=0]#np.array(pyAll).ravel()
    pz=pzAll[pxAll !=0]#np.array(pzAll).ravel()
    pAge=pAgeAll[pxAll !=0]#np.array(pAgeAll).ravel()
    pStellarMass=pStellarMassAll[pxAll !=0]#np.array(pStellarMassAll).ravel()
    pMetallicity=pMetallicityAll[pxAll !=0]#np.array(pMetallicityAll).ravel()
    rAll=rAll[pxAll !=0]#np.array(rAll).ravel()
    pMetallicitylog=pMetallicity#[pMetallicityFinal !=0]#np.log10(pMetallicityFinal[pMetallicityFinal !=0])
    dx2=(Gx-px)**2.
    dy2=(Gy-py)**2.
    dz2=(Gz-pz)**2.
    rALl=np.sqrt(dx2+dy2+dz2)
    #print(pMetallicityFinal)
    #rsquared=(Gx-pxAll)**2.+(Gy-pyAll)**2.+(Gz-pzAll)**2.
    #pR=rAllFinal#np.sqrt(np.array(rsquared))*1000.
    #pR=pR[pMetallicityAll !=0]
    r2=rAll#pR[pMetallicityFinal !=0]*1000
    #pMetallicitylog=pMetallicitylog[pMetallicitylog>-3]
    if px:
        fig0=plt.figure(0)
        ax01=fig0.add_subplot(221)
        #ax01.plot(np.log10(Rs),np.log10(Rho))
        ax01.set_xlabel("$log(R(kpc))$")
        ax01.set_ylabel("$log(\\rho) [M_\\odot /kpc^{-3}]$")
        ax02=fig0.add_subplot(222)
        ax02.hist(pAge,linewidth=2, bins=10,weights=pStellarMass, log=True,cumulative=False, histtype='step', alpha=0.9,color='blue',label='age')
        #ax02.hist(pAgeAll,linewidth=2, bins=10, log=False,cumulative=False, histtype='step', alpha=0.9,color='blue',label='age')
        ax02.set_xlabel("Age")
        ax03=fig0.add_subplot(223)
        ax03.hist(pMetallicitylog,linewidth=2, bins=10,weights=pStellarMass, log=True,cumulative=False, histtype='step', alpha=0.9,color='blue',label='metallicity')
        #ax03.hist(pMetallicitylog,linewidth=2, bins=10, log=True,cumulative=False, histtype='step', alpha=0.9,color='blue',label='metallicity')
        ax03.set_xlabel("Metallicity$(Log(Z/Z_{\\odot}))$")
        ax04=fig0.add_subplot(224)
        ax04.scatter(r2,pMetallicitylog,s=1,c='black')
        #for i in range(0,len(Idh)):
        #metalicity-halo mass dependence
        #metalicity of the halo is the average metalicity
        print(GMv,GRv,GRd,GSM)
        print(np.sum(pStellarMass))
        plt.show()