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
    #address='/media/shahram/ShahramWD1/AllTagsPosFixed/test/*.h5'
    #address='/media/shahram/ShahramWD1/m12i_SH/accretedsmooth/*.h5'
    #address='/media/shahram/JB3/m12i_SH/accretedsmooth/*.h5'
    #address='/media/shahram/JB3/m12f_SH/*.h5'
    address='/media/shahram/JB3/m12b_SH/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    #m12i
    #gx=29.3575
    #gy=31.0276
    #gz=32.4926
    #Rv=0.139977
    #m12f
    # gx=27.1756
    # gy=33.4577
    # gz=32.8438
    # Rv=0.158065
    #m12b
    gx=27.5708
    gy=29.1913
    gz=27.5166
    Rv=0.15109
    fig1= plt.figure(1)
    MetallicityT=[]
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
            f.close()
            print("finished finding particles in Rv")
            xT.extend(x)
            zT.extend(z)
            MetallicityT.extend(Metallicity)
            rBin=r[r<Rv]
            Rss.extend(rBin)
            #print(ID0)
            #print(float(Vx0[ID0==313488]))
            #if(len(ID0[ID0==313488])>0):
            #    print("yay!")
    #plt.scatter(xT,zT , c=np.log10(MetallicityT),cmap = 'gist_earth', s =2, alpha =0.8)
    #cmap = 'gist_earth'    'viridis'  gnuplot   YlGn
    #cbar = plt.colorbar()
    #plt.clim(-7, -2)
    #cbar.set_label('Metallicity')
    plt.scatter(gx,gz,c='r',marker='+',alpha=0.5,s=25)
    plt.title("metallicity")
    plt.xlabel('x (Mpc)')
    plt.ylabel('z (Mpc)')
    #plt.savefig('Metallicity.png')
    #plt.hist(np.log10(MetallicityT),linewidth=2, bins=100, log=False,cumulative=-1, histtype='step', alpha=0.9,color='blue',label='DMO')
    #plt.hist(np.log10(MetallicityT[MetallicityT!=0]),linewidth=2, bins=100, log=False, histtype='step', alpha=0.9,color='blue',label='Z')
    fig2= plt.figure(2)
    MetallicityT=np.array(MetallicityT)
    logzz=np.log10(MetallicityT/0.019)
    Rss=np.array(Rss)
    logzz=np.array(logzz)
    logzz[np.isnan(logzz)]=0.
    logzz[np.isinf(logzz)]=0.
    Rss2=Rss[(logzz > -2.5)& (logzz < -0.1)]
    logzz2=logzz[(logzz > -2.5)& (logzz < -0.1)]
    zHeat,xedge, yedge=np.histogram2d(Rss2,logzz2,bins=[50,50])#, weights=SM)
    #ext=[xedge[0],xedge[-1],-5,yedge[-1]]#[0,12,-5,0]#[xedge[0],xedge[-1],yedge[0],yedge[-1]]
    Xmesh, Ymesh = np.meshgrid(xedge, yedge)
    plt.title("$Metallicity-R$")
    plt.xlabel('$R$')
    plt.ylabel('$Metallicity$')
    #plt.imshow(VrHeat.T, origin='lower')
    plt.pcolormesh(Xmesh, Ymesh, zHeat.T,cmap='gist_earth')#,extent=ext)# cmap='gist_earth')
    cbar = plt.colorbar()
    #cross section cuts
    fig3=plt.figure(3)
    #ax1=fig3.add_subplot(111)
    #colors=ax1.pcolormesh(xT,zT,np.log10(MetallicityT), cmap='RdBu')#, vmin=np.min(np.log10(MetallicityT)), vmax=np.max(np.log10(MetallicityT)))
    #ax1.set_title('Metallicity')
    #ax1.axis([xT.min(), xT.max(), zT.min(),zT.max()])
    #fig3.colorbar(colors,ax=ax1)
    pCount,xHist, zHist=np.histogram2d(xT,zT,bins=[60,60])
    XHistmesh, ZHistmesh = np.meshgrid(xHist, zHist)
    print(pCount.shape)
    print(logzz.shape)
    logzHist=pCount*0.0 #just to create an empty array with the same shape
    print(zHist)
    xT=np.array(xT)
    zT=np.array(zT)
    logzz=np.array(logzz)
    zCutUp=0
    zCutDown=-3.0
    xRange=np.array(xT[(logzz > zCutDown)& (logzz < zCutUp)])
    zRange=np.array(zT[(logzz >zCutDown)& (logzz < zCutUp)])
    logzRange=np.array(logzz[(logzz > zCutDown)& (logzz < zCutUp)])
    for ii in range(0,len(xHist)-1):
        for jj in range(0,len(zHist)-1):
            xMin=xHist[ii]
            xMax=xHist[ii+1]
            zMin=zHist[jj]
            zMax=zHist[jj+1]
            xRangeNew=xRange[(xRange>xMin) & (xRange<xMax)]
            zRangeNew=zRange[(xRange>xMin) & (xRange<xMax)]
            logzRangeNew=logzRange[(xRange>xMin) & (xRange<xMax)]
            #print(len(xRangeNew),len(zRangeNew),len(logzRangeNew))
            #print(len(zMin))
            #print(len(xRangeNew[(zRangeNew>zMin)]))
            xRangeNew2=xRangeNew[(zRangeNew>zMin) & (zRangeNew<zMax)]
            zRangeNew2=zRangeNew[(zRangeNew>zMin) & (zRangeNew<zMax)]
            logzRangeNew2=logzRangeNew[(zRangeNew>zMin) & (zRangeNew<zMax)]
            #print(xHist[ii])
            #print(xHist[ii+1])
            #print(np.min([xHist[ii],xHist[ii+1]]))
            #logzRangeNew=logzRange[(xRange>xHist[ii]) & (xRange<xHist[ii+1]) & (zRange>zHist[jj]) & (zRange<zHist[jj+1])]
            #print(zRange[xRange>xHist[ii]])
            #print(xRange<xHist[ii+1])
            #zRangeNew=zRange[(xRange>xHist[ii]) & (xRange<xHist[ii+1])]
            #logzRangeNew=logzRangeNew[(zRangeNew<zHist[jj]) & (zRangeNew<zHist[jj+1])]
            logzHist[ii,jj]=np.mean(logzRangeNew2)
    #zz=logzz.reshape(len(zHist),len(xHist))
    plt.pcolormesh(XHistmesh, ZHistmesh, logzHist,cmap='jet')#,extent=ext)# cmap='gist_earth')
    cbar = plt.colorbar()
    plt.title("$Metallicity$")# (-1 > 0)$")
    plt.xlabel('$x [Mpc h^{-1}]$')
    plt.ylabel('$z[Mpc h^{-1}]$')
    plt.show()
