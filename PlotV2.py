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
#from scipy.constants import G
from scipy.optimize import minimize
from scipy.optimize import curve_fit
plt.rcParams["font.size"] =13
from scipy.interpolate import UnivariateSpline
import math
from scipy.spatial.transform import Rotation as R

def FitMass(variables,a,b):
	#a,b,c=params
	#R,z=variables
    R=variables
	#return G*((a*(R)+ (b*z**1.5))*kn(0,(c*R)/(2*R0))*np.exp((-e*z)/z0)+d)
    return a*np.log10(b*R)




G=6.6E-6
######################################
#mass unit is M_sun
#m12b
#massSnap=[10432139781.406046, 34359768978.44607, 63591019457.640594, 95161007241.81418, 126847255681.30275, 158054225500.92142, 189114090001.5701, 219339724983.67294, 248479153271.83298, 275933276591.2128, 302311395682.3204, 327857895187.9606, 351536394409.3256, 373782040370.5945, 395120141420.476, 415560425491.3537, 435510042143.47345, 455007903106.3693, 474255685027.01624, 493589357472.6101, 513603510651.76184, 534326853828.33496, 554744834834.7236, 574825422674.0554, 594199667715.6877, 612103334100.9073, 628781755539.0933, 644686213186.1823, 660182809473.0948, 675465391247.5704, 690570260306.0646, 705481519783.4628, 720641671969.8925, 735633839372.7241, 749986573171.3235, 763828709153.1108, 776967824070.2079, 789739650223.167, 801985930760.7737, 0.0]
#RSnap=[0.001937051282051282, 0.005811153846153846, 0.00968525641025641, 0.013559358974358974, 0.017433461538461537, 0.021307564102564104, 0.025181666666666665, 0.029055769230769232, 0.03292987179487179, 0.03680397435897435, 0.04067807692307693, 0.04455217948717949, 0.04842628205128205, 0.05230038461538461, 0.056174487179487184, 0.060048589743589745, 0.0639226923076923, 0.06779679487179487, 0.07167089743589744, 0.075545, 0.07941910256410256, 0.08329320512820512, 0.0871673076923077, 0.09104141025641026, 0.09491551282051282, 0.09878961538461538, 0.10266371794871795, 0.10653782051282051, 0.11041192307692307, 0.11428602564102563, 0.11816012820512821, 0.12203423076923077, 0.12590833333333334, 0.12978243589743588, 0.13365653846153847, 0.137530641025641, 0.1414047435897436, 0.14527884615384618, 0.1491529487179487, 0.0]
######################################
#m12f
massSnap=[1507829519.4422176, 12861987477.48909, 43178732850.695915, 80348925221.39896, 117968975148.17792, 155190654380.52045, 191671112618.56912, 227105699492.06998, 261362494747.23706, 295327926564.1862, 330219172824.0722, 365088353286.1127, 396503880485.3678, 425143150689.0299, 452112063722.45233, 478314130964.08295, 504046884784.8706, 528823216765.70996, 552138460755.8766, 574698959553.0698, 596131255460.4208, 616691359186.2659, 636640738571.7422, 656097315138.6691, 675204397939.7205, 693979544706.5155, 712821126133.4905, 731507772442.26, 749748357461.2496, 767322223212.0044, 784426775541.9161, 801249692365.9935, 817506253712.0365, 833399085293.594, 849164267420.9485, 864224055817.1776, 879391800536.1992, 894984964347.0164, 910680864614.4689, 0.0]

RSnap=[0.0020264743589743592, 0.006079423076923078, 0.010132371794871796, 0.014185320512820515, 0.01823826923076923, 0.022291217948717953, 0.02634416666666667, 0.03039711538461539, 0.034450064102564106, 0.03850301282051283, 0.04255596153846154, 0.046608910256410264, 0.05066185897435898, 0.0547148076923077, 0.058767756410256416, 0.06282070512820515, 0.06687365384615385, 0.07092660256410258, 0.07497955128205129, 0.0790325, 0.08308544871794873, 0.08713839743589745, 0.09119134615384616, 0.0952442948717949, 0.0992972435897436, 0.10335019230769232, 0.10740314102564104, 0.11145608974358975, 0.11550903846153848, 0.1195619871794872, 0.12361493589743591, 0.12766788461538464, 0.13172083333333334, 0.13577378205128207, 0.1398267307692308, 0.1438796794871795, 0.1479326282051282, 0.15198557692307696, 0.15603852564102566, 0.0]


######################################

#m12i
######################################
#massSnap=[5261387819.611515, 19903586923.280773, 38571489167.213776, 58693124135.87101, 78927698026.83423, 98870671212.93599, 118831202130.65675, 139281688667.20493, 160142646497.00803, 181023771991.5086, 201582926650.77975, 221800417343.41098, 241588455411.30728, 261234370759.7471, 280149504846.20703, 298594613711.89355, 317272955688.14056, 335762907949.4486, 353906213644.6032, 371334634942.89246, 388686656381.97473, 405790259645.3125, 423030291228.6628, 440098304495.4756, 456510750026.29803, 472389108342.60205, 487611661656.2722, 502349716483.91547, 517217556165.55316, 531571713564.0135, 545617187059.28033, 558949199757.5428, 571833490299.5209, 583946664250.1241, 595438199594.7574, 606570987774.4188, 617331030057.1418, 627663517848.278, 637733114198.416, 0.0]
#RSnap=[0.001794576923076923, 0.005383730769230769, 0.008972884615384615, 0.012562038461538461, 0.016151192307692307, 0.019740346153846153, 0.0233295, 0.026918653846153846, 0.03050780769230769, 0.034096961538461534, 0.03768611538461539, 0.041275269230769226, 0.04486442307692308, 0.04845357692307692, 0.05204273076923077, 0.05563188461538461, 0.059221038461538464, 0.0628101923076923, 0.06639934615384616, 0.0699885, 0.07357765384615385, 0.07716680769230769, 0.08075596153846154, 0.08434511538461538, 0.08793426923076923, 0.09152342307692307, 0.09511257692307692, 0.09870173076923076, 0.10229088461538462, 0.10588003846153846, 0.10946919230769231, 0.11305834615384615, 0.1166475, 0.12023665384615384, 0.1238258076923077, 0.12741496153846155, 0.13100411538461537, 0.13459326923076922, 0.13818242307692308]
#####################################



if __name__ == "__main__":
    #address='*.h5'
    #address='/media/shahram/SD/Sample100Mpc/m12i/tags/rem/AllTags_161.h5'
    #address='/media/shahram/JB3/2021/AllTags/AllTags_262.h5'
    #AllTags_264.h5
    #address='/media/shahram/JB3/2021/AllTagsPosFixed/test/*.h5'
    #address='/media/shahram/ShahramWD1/AllTagsPosFixed/rem/*.h5'
    address='/media/shahram/JB3/m12f_SH/*.h5'
    #address='/media/shahram/JB3/m12b_SH/*.h5'
    #address='/media/shahram/SD/Sample100Mpc/m12b/m12b_SH/*.h5'
    #address='/media/shahram/JB3/m12i_SH/accretedsmooth/*.h5'
    #AllTagsPosFixedPosFixed_194.h5
    #Halo info from Rockstar
    #m12i
    #m12i: 29.3575,31.0276,32.4926,6.37801e+11,0.139977,264,5.1374e+10,0,0,0,3.99904e+
    #6440 2389286 7.631e+11 7.631e+11 188.977402 149.599915 54.658009 143.317734 29.357544 31.029387 32.492146 -52.824707 72.188049 100.574295 1.37545e+10 7.13653e+10 -1.10587e+11 -7.5704e+15 0.00526233 0.000045 0.274646 -46.014793 73.618195 95.998154 0.079915 354301 8.420496e+11 6.376857e+11 4.718821e+11 2.158707e+11 4.269987 8.328306 0.004924 0.628998 0.506087 1.346174 11.681837 -3.791546 0.619734 0.468343 0.902458 8.931920 -2.634730 24.263170 25.247494 0.509286 8.279805e+11 4.877253e+11 73.982613 41996 -1 -1 3876024 10000000000.000000
    #m12b
    #m12b: 27.5708,29.1913,27.5166,8.01988e+11,0.15109,41,3.11153e+10,0,0,0,470366
    #4743 3092690 9.359e+11 9.359e+11 202.266769 175.132599 36.814861 167.140518 27.570570 29.191301 27.516590 -80.757416 -28.606602 -103.005913 3.08601e+10 -6.20163e+11 1.37631e+12 -1.12312e+16 0.0439077 0.000031 0.272087 -77.683327 -28.382364 -103.024040 0.084158 563571 1.013820e+12 8.018820e+11 6.181135e+11 3.439163e+11 7.166698 3.082309 0.039981 0.678771 0.597124 9.754769 2.045950 -2.318844 0.614595 0.531647 7.884814 1.441257 -2.065913 18.070126 17.765816 0.537297 1.067220e+12 4.676613e+11 72.334473 39743 -1 -1 4465312 10000000000.000000
    #m12f
    ##m12f: 27.1756,33.4577,32.8438,9.1826e+11,0.158065
    #10156 3370251 1.106e+12 1.106e+12 213.875168 180.465057 49.707188 176.831360 27.169409 33.463512 32.848351 -154.445251 167.958664 121.871223 -1.07745e+12 -2.29755e+12 -2.17838e+11 -1.31142e+16 0.0526631 0.000039 0.276490 -150.686554 162.742706 109.950363 0.081887 576036 1.234221e+12 9.089681e+11 7.032026e+11 3.860202e+11 10.361948 13.544039 0.051022 0.513396 0.446613 10.009308 -5.450640 9.133945 0.467729 0.431883 7.850304 -4.099245 7.586720 21.718937 21.017351 0.568074 1.221494e+12 6.061415e+11 77.199913 59149 -1 -1 5614381 10000000000.000000
    ######################################
    #m12i
    #gpos=np.array([29.3575,31.0276,32.4926])
    #gvel=np.array([-52.883121,72.168541,100.636299])
    #gJ=np.array([1.37545e+10,7.13653e+10,-1.10587e+11])
    #Rv=0.139977
    #m12b
    #gpos=np.array([27.5708,29.1913,27.5166])
    #gvel=np.array([-80.757416,-28.606602,-103.005913])
    #gJ=np.array([3.08601e+10,-6.20163e+11,1.37631e+12])
    #Rv=0.15109
    #m12f
    gpos=np.array([27.1756,33.4577,32.8438])
    gvel=np.array([-154.445251,167.958664,121.871223])
    gJ=np.array([-1.07745e+12,-2.29755e+12,-2.17838e+11])
    Rv=0.158065
    ## bins number
    NBins=10
    ##############
    gx=gpos[0]
    gy=gpos[1]
    gz=gpos[2]
    VxH=gvel[0]
    VyH=gvel[1]
    VzH=gvel[2]
    Jx=gJ[0]
    Jy=gJ[1]
    Jz=gJ[2]
    JNorm=np.sqrt(Jx**2.+Jy**2.+Jz**2.)
    JxyNorm=np.sqrt(Jx**2.+Jy**2.)
    rotationAxis= np.array([Jy/JxyNorm,-Jx/JxyNorm,0])
    rotationAngle=np.arccos(Jz/JNorm)
    rotationVector=rotationAngle*rotationAxis
    rotationJ=R.from_rotvec(rotationVector)
    #change them to a list when you're done with tests!
    fig1= plt.figure(1)
    VxT=[]
    VyT=[]
    VzT=[]
    rT=[]
    xT=[]
    yT=[]
    zT=[]
    SMT=[]
    for h5name in glob.glob(address):
        #print(h5name)
        with h5.File(h5name, "r") as f:
            # List all groups
            print("Keys: %s" % f.keys())
            f_key=list(f.keys())
            print("Read keys")
            a_group_key = f_key[0]
            age0 =np.array(f[a_group_key])
            #print("0")
            a_group_key = f_key[1]
            GID0 = np.array(f[a_group_key])
            #print("1")
            a_group_key = f_key[2]
            HID0 = np.array(f[a_group_key])
            #print("2")
            a_group_key = f_key[3]
            ID0 =np.array(f[a_group_key])
            #print("3")
            a_group_key = f_key[4]
            Metallicity0=np.array(f[a_group_key])
            #print("4")
            a_group_key = f_key[5]
            StellarMass0 =np.array(f[a_group_key])
            #print("5")
            a_group_key = f_key[6]
            Vx0 =np.array(f[a_group_key])
            #print("6")
            a_group_key = f_key[7]
            Vy0 =np.array(f[a_group_key])
            #print("7")
            a_group_key = f_key[8]
            Vz0=np.array(f[a_group_key])
            #print("8")
            a_group_key = f_key[9]
            x0 =np.array(f[a_group_key])
            #print("9")
            a_group_key = f_key[10]
            y0 =np.array(f[a_group_key])
            #print("10")
            a_group_key = f_key[11]
            z0 =np.array(f[a_group_key])
            #print("11")
            #print(z0[z0 !=0])
            #print("finished reading%s"%h5name)
            dx=x0-gx
            dy=y0-gy
            dz=z0-gz
            r=np.sqrt(dx*dx+dy*dy+dz*dz)#**0.5
            print("separation is done")
            #print("len age=%d, len r=%d"%(len(age0),len(r)))
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
            x-=gx
            y-=gy
            z-=gz
            rr=r[r<Rv]
            f.close()
            print("finished finding particles in Rv")
            xT.extend(x)
            yT.extend(y)
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
    #Wr better rotate these to disk coordinates
    for pp in range(0,len(xT)):
        pos=[xT[pp],yT[pp],zT[pp]]
        vel=[VxT[pp],VyT[pp],VzT[pp]]
        rotPos=rotationJ.apply(pos)
        rotVel=rotationJ.apply(vel)
        xT[pp]=rotPos[0]
        yT[pp]=rotPos[1]
        zT[pp]=rotPos[2]
        VxT[pp]=rotVel[0]
        VyT[pp]=rotVel[1]
        VzT[pp]=rotVel[2]
    plt.title("Vz-R")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('Vz $(Km s^{-1})$')
    #now binning
    fig2= plt.figure(2)
    #NBins=10
    rBins=np.linspace(0,Rv,NBins+1)
    Rs=[0.]*NBins
    SigmaV=[0.]*NBins
    Density=[0.]*NBins
    Vc=[0.]*NBins
    MIn=[0.]*NBins
    Vr2Mean=[0.]*NBins
    Vfi2Mean=[0.]*NBins
    Vtheta2Mean=[0.]*NBins
    beta=[0.]*NBins
    vr=[0.]*NBins
    vt=[0.]*NBins
    vf=[0.]*NBins
    Sigvr=[0.]*NBins
    Sigvt=[0.]*NBins
    Sigvf=[0.]*NBins
    #for i in range(0,NBins):
    #    Rs[i]=(Rbins[i]+Rbins[i+1])/2.
    rT=np.array(rT)
    xT=np.array(xT)
    yT=np.array(yT)
    zT=np.array(zT)
    VxT=np.array(VxT)
    VyT=np.array(VyT)
    VzT=np.array(VzT)
    SMT=np.array(SMT)#*2.2.3726664349995558E+005
    #mass resolution in Fire 2.3726664349995558E-005 and mass unit is E019 M_sun in Gadget, h is missing here
    #massConversion=(2.3726664349995558E+005)/0.702
    #SMT*=massConversion
    mm=0.
    VrAll=[]
    Rss=[]
    for i in range(0,NBins):
        Rin=rBins[i]
        Rout=rBins[i+1]
        Rs[i]=(rBins[i]+rBins[i+1])/2.
        xBin=xT[(rT>Rin) & (rT<Rout)]
        yBin=yT[(rT>Rin) & (rT<Rout)]
        zBin=zT[(rT>Rin) & (rT<Rout)]
        rBin=rT[(rT>Rin) & (rT<Rout)]
        VxBin=VxT[(rT>Rin) & (rT<Rout)]
        VyBin=VyT[(rT>Rin) & (rT<Rout)]
        VzBin=VzT[(rT>Rin) & (rT<Rout)]
        SMBin=SMT[(rT>Rin) & (rT<Rout)]
        SMBin[np.isnan(SMBin)]=0.
        SMBin[np.isinf(SMBin)]=0.
        #print("SMBin:")
        #print(SMBin[np.isinf(SMBin)])
        totMIn=np.sum(SMBin,dtype=np.float64) #mass unit is tha same as tag unit
        #for j in range(0,len(SMBin)):
        #    mm+=SMBin[j]
        mm+=totMIn
        MIn[i]=mm
        dV=(4./3.)*3.1415*((Rout*1000.)**3.-(Rin*1000.)**3.) #unit is kpc^3
        rho=totMIn/dV
        Density[i]=rho
        print("m:%g"%mm)
        #print(SMBin)
        #print("SM:%g"%np.nansum(SMBin,dtype=np.float64))
        Vbin=(VxBin**2.+VyBin**2.+VzBin**2.)**0.5
        #SigmaV[i]=np.std(Vbin)
        #1st calc
        # thetaBin=np.arctan2(VyBin,VxBin)
        # fiBin=np.arccos(VzBin/Vbin)
        # VrBin=VxBin*np.sin(fiBin)*np.cos(thetaBin)+VyBin*np.sin(fiBin)*np.sin(thetaBin)+VzBin*np.cos(fiBin)
        # VthetaBin=VxBin*(-np.sin(thetaBin))+VyBin*np.cos(thetaBin)
        # VfiBin=VxBin*np.cos(fiBin)*np.cos(thetaBin)+VyBin*np.cos(fiBin)*np.sin(thetaBin)-VzBin*np.sin(fiBin)
        #2nd calc
        fiBin=np.arctan2(yBin,xBin)
        thetaBin=np.arccos(zBin/rBin)
        VrBin=VxBin*np.cos(fiBin)*np.sin(thetaBin)+VyBin*np.sin(fiBin)*np.sin(thetaBin)+VzBin*np.cos(thetaBin)
        VthetaBin=VxBin*np.cos(fiBin)*np.cos(thetaBin)+VyBin*np.sin(fiBin)*np.cos(thetaBin)-VzBin*np.sin(thetaBin)
        VfiBin=VxBin*(-np.sin(fiBin))+VyBin*np.cos(fiBin)
        vr[i]=statistics.mean(VrBin)
        vt[i]=statistics.mean(VthetaBin)
        vf[i]=statistics.mean(VfiBin)
        Sigvr[i]=np.std(VrBin)
        Sigvt[i]=np.std(VthetaBin)
        Sigvf[i]=np.std(VfiBin)
        #
        #Vcirc=np.sqrt(VthetaBin**2.+VfiBin**2.)
        Vcirc=np.sqrt(vt[i]**2.+vf[i]**2.)
        SigmaV[i]=np.std(VrBin)
        Vr2Mean[i]=statistics.mean(VrBin**2.)
        Vfi2Mean[i]=statistics.mean(VfiBin**2.)
        Vtheta2Mean[i]=statistics.mean(VthetaBin**2.)
        beta[i]=1-(Vtheta2Mean[i]+Vfi2Mean[i])/(2.*Vr2Mean[i])
        VrAll.extend(VrBin)
        #rss=np.sqrt(xBin**2.+yBin**2.+zBin**2.)
        #print(rBin)
        Rss.extend(rBin) #Mpc
        #beta[i]/=2.
        #massBin=np.nansum(SMBin)
        #Density[i]=np.sum(rho)#massBin/dV
        Vc[i]=Vcirc#statistics.mean(Vcirc)
    #Now fitting V_sigma
    #
    SigmaVFit = UnivariateSpline(Rs, SigmaV)
    Vr2MeanFit= UnivariateSpline(Rs, Vr2Mean)
    DensityFit = UnivariateSpline(Rs, Density,s=2)
    # Log fit
    SigmaVFit2 = UnivariateSpline(Rs, SigmaV)# np.log(SigmaV))
    Vr2MeanFit20=UnivariateSpline(Rs,np.log10(Vr2Mean))
    Vr2MeanFit2=10.**(Vr2MeanFit20(Rs))
    DensityFit20 = UnivariateSpline(Rs, np.log10(Density),s=2)
    DensityFit2=10.**(DensityFit20(Rs))
    #
    Rss=np.array(Rss)
    VrAll=np.array(VrAll)
    #
    #Now Jeans modeling
    M_enclosed0=[0.]*NBins
    M_enclosedbeta=[0.]*NBins
    r=[0.]*NBins
    M_enclosed02=[0.]*NBins
    M_enclosedbeta2=[0.]*NBins
    dLnr2=np.gradient(Rs)
    print(dLnr2)
    dLnVr22=np.gradient(Vr2MeanFit2)
    print(dLnVr22)
    dLnRho2=np.gradient(DensityFit2)
    for i in range(0,NBins-1):
        #Dsigma=SigmaV[i+1]-SigmaV[i]
        #Dsigma=DensityFit(Rs[i+1])*SigmaVFit(Rs[i+1])-DensityFit(Rs[i])*SigmaVFit(Rs[i])
        #Dsigma=DensityFit(Rs[i+1])*Vr2MeanFit(Rs[i+1])-DensityFit(Rs[i])*Vr2MeanFit(Rs[i])
        #DR=Rs[i+1]-Rs[i]
        #DR/=DensityFit(Rs[i+1])
        r[i]=(Rs[i+1]+Rs[i])/2.
        dLnr=np.log(1000.*Rs[i+1])-np.log(1000.*Rs[i])
        dLnVr2=np.log(Vr2MeanFit(Rs[i+1]))-np.log(Vr2MeanFit(Rs[i]))
        dLnRho=np.log(DensityFit(Rs[i+1]))-np.log(DensityFit(Rs[i]))
        #M_enclosed[i]=np.abs((Dsigma/DR)*(r[i]**2.)/G)
        M_enclosed0[i]=np.abs(((dLnRho/dLnr)+(dLnVr2/dLnr))*(r[i]*1000.*Vr2MeanFit(Rs[i+1]))/G)
        #M_enclosed0[i]*=0.7
        M_enclosedbeta[i]=np.abs(((dLnRho/dLnr)+(dLnVr2/dLnr)+2*beta[i+1])*(r[i]*1000.*Vr2MeanFit(Rs[i+1]))/G)
        #M_enclosedbeta[i]*=0.7 #damn little h
        print(i)
        M_enclosed02[i]=np.abs(((dLnRho2[i]/dLnr2[i])+(dLnVr22[i]/dLnr2[i]))*(Rs[i]*Vr2MeanFit(Rs[i]))/G)
        M_enclosedbeta2[i]=np.abs(((dLnRho2[i]/dLnr2[i])+(dLnVr22[i]/dLnr2[i])+2*beta[i])*(Rs[i]*Vr2MeanFit(Rs[i]))/G)#np.abs(((dLnRho/dLnr)+(dLnVr2/dLnr)+2*beta[i+1])*(r[i]*Vr2MeanFit(Rs[i+1]))/G)
    # dLnr2=np.gradient(Rs)
    # dLnVr22=np.gradient(Vr2MeanFit2)
    # dLnRho2=np.gradient(DensityFit2)
    #M_enclosed02=np.abs(((dLnRho2/dLnr2)+(dLnVr22/dLnr2))*(Rs*Vr2MeanFit2(Rs))/G)#UnivariateSpline(Rs,np.abs(((dLnRho2/dLnr2)+(dLnVr22/dLnr2))*(Rs*Vr2MeanFit2(Rs))/G))
    #M_enclosedbeta2=np.abs(((dLnRho2/dLnr2)+(dLnVr22/dLnr2)+2*beta)*(Rs*Vr2MeanFit2(Rs))/G)#UnivariateSpline(Rs,np.abs(((dLnRho2/dLnr2)+(dLnVr22/dLnr2)+2*beta)*(Rs*Vr2MeanFit2(Rs))/G))
    #M_enclosedbeta2=np.abs(((dLnRho/dLnr)+(dLnVr2/dLnr)+2*beta[i+1])*(r[i]*Vr2MeanFit(Rs[i+1]))/G)
    #plt.scatter(xT,zT , c=ageT,cmap = 'gist_earth', s =1, alpha =0.3) # gist_earth YlGn
    #Spline fit
    M_enclosed0Fit = UnivariateSpline(Rs, M_enclosed0)
    M_enclosedbetaFit=UnivariateSpline(Rs, M_enclosedbeta)
    #non-spline fit
    initial_guess = [2.0, 2.0]
    data=r#[r,z]
    #parameters,param_covariance=curve_fit(FitForce,(R,z),Fz,initial_guess)
    #parameters,param_covariance=curve_fit(FitMass,(r),M_enclosedbeta,initial_guess)
    #a,b=parameters
    #M_enclosedbetaFit2=FitMass(data,a,b)
    #plot
    plt.plot(Rs,SigmaV, c='red', linestyle=':',label='Data')
    plt.plot(Rs,SigmaVFit(Rs),c='blue',label='Fit')
    plt.plot(Rs,SigmaVFit2(Rs),c='blue',linestyle='-.',label='Fit2')
    plt.title("$\\sigma_v $")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('$\\sigma_v (Km s^{-1})$')
    plt.legend()
    print(Density)
    #plt.savefig('Age.png')
    fig3= plt.figure(3)
    plt.plot(np.log10(Rs),np.log10(Density), c='red', linestyle=':',label='Data')
    plt.plot(np.log10(Rs),np.log10(DensityFit(Rs)), c='blue',linestyle='-.',label='Fit')
    #plt.plot(Rs,np.log10(np.exp(DensityFit2(Rs))), c='black',linestyle='-.',label='Fit2')
    #plt.scatter(rT,VzT, s =1,c='black', alpha =0.5) # gist_earth YlGn
    #cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #plt.scatter(gx,gz,c='r',marker='+',alpha=0.5)
    plt.title("$ Log(\\rho) -R$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$ Log( \\rho) $")
    plt.legend()
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
    r1 = r[:-1]
    #M_enclosed0Fit1 = M_enclosed0Fit[:-1]
    #M_enclosedbetaFit1=M_enclosedbetaFit[:-1]
    fig5= plt.figure(5)
    plt.plot(r,np.log10(M_enclosed0), c='black',linestyle='--', label='Jeans, $\\beta =0$')
    plt.plot(r,np.log10(M_enclosedbeta), c='black', label='Jeans, $\\beta$')
    plt.plot(r1,np.log10(M_enclosed0Fit(r1)), c='green',linestyle='--', label='Jeans, $\\beta =0$, fit')
    plt.plot(r1,np.log10(M_enclosedbetaFit(r1)), c='green', label='Jeans, $\\beta$, fit')
    plt.plot(RSnap,np.log10(massSnap), c='grey',linestyle=':',label='DM halo')
    plt.plot(Rs,np.log10(MIn), c='grey',linestyle='-.',label='SM')
    #plt.plot(r,np.log10(M_enclosed02), c='gray',linestyle='-.', label='Jeans2, $\\beta =0$')
    #plt.plot(r,np.log10(M_enclosedbeta2), c='gray', label='Jeans2, $\\beta$')
    plt.title("$M_{enclosed}$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$Log(M_{enclosed})$")
    plt.legend()
    #
    fig6= plt.figure(6)
    plt.plot(Rs,beta, c='black', label='$\\beta$')
    #plt.plot(RSnap,np.log10(massSnap), c='grey',linestyle=':',label='DM halo')
    #plt.plot(Rs,np.log10(MIn), c='grey',linestyle='-.',label='SM')
    plt.title("$\\beta$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$\\beta$")
    plt.legend()
    #
    fig7= plt.figure(7)
    plt.plot(Rs,np.abs(vr), c='blue', linestyle='-', label='$V_r$')
    plt.plot(Rs,np.abs(vt), c='blue', linestyle='-.',label='$V_{\\theta}$')
    plt.plot(Rs,np.abs(vf), c='blue',linestyle='--', label='$V_{\\phi}$')
    plt.plot(Rs,np.abs(Sigvr), c='red', linestyle='-', label='$\\sigma_{V_r}$')
    plt.plot(Rs,np.abs(Sigvt), c='red', linestyle='-.',label='$\\sigma_{V_{\\theta}}$')
    plt.plot(Rs,np.abs(Sigvf), c='red',linestyle='--', label='$\\sigma{V_{\\phi}}$')
    #plt.plot(RSnap,np.log10(massSnap), c='grey',linestyle=':',label='DM halo')
    #plt.plot(Rs,np.log10(MIn), c='grey',linestyle='-.',label='SM')
    plt.title("$V$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$V$")
    plt.legend()
    #
    fig8= plt.figure(8)
    #print(Rss)
    VrHeat,xedge, yedge=np.histogram2d(Rss,VrAll,bins=80)
    ext=[0,0.2,0,60]#[xedge[0],xedge[-1],yedge[0],yedge[-1]]
    Xmesh, Ymesh = np.meshgrid(xedge, yedge)
    plt.title("$V_r-R$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel('$V_r(Km s^{-1})$')
    #plt.imshow(VrHeat.T, origin='lower')
    plt.pcolormesh(Xmesh, Ymesh, VrHeat.T, cmap='gist_earth')
    cbar = plt.colorbar()
    #cbar.set_label('Age (Gyr)')
    #print(xedge)
    #print(yedge)
    #plt.pcolor()
    fig9= plt.figure(9)
    plt.plot(Rs,np.log10(np.abs(Vr2MeanFit(Rs))), c='black', linestyle='-', label='$V^2_r$')
    plt.title("$V$")
    plt.xlabel('d ($Mpc h^{-1}$)')
    plt.ylabel("$Log(V^2_r)$")
    plt.show()
