#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
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

#How to use: $python TagAnalysis.py HDf_tag_file
#example: python TagAnalysis.py StellarHalo.h5


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("TagFile", type=str)
    #f=h5.File("StellarHalo.h5","r")
    f=h5.File(TagFile,"r")
    args = parser.parse_args()
    #
    datasetNames = [n for n in f.keys()]
    for n in datasetNames:
        print(n)
    halo=f['FinalTag']
    age=halo['Age']
    StellarMass=halo['StellarMass']
    metallicity=halo['ZZ']
    print(halo.shape)
    x=halo['X']
    y=halo['Y']
    z=halo['Z']
    #
    #
    print(len(x))
    print(metallicity)
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
    plt.scatter(x,z , c=np.log10(metallicity),cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.title("metallicity")
    fig4=plt.figure(4)
    plt.scatter(x,z , c=StellarMass,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.title("StellarMass")
    fig5=plt.figure(5)
    plt.scatter(x,z , c=age,cmap = 'gist_earth', s =2, alpha =0.8)
    cbar = plt.colorbar()
    plt.title("age")
    plt.show()
