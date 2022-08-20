import numpy as np
import scipy as sp 
import math
import sys
import os
import scipy.spatial.distance as dt
import scipy.sparse.linalg as la
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D
from scipy import ndimage
import time


def Shape(geometry,diel,d,colormax=1,shapebarrier=0.5,plotDpi=100,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the shape of object as dot matrix.
    #Input:
    # --SC                                                         SolutionClass
    #   Solution Class.
    # --FullLattice   Boolean
    #   If true. Geometry is a full n*n*n matrix. diel=0 for point seen as air.
    """
    #d = SC.Get_d()
    N=round(np.shape(geometry)[0]/3)
    if(N!=np.shape(diel)[0]/3):
        print("size not equal!")
    geometry=np.reshape(geometry,(N,3))
    #geometry = SC.Get_geometry()
    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    [X,Y,Z] = map(int,list(np.amax(geometry,axis=0)+1))
    Axis_max = max(X,Y,Z)*1.2*d

    diel=np.reshape(diel,(N,3))
    diel=diel[:,0]
    #diel = SC.Get_diel()[:,0]


    cmaparg = 'Spectral_r'
    minn, maxx = 0, colormax
    norm = matplotlib.colors.Normalize(minn, maxx)
    colorset = cm.ScalarMappable(norm=norm, cmap=cmaparg)
    colorset.set_array([])
    if FullLattice:
        index = np.where(diel>shapebarrier)
        diel = diel[index]
        geometry = geometry[index]
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    geo_dic = set()
    surf_list = []
    surf_color = []
    x_grid, y_grid, z_grid = (np.indices((X+1,Y+1,Z+1))-0.5)*d
    filled, colors = np.zeros((X,Y,Z)), np.zeros((X,Y,Z))
    for i,pos in enumerate(geometry):
        geo_dic.add((pos[0],pos[1],pos[2]))
    for i,pos in enumerate(geometry):
        if (pos[0]+1,pos[1],pos[2]) not in geo_dic or (pos[0]-1,pos[1],pos[2]) not in geo_dic or\
           (pos[0],pos[1]+1,pos[2]) not in geo_dic or (pos[0],pos[1]-1,pos[2]) not in geo_dic or\
           (pos[0],pos[1],pos[2]+1) not in geo_dic or (pos[0],pos[1],pos[2]-1) not in geo_dic:
            filled[pos[0],pos[1],pos[2]] = 1
            colors[pos[0],pos[1],pos[2]] = diel[i]
    surf_list = np.array(surf_list)
    surf_color = np.array(surf_color)
    # print(x_grid.shape,y_grid.shape,z_grid.shape,filled.shape,colors.shape)
    colors = cm.Spectral_r(norm(colors))
    ln=ax.voxels(x_grid,y_grid,z_grid,filled.astype(bool),facecolors=colors,edgecolor ='white',linewidth=0.2)
    #ax2.scatter(geometry[:,0]*d,geometry[:,1]*d,geometry[:,2]*d,c=E_tot_abs, s=15, cmap=cmaparg)
    ax.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax.set_xlabel("x[nm]")
    ax.set_ylabel("y[nm]")
    ax.set_zlabel("z[nm]")
    ax.grid(False)
    fig.colorbar(colorset, shrink=0.9, aspect=10)
     
    #plt.show()
    #plt.savefig("./optimization_geometries/Iteration{}.png".format(it))
    if iteration==-1:
        plt.savefig("Space.png")
    else:
        fig.suptitle("iteration{}".format(iteration))
        plt.savefig(position+"{}Space.png".format(str(iteration).zfill(decimal)),dpi=plotDpi) 
    #plt.show()
def EField(geometry,diel,d,k_dir,E_dir,E_tot,colormax=1,iteration=-1,position="./",decimal=0,FullLattice=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_abs = np.abs(np.sqrt(np.sum((E_tot)**2,axis=1)))
    
    """
    if FullLattice:
        index = np.where(diel!=0)                                      #disabled all codes with "diel"
        geometry = geometry[index]
        E_tot = E_tot[index]
        E_tot_abs = E_tot_abs[index]
    """
    cmaparg = 'Spectral_r'
    minn, maxx = E_tot_abs.min(), E_tot_abs.max()
    print(minn,maxx,np.argmax(E_tot_abs))
    norm = matplotlib.colors.Normalize(minn, maxx)
    colorset = cm.ScalarMappable(norm=norm, cmap=cmaparg)
    colorset.set_array([])
    
    
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.quiver(geometry[:,0]*d,geometry[:,1]*d,geometry[:,2]*d,np.real(E_tot[:,0]),np.real(E_tot[:,1]),np.real(E_tot[:,2]),
                length=10, lw=1)
    ax1.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax1.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax1.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax1.set_xlabel("x[nm]")
    ax1.set_ylabel("y[nm]")
    ax1.set_zlabel("z[nm]")
    ax1.grid(False)
    fig1.suptitle("E field - Arrow plot\n {}".format(E_dir))
    plt.savefig(position+"{}E_field_arrow.png".format(str(iteration).zfill(decimal)))

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    geo_dic = set()
    surf_list = []
    surf_color = []
    x_grid, y_grid, z_grid = (np.indices((X+1,Y+1,Z+1))-0.5)*d
    filled, colors = np.zeros((X,Y,Z)), np.zeros((X,Y,Z))
    for i,pos in enumerate(geometry):
        geo_dic.add((pos[0],pos[1],pos[2]))
    for i,pos in enumerate(geometry):
        if (pos[0]+1,pos[1],pos[2]) not in geo_dic or (pos[0]-1,pos[1],pos[2]) not in geo_dic or\
           (pos[0],pos[1]+1,pos[2]) not in geo_dic or (pos[0],pos[1]-1,pos[2]) not in geo_dic or\
           (pos[0],pos[1],pos[2]+1) not in geo_dic or (pos[0],pos[1],pos[2]-1) not in geo_dic:
            filled[pos[0],pos[1],pos[2]] = 1
            colors[pos[0],pos[1],pos[2]] = E_tot_abs[i]
    surf_list = np.array(surf_list)
    surf_color = np.array(surf_color)
   

    colors = cm.Spectral_r(norm(colors))
    ax2.voxels(x_grid,y_grid,z_grid,filled.astype(bool),facecolors=colors,linewidth=0.5)
    

    ax2.set_xlim(-(Axis_max-X*d)/2,(Axis_max+X*d)/2)
    ax2.set_ylim(-(Axis_max-Y*d)/2,(Axis_max+Y*d)/2)
    ax2.set_zlim(-(Axis_max-Z*d)/2,(Axis_max+Z*d)/2)
    ax2.set_xlabel("x[nm]")
    ax2.set_ylabel("y[nm]")
    ax2.set_zlabel("z[nm]")
    ax2.grid(False)
    fig2.colorbar(colorset, shrink=0.9, aspect=10)
    if iteration==-1:
        fig2.suptitle("E field - Scatter plot\n, {}".format(k_dir)) 
        plt.savefig("E_field.png")
    else:
        fig2.suptitle("E field - Scatter plot\n, {} - iteration{}".format(k_dir,iteration)) 
        plt.savefig(position+"{}E_field.png".format(str(iteration).zfill(decimal)))
    plt.show()


    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """

def EField_slice(geometry,diel,d,k_dir,E_dir,E_tot,plotDpi=100, Elimit=False, Elimitlow=-1,Elimithigh=-1,Xslice=-1,Yslice=-1,Zslice=-1,iteration=-1,position="./",decimal=0,FullLattice=False,plotlimit=False):
    """Plot the E field of object as arrow matrix.
    # Input:
    # --SC         SolutionClass
    #   Solved Solution Class.
    # --idx1,idx2  int
    #   Indexs of the target instance.
    """
    #print(geometry)

    N=round(np.shape(geometry)[0]/3)
    print(N)
    geometry=np.reshape(geometry,(N,3))
    
    print(geometry.shape)
    #diel=np.reshape(diel,(N,3))
    #diel=diel[:,0]
    

    for i in range(3):
        geometry[:,i] -= np.amin(geometry[:,i])
    
    #print(geometry)
    [X,Y,Z] = list(np.amax(geometry,axis=0)+1)
    Axis_max = max(X,Y,Z)*1.2*d
    
    print("{}, {}, {}".format(X,Y,Z))

    E_tot = E_tot.reshape(int(E_tot.size/3),3)
    E_tot_real=E_tot.real
    E_tot_imag=E_tot.imag
    E_tot_abs = np.sqrt(np.sum(E_tot_real**2+E_tot_imag**2,axis=1))
    
    slicedim=-1
    

    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Y, Z))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((Z, X))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((X, Y))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-2]][geometry[i][slicedim-1]]=ele

    """
    if Xslice!=-1:
        slicedim=0
        Eslice=np.zeros((Z, Y))
    if Yslice!=-1:
        slicedim=1
        Eslice=np.zeros((X, Z))
    if Zslice!=-1:
        slicedim=2
        Eslice=np.zeros((Y, X))
    
    slicepos= max([Xslice, Yslice, Zslice])
    for i, ele in enumerate(E_tot_abs):
        if geometry[i][slicedim]==slicepos:
            Eslice[geometry[i][slicedim-1]][geometry[i][slicedim-2]]=ele
    """
    rotated_img = ndimage.rotate(Eslice, 90)
    
    fig1 = plt.figure(figsize=(10, 10))
    plt.imshow(rotated_img, cmap='jet', interpolation='bilinear')
    if Elimit:
        plt.clim(Elimitlow, Elimithigh)
    plt.colorbar()
    plt.savefig(position+"Model{} E_slice_{}at{}.png".format(iteration, (["X", "Y", "Z"])[slicedim], slicepos), dpi=plotDpi) 
