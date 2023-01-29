import numpy as np
from path import Path
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
import configparser
import plot_func
config1 = configparser.ConfigParser()
config1.read('task.ini')
print(config1['Model Name']['name'])
modelToPath={'DDA verify':'DDA_verify_path', 'DDA input':'DDA_input_path', 'EvoOpt 2D input':'EvoOpt_2D_input_path', 'EvoOpt 2D input periodic':'EvoOpt_2D_input_periodic_path', 'NN data generate':'NN_data_generate_path', 'NN data generate 2D':'NN_data_generate_2D_path' }
config=configparser.ConfigParser()
config.read(config1["Path"][modelToPath[config1['Model Name']['name']]])

# print(config1["Path"][modelToPath[config1['Model Name']['name']]])

if config1['Model Name']['name']=="DDA input" and config["Plot"]["plot"]=="True":

    objective_number = 1
    pos=config["Output"]["saveDir"]
    print(pos)
    itStart = config["Plot"]["itStart"]
    itEnd = config["Plot"]["itEnd"]
    
    datacommon=np.genfromtxt(pos+"commondata.txt")
    N=int(np.real(datacommon[3]))
    print(N)
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    xPlot=config["Plot"]["xPlot"]
    xslice=int(config["Plot"]["x"])
    yPlot=config["Plot"]["yPlot"]
    yslice=int(config["Plot"]["y"])
    zPlot=config["Plot"]["zPlot"]
    zslice=int(config["Plot"]["z"])
    shapeSolidPlot=config["Plot"]["shapeSolidPlot"]
    shapePlot=config["Plot"]["shapePlot"]
    colorMax=float(config["Plot"]["colorMax"])
    shapeBarrier=float(config["Plot"]["shapeBarrier"])
    plotDpi=int(config["Plot"]["plotDpi"])
    ELimit=False
    azim=float(config["Plot"]["azim"])
    elev=float(config["Plot"]["elev"])
    if config["Plot"]["ELimit"]=="True":
        ELimit=True
    ELimitLow=config["Plot"]["ELimitLow"]
    ELimitHigh=config["Plot"]["ELimitHigh"]

    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure\\","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output\\","Model_results0"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            
            if(nameit >= int(itStart) and nameit <= int(itEnd)):
                if shapeSolidPlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"ShapeSolid\\", decimal=dec, FullLattice=True)
                if shapePlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"Shape\\", decimal=dec, FullLattice=False)
                if xPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Xslice=xslice,position=pos+"E-field\\")
                if yPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Yslice=yslice,position=pos+"E-field\\")
                if zPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Zslice=zslice,position=pos+"E-field\\")

if config1['Model Name']['name']=="EvoOpt 2D input" and config["Plot"]["plot"]=="True":

    objective_number = 1
    pos=config["Output"]["saveDir"]
    print(pos)
    itStart = config["Plot"]["itStart"]
    itEnd = config["Plot"]["itEnd"]
    
    datacommon=np.genfromtxt(pos+"commondata.txt")
    N=int(np.real(datacommon[3]))
    print(N)
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    xPlot=config["Plot"]["xPlot"]
    xslice=int(config["Plot"]["x"])
    yPlot=config["Plot"]["yPlot"]
    yslice=int(config["Plot"]["y"])
    zPlot=config["Plot"]["zPlot"]
    zslice=int(config["Plot"]["z"])
    shapeSolidPlot=config["Plot"]["shapeSolidPlot"]
    shapePlot=config["Plot"]["shapePlot"]
    colorMax=float(config["Plot"]["colorMax"])
    shapeBarrier=float(config["Plot"]["shapeBarrier"])
    plotDpi=int(config["Plot"]["plotDpi"])
    azim=float(config["Plot"]["azim"])
    elev=float(config["Plot"]["elev"])
    ELimit=False
    if config["Plot"]["ELimit"]=="True":
        ELimit=True
    ELimitLow=config["Plot"]["ELimitLow"]
    ELimitHigh=config["Plot"]["ELimitHigh"]

    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure\\","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output\\","Model_results0"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            
            if(nameit >= int(itStart) and nameit <= int(itEnd)):
                if shapeSolidPlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"ShapeSolid\\", decimal=dec, FullLattice=True)
                if shapePlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"Shape\\", decimal=dec, FullLattice=False)
                if xPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Xslice=xslice,position=pos+"E-field\\")
                if yPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Yslice=yslice,position=pos+"E-field\\")
                if zPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Zslice=zslice,position=pos+"E-field\\")

if config1['Model Name']['name']=="EvoOpt 2D input periodic" and config["Plot"]["plot"]=="True":

    objective_number = 1
    pos=config["Output"]["saveDir"]
    print(pos)
    itStart = config["Plot"]["itStart"]
    itEnd = config["Plot"]["itEnd"]
    
    datacommon=np.genfromtxt(pos+"commondata.txt")
    N=int(np.real(datacommon[3]))
    print(N)
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    xPlot=config["Plot"]["xPlot"]
    xslice=int(config["Plot"]["x"])
    yPlot=config["Plot"]["yPlot"]
    yslice=int(config["Plot"]["y"])
    zPlot=config["Plot"]["zPlot"]
    zslice=int(config["Plot"]["z"])
    shapeSolidPlot=config["Plot"]["shapeSolidPlot"]
    shapePlot=config["Plot"]["shapePlot"]
    colorMax=float(config["Plot"]["colorMax"])
    shapeBarrier=float(config["Plot"]["shapeBarrier"])
    plotDpi=int(config["Plot"]["plotDpi"])
    azim=float(config["Plot"]["azim"])
    elev=float(config["Plot"]["elev"])
    ELimit=False
    if config["Plot"]["ELimit"]=="True":
        ELimit=True
    ELimitLow=config["Plot"]["ELimitLow"]
    ELimitHigh=config["Plot"]["ELimitHigh"]

    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure\\","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output\\","Model_results0"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            
            if(nameit >= int(itStart) and nameit <= int(itEnd)):
                if shapeSolidPlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"ShapeSolid\\", decimal=dec, FullLattice=True)
                if shapePlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"Shape\\", decimal=dec, FullLattice=False)
                if xPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Xslice=xslice,position=pos+"E-field\\")
                if yPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Yslice=yslice,position=pos+"E-field\\")
                if zPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Zslice=zslice,position=pos+"E-field\\")

if config1['Model Name']['name']=="NN data generate" and config["Plot"]["plot"]=="True":

    objective_number = 1
    pos=config["Output"]["saveDir"]
    print(pos)
    itStart = config["Plot"]["itStart"]
    itEnd = config["Plot"]["itEnd"]
    
    datacommon=np.genfromtxt(pos+"commondata.txt")
    N=int(np.real(datacommon[3]))
    print(N)
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    xPlot=config["Plot"]["xPlot"]
    xslice=int(config["Plot"]["x"])
    yPlot=config["Plot"]["yPlot"]
    yslice=int(config["Plot"]["y"])
    zPlot=config["Plot"]["zPlot"]
    zslice=int(config["Plot"]["z"])
    shapeSolidPlot=config["Plot"]["shapeSolidPlot"]
    shapePlot=config["Plot"]["shapePlot"]
    colorMax=float(config["Plot"]["colorMax"])
    shapeBarrier=float(config["Plot"]["shapeBarrier"])
    plotDpi=int(config["Plot"]["plotDpi"])
    azim=float(config["Plot"]["azim"])
    elev=float(config["Plot"]["elev"])
    ELimit=False
    if config["Plot"]["ELimit"]=="True":
        ELimit=True
    ELimitLow=config["Plot"]["ELimitLow"]
    ELimitHigh=config["Plot"]["ELimitHigh"]

    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure\\","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output\\","Model_results"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            
            if(nameit >= int(itStart) and nameit <= int(itEnd)):
                if shapeSolidPlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"ShapeSolid\\", decimal=dec, FullLattice=True)
                if shapePlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"Shape\\", decimal=dec, FullLattice=False)
                if xPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Xslice=xslice,position=pos+"E-field\\")
                if yPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Yslice=yslice,position=pos+"E-field\\")
                if zPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Zslice=zslice,position=pos+"E-field\\")

if config1['Model Name']['name']=="NN data generate 2D" and config["Plot"]["plot"]=="True":

    objective_number = 1
    pos=config["Output"]["saveDir"]
    print(pos)
    itStart = config["Plot"]["itStart"]
    itEnd = config["Plot"]["itEnd"]
    
    datacommon=np.genfromtxt(pos+"commondata.txt")
    N=int(np.real(datacommon[3]))
    print(N)
    geometry=np.real(datacommon[4:(3*N+4)]).astype(int)
    d=np.real(datacommon[3*N+4])
    E_dir=np.real(datacommon[(3*N+5):(3*N+8)])
    k_dir=np.real(datacommon[(3*N+8):(3*N+11)])
    
    dec = 5
    cutnumber=13
    xPlot=config["Plot"]["xPlot"]
    xslice=int(config["Plot"]["x"])
    yPlot=config["Plot"]["yPlot"]
    yslice=int(config["Plot"]["y"])
    zPlot=config["Plot"]["zPlot"]
    zslice=int(config["Plot"]["z"])
    shapeSolidPlot=config["Plot"]["shapeSolidPlot"]
    shapePlot=config["Plot"]["shapePlot"]
    colorMax=float(config["Plot"]["colorMax"])
    shapeBarrier=float(config["Plot"]["shapeBarrier"])
    plotDpi=int(config["Plot"]["plotDpi"])
    azim=float(config["Plot"]["azim"])
    elev=float(config["Plot"]["elev"])
    ELimit=False
    if config["Plot"]["ELimit"]=="True":
        ELimit=True
    ELimitLow=config["Plot"]["ELimitLow"]
    ELimitHigh=config["Plot"]["ELimitHigh"]

    for filename in sorted(os.listdir(pos+"CoreStructure"), key = lambda x: int(x[cutnumber:x.index(".txt")])):
        if filename.endswith(".txt"):
            nameit=int((filename[cutnumber:])[:-4])
            print(nameit)
            CoreStructure=np.genfromtxt(os.path.join(pos+"CoreStructure\\","CoreStructure"+str(nameit)+".txt"),dtype=complex)
            Modelresults=np.genfromtxt(os.path.join(pos+"Model_output\\","Model_results"+"it"+str(nameit)+".txt"),dtype=complex)
        
            diel=np.real(CoreStructure[(0):(3*N)])
            E_tot=(Modelresults[(0):(3*N)])
            P_tot=(Modelresults[(3*N):(6*N)])
            
            if(nameit >= int(itStart) and nameit <= int(itEnd)):
                if shapeSolidPlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"ShapeSolid\\", decimal=dec, FullLattice=True)
                if shapePlot=="True":
                    plot_func.Shape(geometry, diel, d, azim, elev, colormax=colorMax,shapebarrier=shapeBarrier,plotDpi=plotDpi, iteration=nameit, position=pos+"Shape\\", decimal=dec, FullLattice=False)
                if xPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Xslice=xslice,position=pos+"E-field\\")
                if yPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Yslice=yslice,position=pos+"E-field\\")
                if zPlot=="True":
                    plot_func.EField_slice(geometry, diel, d, k_dir, E_dir, E_tot, Elimit=ELimit, Elimitlow=ELimitLow, Elimithigh=ELimitHigh, plotDpi=plotDpi, iteration=nameit, Zslice=zslice,position=pos+"E-field\\")

  
    
    
    
        
    
                


