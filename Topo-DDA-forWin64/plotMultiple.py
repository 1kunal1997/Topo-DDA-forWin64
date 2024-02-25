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
modelToPath={'DDA verify':'DDA_verify_path', 'DDA input':'DDA_input_path', 'EvoOpt 2D input':'EvoOpt_2D_input_path', 'EvoOpt 2D input periodic':'EvoOpt_2D_input_periodic_path', 'NN data generate':'NN_path' }
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

    weightArray = [0.2]
    stepArray = [0.1]
    penaltytypeArray = ["piecewise"]
    for i in range(1):
    #    weightArray[i] = round(weightArray[i],1)
        stepArray[i] = round(stepArray[i],1)
    #    print(weightArray[i])
        print(stepArray[i])
    
    #stepArray[0] = round(stepArray[0],1)
    #print(stepArray[0])
    
    for step in stepArray:
        for weight in weightArray:
            for penaltytype in penaltytypeArray:

                #pos = ".\\cylinder_it300_lam542_sym_epsilon_0.1_penaltytype_" + penaltytype + "_absolute_0.0to0.5\\"
                pos = "..\\Calculations\\SmallerGridCheck12x12x4\\"

                plt.figure(1)
                objectfunc = np.loadtxt(pos + "\\convergence.txt")
                #objectfuncwpenalty = np.loadtxt('cylinder_it200_sym_epsilon_0.100000_weight_0.300000/convergenceWithPenalty.txt')
                plt.plot(objectfunc, label = 'E')
                #plt.plot(objectfuncwpenalty, label = 'E - p')
                plt.legend(loc='lower right')
                plt.title('Object Function Plot')
                plt.ylabel('Object Function')
                plt.xlabel('Iteration #')
                #plt.xlim([0,200])
                #plt.yscale('symlog')
                #plt.ylim([1, 10**(6)])
                plt.rc('axes', titlesize=14)     # fontsize of the axes title
                plt.rc('axes', labelsize=12)    # fontsize of the x and y labels
                plt.savefig(pos + "\\Object Function.png")
                #plt.show()
                #plt.close(1)
                print("--------------Step is: " + str(step) + ", weight is, " + str(weight) + "-----------------")

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

    
    
    
        
    
                


