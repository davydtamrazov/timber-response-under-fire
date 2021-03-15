#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import FiberSection as fb
import OpenSeesAnalysis as ops
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------------
# NOTE:
#
# The following imports developed OpenSeesAnalysis and FiberSection files and
# displays the following outputs:
#   1. Section temperature profile under one-sided fire exposure
#   2. Section temperature profile under four-sided fire exposure
#   3. Applied total load vs lateral displacement plot (pin-pin column)
#   4. Critical load vs time of exposure to fire plot (pin-pin column)
#
#   * Please allow several minutes for the code to run.
#   ** Convergence issue displayed means that limit point has been reached for
#      that particular analysis and the warning message can be ignored.
# ----------------------------------------------------------------------------
#

#%%# Input Parameters
 
# Mesh size
mSize = 5
 
# Section dimensions (mm)
b = 400
h = 400
   
# Temperature distribution at distance x from the char front
Ti = 20     # Celcius
Tp = 300    # Celcius
a = 35      # mm

# Column height and eccentricity
L = 4000 #mm
e = L/1000

# Applied Load
P = 15000e3
numIncr = 1000


    
        
# m = StructuralModel(texp=1, w=100, h=200, mesh_size=[5,5])
# m.define_material()
# m.define_section()

# nodes = [[1, 0, 0],
#           [2, 0, 5000],
#           [3, 5000, 5000],
#           [4, 10000, 5000],
#           [5, 10000, 0]]

# elements = [[0, 1, 2],
#             [1, 2, 3],
#             [2, 3, 4],
#             [3, 4, 5]]

# fixities = [[1, 1, 1, 0],
#             [5, 0, 1, 0]]

# # nodes = [[1, 0, 0],
# #          [2, 10, 5000],
# #          [3, 0, 10000]]

# # elements = [[0, 1, 2],
# #             [1, 2, 3],
# #             [2, 3, 4]]

# # fixities = [[1, 1, 1, 1]]

# m.define_geometry(nodes, elements, fixities)

# loads = [[3, 0, -200, 0]]
# m.define_loads(loads)

# recorders = [['disp_node3.txt', 3, 2, 'disp']]
# m.define_recorders(recorders)

# # m.material.plot(5)
# m.section.plot()
# # m.plot_structure()

# m.analyse(1000)
# m.plot_deformed_shape()
# m.wipe_model()


#%% Example 1: Plot fiber section at time t
plot = True
disp_fiber = True

# Time of exposure (minutes)
t = 120

# One-sided exposure
[FiberList, MatList] = fb.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, 1)

# Four-sided exposure
[FiberList, MatList] = fb.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, 4)



#%% Example 2: Structural response of a pin-pin column under four-side fire exposure
n = 25
tList = np.linspace(0,120,n)
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.coolwarm(np.linspace(0,1,n)))

Pmax = []

# Calculate critical load at ambient temperature
force,_ = ops.RunAnalysis(0,b,h,L,Ti,Ti,a,mSize,e,P,numIncr)
P0 = max(force)

# Set up plot
fig, ax = plt.subplots(dpi=175)
fig.set_size_inches(7,7)
ax.set( xlabel = 'Lateral displacement of middle node (mm)', ylabel = 'Applied total load (kN)')
ax.grid('on')
ax.minorticks_on()
ax.grid(b = True, which = 'major', linestyle='-', linewidth=0.5, alpha=0.75)
ax.grid(b = True, which = 'minor', color = 'gray', linestyle='--', linewidth=0.25, alpha=0.5)

# Loop through time steps
for t in tList:
    [force, disp] = ops.RunAnalysis(t,b,h,L,Ti,Tp,a,mSize,e,P,numIncr)
    Pmax.append(max(force)/1000)
    ax.plot(disp,force/1000,linewidth = 2, label = 't = ' + str(t) + ' min')

fig.legend(framealpha=1)
fig.tight_layout()

#%% Plot Critical Load Reduction Factor
fig, ax = plt.subplots(dpi=175)
fig.set_size_inches(7,7)
ax.set(xlabel = 'Time of exposure to fire / mins', ylabel = 'Critical Load Reduction Factor',
        ylim = [0,0.85], yticks = np.arange(0,1.05,.05), xticks = np.arange(0,200,20))
ax.grid('on')
ax.minorticks_on()
ax.grid(b = True, which = 'major', linestyle='-', linewidth=0.5, alpha=0.75)
ax.grid(b = True, which = 'minor', color = 'gray', linestyle='--', linewidth=0.25, alpha=0.5)
ax.plot(tList,[i/P0*1000 for i in Pmax], color = 'black', linewidth = 2)
fig.tight_layout()

