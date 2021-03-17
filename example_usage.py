#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from opensees_analysis import StructuralModel
import matplotlib.pyplot as plt
import numpy as np
 
# Section dimensions (mm)
w = 200
h = 200
mesh_size = [5,5]

# Applied loads
Py = -1000
num_incr = 1000

nodes = [
    [1, 0, 0],
    [2, 0, 5000],
    [3, 5000, 5000],
    [4, 10000, 5000],
    [5, 10000, 0]
    ]

elements = [
    [0, 1, 2, 1],
    [1, 2, 3, 1],
    [2, 3, 4, 1],
    [3, 4, 5, 1]]

fixities = [[1, 1, 1, 0],
            [5, 0, 1, 0]]

loads = [[3, 0, Py, 0]]


def plot_displacement(time, disp, save, xlim, ylim, marker='o', c='k', size=10):
    fig, ax = plt.subplots(dpi=75)
    ax.grid(True, which='both', alpha=0.5)
    ax.set(xlabel='Exposure time (min)', ylabel='Vertical displacement (mm)',
            xlim=xlim, ylim=ylim)
    ax.plot(time, disp, c='k', zorder=0)
    ax.scatter(time[-1], displacement[-1], marker=marker, s=size, c=c, zorder=1)
    fig.savefig(save, transparent=True)
    

displacement = []

t_ig = 2 # 2 minutes to ignition
texp = np.arange(0, 100, 0.25)

T_range = np.linspace(25,300, num=texp[np.where(texp<=t_ig)].shape[0])

for i, t in enumerate(texp):
    print(t)
    if t <= t_ig:
        m = StructuralModel(0, w, h, T_range[i], mesh_size, nodes, elements, fixities, loads)
    else:
        m = StructuralModel(t, w, h, T_range[-1], mesh_size, nodes, elements, fixities, loads)
        
    m.define_material()
    m.define_fiber_section()
    m.define_geometry(num_integ=5)
    m.define_loads()
 
    ok, disp = m.analyse(num_incr, 3)
    
    
    if ok == 0:
        displacement.append(disp[1])
        m.section.plot(factor=50, levels=40, save=f'section/exposure{round(t*100)}')
        m.plot_deformed_shape(xlim=[-1000,12000], ylim=[-500,5000], scale=4, 
                              arrow_len=1, arrow_width=1.5, save=f'deformation/exposure{round(t*100)}')
        plot_displacement(texp[:i+1], displacement, save=f'displacement/exposure{round(t*100)}',
                          xlim=[0,20], ylim=[-75,-15], marker='o', c='k')
        m.wipe_model()
    else:
        displacement.append(-70)
        m.section.plot(factor=50, levels=40, save=f'section/exposure{round(t*100)}')
        plot_displacement(texp[:i+1], displacement, save=f'displacement/exposure{round(t*100)}',
                          xlim=[0,20], ylim=[-75,-15], marker='x', c='orangered', size=20)
        print(f'Structure failed in {t} minutes')
        break
    
#%%

# #%% Example 1: Plot fiber section at time t
# plot = True
# disp_fiber = True

# # Time of exposure (minutes)
# t = 120

# # One-sided exposure
# [FiberList, MatList] = fb.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, 1)

# # Four-sided exposure
# [FiberList, MatList] = fb.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, 4)



# #%% Example 2: Structural response of a pin-pin column under four-side fire exposure
# n = 25
# tList = np.linspace(0,120,n)
# plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.coolwarm(np.linspace(0,1,n)))

# Pmax = []

# # Calculate critical load at ambient temperature
# force,_ = ops.RunAnalysis(0,b,h,L,Ti,Ti,a,mSize,e,P,numIncr)
# P0 = max(force)

# # Set up plot
# fig, ax = plt.subplots(dpi=175)
# fig.set_size_inches(7,7)
# ax.set( xlabel = 'Lateral displacement of middle node (mm)', ylabel = 'Applied total load (kN)')
# ax.grid('on')
# ax.minorticks_on()
# ax.grid(b = True, which = 'major', linestyle='-', linewidth=0.5, alpha=0.75)
# ax.grid(b = True, which = 'minor', color = 'gray', linestyle='--', linewidth=0.25, alpha=0.5)

# # Loop through time steps
# for t in tList:
#     [force, disp] = ops.RunAnalysis(t,b,h,L,Ti,Tp,a,mSize,e,P,numIncr)
#     Pmax.append(max(force)/1000)
#     ax.plot(disp,force/1000,linewidth = 2, label = 't = ' + str(t) + ' min')

# fig.legend(framealpha=1)
# fig.tight_layout()

# #%% Plot Critical Load Reduction Factor
# fig, ax = plt.subplots(dpi=175)
# fig.set_size_inches(7,7)
# ax.set(xlabel = 'Time of exposure to fire / mins', ylabel = 'Critical Load Reduction Factor',
#         ylim = [0,0.85], yticks = np.arange(0,1.05,.05), xticks = np.arange(0,200,20))
# ax.grid('on')
# ax.minorticks_on()
# ax.grid(b = True, which = 'major', linestyle='-', linewidth=0.5, alpha=0.75)
# ax.grid(b = True, which = 'minor', color = 'gray', linestyle='--', linewidth=0.25, alpha=0.5)
# ax.plot(tList,[i/P0*1000 for i in Pmax], color = 'black', linewidth = 2)
# fig.tight_layout()

