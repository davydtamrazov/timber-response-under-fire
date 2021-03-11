#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Function to compute temperature profile across the section for a given 
# exposure time (t) and discretize the section into fibers
#
# ----------------------------------------------------------------------------
# INPUT:
#   t           = time of exposure
#   b           = section width (mm)
#   h           = section height (mm)
#   Ti          = ambient temperature (Celsius)
#   Tp          = surface temperature (Celsius)
#   a           = temperature penetration depth (mm)
#   mSize       = mesh size (mm)
#   plot        = boolean to display section temperature profile plot
#   disp_fiber  = boolean to display mesh on the plot
#   side        = number of sides from which fire is applied (either 1 or 4)
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# OUTPUT:
#   patchList   = list of patches in OpenSeesPy format: 
#                 (mattag, ny, nz, coord1, coord2)
#   mattag      = dictionary of unique temperatures in discretized section and
#                 corresponding material tags
# ----------------------------------------------------------------------------
#

def GetFibers(t,b,h,Ti,Tp,a,mSize,plot,disp_fiber,side = 4):
    
    ## Calculate temperature profile across the section
    # Charring depth at time (t) in mm   
    bn = 0.635 #mm/min
    c = 2.58*bn*t**(0.813)

    # Generate mesh within section domain
    x, y = np.meshgrid(np.linspace(0,b,int(b/mSize)+1), np.linspace(0,h,int(h/mSize)+1))
    
    # Calculate temperature profile in x and y directions (depending on side
    # parameter) and find overall maximum temperature profile section-wide
    if side == 4:
        Tx = np.maximum(Temp(x,a,c,Ti,Tp),Temp(b-x,a,c,Ti,Tp))
        Ty = np.maximum(Temp(y,a,c,Ti,Tp),Temp(h-y,a,c,Ti,Tp))
        z = np.maximum(Tx,Ty)
    elif side == 1:
        z = Temp(x,a,c,Ti,Tp)
        
    # Map temperature (z) values to meshgrid format
    mask = np.zeros_like(z,dtype=bool)
    z = np.ma.array(z,mask=mask)
    
    # Generate  contour plot of the temperature profile across the section
    if plot:
        fig, ax = plt.subplots(dpi=175)
        factor = round(max(b,h)/10,0)
        fig.set_size_inches(b/factor,h/factor)

        ax.set(xlim = [0,b], ylim = [0,h])
        
        if disp_fiber:
            plt.plot(x,y,marker='.', color='k', linestyle='-', markersize = 3, linewidth = 0.75)
            plt.plot(y,x,color='k', linestyle='-', linewidth = 0.75)
            plt.plot(y,h-x,color='k', linestyle='-', linewidth = 0.75)  
        else: plt.contour(x,y,z,colors='gray', linewidths=0.75,linestyles='--')      
        
        plt.contourf(x,y,z,cmap='YlOrRd',vmin = Ti, vmax = 300, zorder=2)
        
        ax.set_axis_off()
        ax.add_patch(Rectangle((0,0),b,h,fill=True, facecolor = 'black', alpha = 0.75, hatch='/'))     
    
    
    ## Calculate average temperature of each mesh element (fiber) and group 
    ## into patches based on those temperatures
    
    # Mesh size
    mDim = z.shape
    
    # Initialize an empty list of fiber temperatures
    fbT = []
    
    # Loop through the each point in the meshgrid and calculate the average 
    # temperature at four adjacent points that define that element
    for idR in range(mDim[0]-1):
        tmp = []
        for idC in range(mDim[1]-1):
            Tmean = np.mean(z[idR: idR + 2, idC: idC + 2])
            tmp.append(round(Tmean,12))        
        fbT.append(tmp)
    
    # Replace NaN values (char layer) with the burning temperature
    fbT = np.where(np.isnan(fbT), Tp, fbT)
    
    # Define dictionary of unique fiber temperatures and the corresponding
    # material tags 
    mattag = np.unique(fbT)  
    mattag = {k:(v+1) for v, k in enumerate(mattag)}

    # Calculate number of temperature layers given  mesh size
    Nlayers = int(0.5*min(b,h)/mSize)
    
    # Initialize an empty list of patch parameters
    patchList = []  # (mattag, ny, nz, coord1, coord2)
    
    # Loop through each layer and append patch parameters to the list
    for i in range(Nlayers): 
        
        # Check if we are in the core, which can be defined as one rectangular
        # patch. This is because beyond this point nomore layers are needed,
        # as temperature remains unchanged with depth
        if fbT[i][i] == Ti:
            c1,c2 = [(mSize*(i),mSize*(i)), (b-mSize*(i),h-mSize*(i))]
            ny,nz = [(b-2*mSize*(i))/mSize,(h-2*mSize*(i))/mSize]
            patchList.append([mattag[Ti],ny,nz,c1,c2])
            break
        
        # Otherwise define four corners and two patches in each horintal and
        # vertical directions
        else:
            # Define 4 corners
            ny,nz = [1,1]
            
            coord1,coord2 = [(mSize*(i),mSize*(i)), (mSize*(i+1),mSize*(i+1))]         
            patchList.append([mattag[fbT[i][i]],ny,nz,coord1,coord2])
            
            coord1,coord2 = [(b-mSize*(i+1),mSize*(i)),(b-mSize*(i),mSize*(i+1))]
            patchList.append([mattag[fbT[i][i]],ny,nz,coord1,coord2])
            
            coord1,coord2 = [(mSize*(i),h-mSize*(i+1)), (mSize*(i+1),h-mSize*(i))]         
            patchList.append([mattag[fbT[i][i]],ny,nz,coord1,coord2])
            
            coord1,coord2 = [(b-mSize*(i+1),h-mSize*(i+1)),(b-mSize*(i),h-mSize*(i))]
            patchList.append([mattag[fbT[i][i]],ny,nz,coord1,coord2])
            

            # Define 2 horizontal fibers
            ny,nz = [(b-2*mSize*(i+1))/mSize,1]
            
            c1,c2 = [(mSize*(i+1),mSize*(i)), (b-mSize*(i+1),mSize*(i+1))]
            patchList.append([mattag[fbT[i][i+1]],ny,nz,c1,c2])
            
            c1,c2 = [(mSize*(i+1),h-mSize*(i+1)), (b-mSize*(i+1),h-mSize*(i))]
            patchList.append([mattag[fbT[i][i+1]],ny,nz,c1,c2])
            
            # Define 2 vertical fibers
            ny,nz = [1,(h-2*mSize*(i+1))/mSize]
            
            c1,c2 = [(mSize*(i),mSize*(i+1)), (mSize*(i+1),h-mSize*(i+1))]         
            patchList.append([mattag[fbT[i+1][i]],ny,nz,c1,c2])
            
            c1,c2 = [(b-mSize*(i+1),mSize*(i+1)), (b-mSize*(i),h-mSize*(i+1))]         
            patchList.append([mattag[fbT[i+1][i]],ny,nz,c1,c2])

    return patchList, mattag


# Function to calculate temperature at depth (x) from the surface
#
# ----------------------------------------------------------------------------
# INPUT:
#   x           = depth from the surface
#   a           = temperature penetration depth (mm)
#   c           = char depth (mm)
#   Ti          = ambient temperature (Celsius)
#   Tp          = surface temperature (Celsius)
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# OUTPUT:
#   T           = temperature at depth (x) from the surface (NaN if char)
# ----------------------------------------------------------------------------
#
    
def Temp(x,a,c,Ti,Tp):
    xi = x-c
    xi[xi<0] = np.NaN
    T = Ti+(Tp-Ti)*((1-xi/a).clip(0))**2
    return T