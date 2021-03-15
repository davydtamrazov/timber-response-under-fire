#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from skimage.measure import label, regionprops

class Section():
    '''
    Timber section under temperature loading discretized with fibers.
    
    Attributes:
        w: An integer section width (mm).
        h: An integet section height (mm).
        Ti: An integer ambient temperature (Celsius)
        Tp: An integer surface temperature (Celsius)
        texp: A float fire exposure time (min)
        z: A 2D array of fiber split locations along the vertical axis (mm)
        y: A 2D array of fiber split locations along the horizontal axis (mm)
        temp_profile: A 2D array of temperatures at each cross split (Celcius)
        temp_dict: A dictionary of indexed unique temperatures in fibers
        patch_list: A list of patches in OpenSeesPy format
        
    Example usage:
        s = Section(100,100)
        s.discretize([5,5])
        s.set_exposure_time(5)
        s.apply_temperature()
        s.plot()
    '''

    NOMINAL_CHAR_RATE = 0.635 # (mm/min)
    
    def __init__(self, width, height, Ti=20, Tp=300):
        '''Initialize timber section parameters'''
        self.w = width
        self.h = height
        self.Ti = Ti
        self.Tp = Tp
        self.texp = 0
        self.z = np.array([[0, height], [0, height]])
        self.y = np.array([[width, width], [0, 0]])
        self.temp_profile = np.zeros((2,2)) + Ti
        self.temp_dict = {self.Ti : 1}
        self.patch_list = []
        
    def __str__(self):
        return (f'{self.w}mm by {self.h}mm timber section.')
        
    def discretize(self, mesh_size):
        '''
        Generate mesh within section domain.
        
        Args:
            mesh_size: A tuple mesh size in along z and y axes
        '''
        fibers = np.meshgrid(
            np.linspace(0, self.w, int(self.w/mesh_size[0])+1),
            np.linspace(0, self.h, int(self.h/mesh_size[1])+1)
            )
        
        self.z = fibers[0]
        self.y = np.flip(fibers[1],0)
        
        self.temp_profile = np.zeros((self.z.shape)) + self.Ti
        
    def set_exposure_time(self, exposure_time):
        '''
        Update exposure time.
        
        Args:
            exposure_time: A float exposure time (min).
        '''
        self.texp = exposure_time
        
    def apply_temperature(self, temp_pdepth=35, side=[1,1,1,1]):
        '''
        Calculate temperature profile in discretised section.
        
        Args:
            temp_pdepth: A float temperature penetration depth (default=35mm)
            side: A list of flags side exposures [left, down, right, up],
                  (default=4-sided exposure)
        '''
        
        cdepth = 2.58*self.NOMINAL_CHAR_RATE*self.texp**(0.813)
        
        # Find maximum of temperature profiles due to fire applied to each side
        self.temp_profile = np.dstack([
            self._temp_bchar(self.z, self.Ti, self.Tp, 
                             temp_pdepth, cdepth, side[0]),
            self._temp_bchar(self.y, self.Ti, self.Tp, 
                              temp_pdepth, cdepth, side[1]),
            self._temp_bchar(self.w-self.z, self.Ti, self.Tp, 
                              temp_pdepth, cdepth, side[2]),
            self._temp_bchar(self.h-self.y, self.Ti, self.Tp, 
                              temp_pdepth, cdepth, side[3])
            ]).max(axis=2)
        
        # Update patch definition
        self._split_fibers()


    def plot(self, lw=1, factor=10, disp_fiber=True, levels=10):
        '''
        Generate contour plot of the temperature profile across the section.
        
        Args:
            lw: Plot linewidth (default=1).
            factor: Plot scale factor wrt section dimensions (default=10)
            disp_fiber: A flag to display fibers on the plot (default=True)
            levels: An integer of contour plot levels (default=10)
        '''
        fig, ax = plt.subplots(figsize=(self.w/factor, self.h/factor), dpi=200)
        ax.set_axis_off()
        ax.set(xlim=[0, self.w], ylim=[0, self.h])
        
        # Display grid of fibers
        if disp_fiber:
            ax.plot(self.z, self.y, c='k', ls='-', lw = lw)
            ax.plot(self.z.T, self.y.T, c='k', ls='-', lw = lw)

        # Display temperature profile
        ax.contourf(self.z, self.y, self.temp_profile, cmap='YlOrRd', 
                      vmin = self.Ti, vmax = self.Tp, levels=levels, zorder=2)
        
        # Display charred area
        ax.add_patch(Rectangle((0,0), self.w, self.h, fill=True, 
                                fc = 'black', alpha = 0.75, hatch='/')) 
        
    
    def _temp_bchar(self, d, Ti, Tp, temp_pdepth, cdepth, exposed=True):
        '''
        Calculate temperature at an effective depth below the char layer
        using Eurocode 5 temperature profile equation.
        
        Args:
            d: A float depth from the section edge (mm)
            Ti: An integer ambient temperature (Celsius)
            Tp: An integer surface temperature (Celsius)
            temp_pdepth: A float temperature penetration depth (default=35mm)
            cdepth: Charred depth from section edge (mm)
        '''
        
        d_eff = (d - exposed*cdepth).clip(max=temp_pdepth)
        d_eff = np.where(d_eff < 0, np.nan, d_eff)
        Tx = Ti + exposed * (Tp - Ti) * (1-d_eff/temp_pdepth)**2
        return Tx
    

    def _split_fibers(self):
        '''
        Split section into fiber patches based on temperature in each fiber.
        '''
        
        # Calculate temperature for each fiber
        fiber_temp = self.temp_profile.copy()
        fiber_temp += np.roll(fiber_temp, shift=1, axis=0)
        fiber_temp += np.roll(fiber_temp, shift=1, axis=1)
        fiber_temp /= 4
        fiber_temp = fiber_temp[1:,1:]
        fiber_temp = fiber_temp[~np.isnan(fiber_temp).all(axis=1), :]
        fiber_temp = fiber_temp[:, ~np.isnan(fiber_temp).all(axis=0)]
        fiber_temp = np.where(np.isnan(fiber_temp), self.Tp, fiber_temp)
        
        if fiber_temp.shape[0]==0 or fiber_temp.shape[1]==0: 
            raise Exception('No remaining section area at the specified exposure time.')
    
        self.temp_dict = {k:(v+1) for v, k in 
                              enumerate(np.unique(fiber_temp))}
        
        # Split section into fiber patches defined in OpenSeesPy format
        self.patch_list = []
        for t in self.temp_dict.keys():
            temp_mask = (fiber_temp==t).astype(int)
            patches = label(temp_mask, connectivity=1)
            for s in regionprops(patches):
                z_max = self.z[:,1:][s.slice].max()
                z_min = self.z[:,:-1][s.slice].min()
                y_max = self.y[:-1,:][s.slice].max()
                y_min = self.y[1:,:][s.slice].min()
                c1, c2 = (y_min, z_min), (y_max, z_max)
                nz, ny = self.z[s.slice].shape
                self.patch_list.append([self.temp_dict[t], ny, nz, c1, c2])
                
                
# t = Section(200,100)
# t.discretize([25,10])
# t.set_exposure_time(20)
# t.apply_temperature(side=[1,0,0,1])
# t.plot()

# print(np.array(t.patch_list))