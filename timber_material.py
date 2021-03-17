#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class Material():
    '''
    Definition of the timber material properties as per Eurocode 5
    
    Attributes:
        f_cu: A float ultimate compressive strength, MPa.
        f_tu: A float ultimate tensile strength, MPa.
        f_r: A float compression stress at rupture, MPa.
        E_wc: A float elastic modulus of elasticity in compression, MPa.
        E_wt: A float elastic modulus of elasticity in tension, MPa.
        E_q: A float softening modulus of elasticity in compression, MPa.
    
    Example usage:
        mat = Material(f_cu=45, f_tu=80, eps_y=4e-3, 
                       eps_r=13e-3, eps_tu=7.1e-3)
        mat.plot_stress_strain(200, 20)
    '''
    
    def __init__(self, f_cu, f_tu, eps_y, eps_r, eps_tu):
        '''
        Initialise timber material properties.
        
        Args:
            f_cu: A float ultimate compressive strength, MPa (default=45MPa).
            f_tu: A float ultimate tensile strength, MPa (default=80MPa).
            eps_y: A float yield compression strain (default=4e-3).
            eps_r: A float compressive strain at rupture (default=13e-3).
            eps_tu: A float tensile strain at rupture (default=7.1d-3).
        '''
        
        self.f_cu = f_cu
        self.f_tu = -f_tu
        self.f_r = self.f_cu*0.75
        
        self.E_wc = self.f_cu/eps_y
        self.E_wt = -self.f_tu/eps_tu
        self.E_q = (self.f_r - self.f_cu)/(eps_r-eps_y)

    def get_stress_strain(self, T, Ti):
        '''
        Update strengths and moduli of elasticity with temperature reduction
        factors as per Eurocode 5 (EN 1995-1-2) and calculate multilinear
        stress-strain points.
        
        Args:
            T: An integer temperature at the surface, Celsius.
            Ti: An integer ambient temperature, Celsius.
            
        Returns:
            stress: An array of stress points along stress-strain curve, MPa.
            strain: An array of strain points along stress-strain curve.
        '''

        # Strength and modulus of elasticity temperature reduction factors 
        # (as per Eurocode 5: Part 1-2)
        
        if T<100:
            R_fc = min(1,1-(T-Ti)*0.35/(100-Ti))
            R_ft = min(1,1-(T-Ti)*0.60/(100-Ti))
            R_Ec = min(1,1-(T-Ti)*0.65/(100-Ti))
            R_Et = min(1,1-(T-Ti)*0.50/(100-Ti))
            
        else:
            R_fc = max(1e-4,0.65-(T-100)*0.65/200)
            R_ft = max(1e-4,0.40-(T-100)*0.40/200)
            R_Ec = max(1e-4,0.35-(T-100)*0.35/200)
            R_Et = max(1e-4,0.50-(T-100)*0.50/200)
    
        # Update modulus of elasticity values
        rE_wc = self.E_wc * R_Ec
        rE_wt = self.E_wt * R_Et
        rE_q = self.E_q * R_Ec
        
        # Update strength values
        rf_cu = self.f_cu * R_fc
        rf_r = self.f_r * R_fc
        rf_tu = self.f_tu * R_ft
        
        # Stree and strain points (t = tension, c = compression)
        e0 = 0
        e1t = rf_tu/rE_wt
        e2t = e1t - 1e-2
        e1c = rf_cu/rE_wc
        e2c = e1c + (rf_r-rf_cu)/rE_q
        
        strain = np.array([e2t, e1t, e0, e1c, e2c])
        stress = np.array([rf_tu, rf_tu, 0, rf_cu, rf_r])
        
        return stress, strain
    
    def plot_stress_strain(self, T, Ti):
        '''
        Plot stress-strain curve at given surface and ambient temperatures.
        
        Args:
            T: An integer temperature at the surface, Celsius.
            Ti: An integer ambient temperature, Celsius.
        '''
        
        stress, strain = self.get_stress_strain(T, Ti)
        
        fig, ax = plt.subplots(figsize=(7,4), dpi=200)
        ax.grid(True, which='both', alpha=0.5)
        ax.axhline(y=0, color='k', lw=1)
        ax.axvline(x=0, color='k', lw=1)
        ax.set(xlabel='Strain ($\epsilon$), $10^{-4}$', 
               ylabel='Stress ($\sigma$), MPa')
        ax.plot(strain*1e4, stress, lw=1, c='k', ls='--')
        plt.show()