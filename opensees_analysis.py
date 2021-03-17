#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openseespy.opensees as ops
from openseespy import test
import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsplt
from openseespy.preprocessing import DiscretizeMember
import matplotlib.pyplot as plt
import numpy as np
from timber_section import Section
from scipy.interpolate import make_interp_spline


class StructuralModel():

    '''
    # Function to perform 2nd-order inelastic analysis of the column in OpenSeesPy
    
    ----------------------------------------------------------------------------
    Args:
      t           = time of exposure
      b           = section width (mm)
      h           = section height (mm)
      L           = column height (mm)
      Ti          = ambient temperature (Celsius)
      Tp          = surface temperature (Celsius)
      a           = temperature penetration depth (mm)
      mSize       = mesh size (mm)
      P           = applied vertical load (N)
      numIncr     = number of load increments
      plot        = boolean to display section temperature profile plot
      disp_fiber  = boolean to display mesh on the plot
      side        = number of sides from which fire is applied (either 1 or 4)
    ----------------------------------------------------------------------------
    
    ----------------------------------------------------------------------------
    Returns:
      appliedLoad = list of the total applied loads at each load step
      nodeDisp    = list of the horizontal middle node displacements at each 
                    load step
    '''
    def __init__(self, texp, w, h, T, mesh_size, nodes, elements, fixs, loads, ndm=2, ndf=3):
                           
        self.section = Section(width=w, height=h)
        self.material = Material()
        self.nodes = nodes
        self.elements = elements
        self.fixs = fixs
        self.loads = loads
        
        self.section.discretize(mesh_size)
        self.section.set_exposure_time(texp)
        self.section.apply_temperature(T=T)
        
        ops.wipe()
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)
        
    def define_material(self):
        '''
        Define unique materials from temperature tags.
        '''
        
        # ops.uniaxialMaterial('Concrete01', 1, 50, 1e-3, 40, 1e-2)

        for T, mattag in self.section.temp_dict.items():
            stress, strain = self.material.get_stress_strain(T)
            ops.uniaxialMaterial('ElasticMultiLinear', mattag, 0.0, '-strain', *strain, '-stress', *stress)

    def define_fiber_section(self):
        '''
        Define fiber section.
        '''
        ops.section('Fiber', 1)
        for p in self.section.patch_list:
            mattag, ny, nz, coord1, coord2 = p
            ops.patch('rect', mattag, ny, nz, *coord1, *coord2)
        
    def define_geometry(self, num_integ=10):
        '''
        Define geometry of the structure.
        
        Args:
            nodes: A list of nodes in a form [node_tag, coord1, coord2]
            elements: A list of elements in a form [ele_tag, node1, node2, disc]
            fixities: A list of fixities in a form [node, x, y, z]
        '''
        
        ops.geomTransf('PDelta', 1)
        ops.beamIntegration('Lobatto', 1, 1, num_integ)
        
        for nd in self.nodes: ops.node(*nd)
            
        for el in self.elements:
            ele_tag, node1, node2, disc = el
            DiscretizeMember.DiscretizeMember(node1, node2, disc, 
                                              'forceBeamColumn', 1, 1,
                                              nodeTag=len(ops.getNodeTags())+1,
                                              eleTag=len(ops.getEleTags())+1)
        
        for fx in self.fixs: ops.fix(*fx)

        
    def define_loads(self):
        '''
        Apply loads
        
        Args:
            loads: A list of point loads (N) in a form [node, Px, Py, Pz]
        '''
        
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        
        for ld in self.loads: ops.load(*ld)
        
        
    def define_recorders(self, node_recorders):
        for nr in node_recorders: 
            name, node, dof, restype = nr
            ops.recorder('Node','-file', name,'-closeOnWrite','-node', node, 
                         '-dof', dof, restype)
            
    def analyse(self, num_incr, return_node_disp=None):
        '''
        Analyse the system.
        '''
        
        # Load step increment (N)

        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('BandGeneral')
        ops.test('NormUnbalance',2e-8, num_incr)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1/num_incr)
        ops.record()
        
        ops.analysis('Static')
        ok = ops.analyze(num_incr)

        # Report analysis status
        if ok == 0: print("Analysis done.")
        else:       print("Convergence issue.")
        
        # print(ops.nodeDisp(3))
        if return_node_disp==None:
            return ok
        else:
            return ok, ops.nodeDisp(return_node_disp)
                
    def wipe_model(self):
        ops.wipe()

    def plot_structure(self):
        fig, ax = plt.subplots(dpi=50)
        ax.set_axis_off()
        opsv.plot_model('nodes','elements')
        plt.show()
        
    def plot_deformed_shape(self, xlim, ylim, scale=1, arrow_len=10, arrow_width=2, save=0):
        
        fig, ax = plt.subplots(dpi=75)
        ax.set_axis_off()
        ax.grid(True, which='both', alpha=0.5)
        ax.axhline(y=0, color='k', lw=1)
        
        opsv.plot_defo(scale, fmt_undefo='k-', fmt_interp='k--')      
        
        ax.axis('equal')
        # ax.set_aspect('auto')
        ax.set(xlim=xlim, ylim=ylim)

        node_list = ops.getNodeTags()
        node_disp = np.array([ops.nodeDisp(n) for n in node_list])
        node_coord = np.array([ops.nodeCoord(n) for n in node_list])
        new_coord = node_disp[:,:-1]*scale + node_coord
        
        for node, Px, Py, M in self.loads:
            c = new_coord[node_list.index(node),:]
            ax.annotate('', xytext = (c[0]+abs(Px)*arrow_len, 
                                      c[1]+abs(Py)*arrow_len), 
                        xy = (c[0], c[1]), 
                        arrowprops = dict(arrowstyle=f'-|>, \
                                          head_width={arrow_width/5},\
                                          head_length={arrow_width/2}',
                                          lw=arrow_width,
                                          fc='orangered',
                                          ec='orangered'))
        
        if save: 
            fig.savefig(save, transparent=True)
            
        plt.show()
        
 
class Material():
    def __init__(self, f_cu=45, f_tu=80, eps_y=4e-3, eps_r=13e-3, eps_tu=7.1e-3):
        self.f_cu = f_cu # Ultimate strength
        self.f_tu = -f_tu # Tensile strength
        self.f_r = self.f_cu*0.75 # Rupture strength
        
        self.E_wc = self.f_cu/eps_y # Modulus of elasticity
        self.E_wt = self.E_wc
        self.E_tc = self.f_tu/eps_tu # Modulus of elasticity
        self.E_q = (self.f_r - self.f_cu)/(eps_r-eps_y) # Softening modulus of elasticity

    def get_stress_strain(self, T, Ti=20):
        '''
        Update strengths and moduli of elasticity with temperature reduction
        factors as per Eurocode 5 (EN 1995-1-2) and calculate stress-strain
        points.
        '''
        
        #array of strain points along stress-strain curve
        
        # Strength temperature reduction factors (Eurocode 5)
        R_fc = min(1,1-(T-Ti)*0.35/(100-Ti)) if T<100 else max(1e-4,0.65-(T-100)*0.65/200)
        R_ft = min(1,1-(T-Ti)*0.60/(100-Ti)) if T<100 else max(1e-4,0.40-(T-100)*0.40/200)
            
        # Modulus of elasticity temperature reduction factors (Eurocode 5)
        R_Ec = min(1,1-(T-Ti)*0.65/(100-Ti)) if T<100 else max(1e-4,0.35-(T-100)*0.35/200)
        R_Et = min(1,1-(T-Ti)*0.50/(100-Ti)) if T<100 else max(1e-4,0.50-(T-100)*0.50/200)
        
        # Update modulus of elasticity values
        rE_wc = self.E_wc * R_Ec
        rE_wt = self.E_wt * R_Et
        rE_q = self.E_q * R_Ec
        
        # Update strength values
        rf_cu = self.f_cu * R_fc
        rf_r = self.f_r * R_fc
        rf_tu = self.f_tu * R_ft
        
        stress = np.array([rf_tu, rf_tu, 0, rf_cu, rf_r])
        
        # Strain points (t = tension, c = compression)
        e0 = 0
        e1t = rf_tu/rE_wt
        e2t = e1t - 1e-2
        e1c = rf_cu/rE_wc
        e2c = e1c + (rf_r-rf_cu)/rE_q
        strain = np.array([e2t, e1t, e0, e1c, e2c])
        
        return stress, strain
    
    def plot(self, T):
        stress, strain = self.get_stress_strain(T)
        
        fig, ax = plt.subplots(figsize=(7,4), dpi=200)
        ax.grid(True, which='both', alpha=0.5)
        ax.axhline(y=0, color='k', lw=1)
        ax.axvline(x=0, color='k', lw=1)
        ax.set(xlabel='Strain ($\epsilon$), $10^{-4}$', ylabel='Stress ($\sigma$), MPa')
        ax.plot(strain*1e4, stress, lw=1, c='k', ls='--')