#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np
from timber_section import Section
from matplotlib.patches import Rectangle


# import OpenSeesPy plotting commands
import openseespy.postprocessing.Get_Rendering as opsplt


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
    def __init__(self, texp, w, h, mesh_size, ndm=2, ndf=3):
                           
        self.section = Section(width=w, height=h)
        self.material = Material()
        
        ops.wipe()
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)
        self.section.discretize(mesh_size)
        self.section.set_exposure_time(0)
        self.section.apply_temperature()
        
    def define_material(self):
        '''
        Define unique materials from temperature tags.
        '''
        
        # ops.uniaxialMaterial('Concrete01', 1, 50, 1e-3, 40, 1e-2)

        for T, mattag in self.section.temp_dict.items():
            stress, strain = self.material.get_stress_strain(T)
            ops.uniaxialMaterial('ElasticMultiLinear', mattag, 0.0, '-strain', *strain, '-stress', *stress)

        
        # opsplt.plot_model()
            # ops.uniaxialMaterial('MultiLinear', mattag, *[i for p in pairs for i in p])
        # Define hysteretic material
        # ops.uniaxialMaterial('Hysteretic',mattag,s1c,e1c,s2c,e2c,s3c,e3c,s1t,e1t,s2t,e2t,s3t,e3t,1,1,1,1)
        
    def define_section(self):
        '''
        Define fiber section.
        '''
        ops.section('Fiber', 1)
        for p in self.section.patch_list:
            mattag, ny, nz, coord1, coord2 = p
            ops.patch('rect', mattag, ny, nz, *coord1, *coord2)
        
    def define_geometry(self, nodes, elements, fixities):
        '''
        Define geometry of the structure.
        
        Args:
            nodes: A list of nodes in a form [node_tag, coord1, coord2]
            elements: A list of elements in a form [elm_tag, node1, node2]
            fixities: A list of fixities in a form [node, x, y, z]
        '''
        
        ops.geomTransf('PDelta', 1)
        ops.beamIntegration('Lobatto', 1, 1, 10)
        
        for nd in nodes: ops.node(*nd)
        for el in elements: ops.element('forceBeamColumn', *el, 1, 1)
        for fx in fixities: ops.fix(*fx)
        
    def define_loads(self, loads):
        '''
        Apply loads
        
        Args:
            loads: A list of point loads (N) in a form [node, Px, Py, Pz]
        '''
        
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        
        for ld in loads: ops.load(*ld)
        
        
    def define_recorders(self, node_recorders):
        for nr in node_recorders: 
            name, node, dof, restype = nr
            ops.recorder('Node','-file', name,'-node', node, 
                         '-dof', dof, restype)
            
    def analyse(self, num_incr):
        '''
        Analyse the system.
        '''
        
        # Load step increment (N)

        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('BandGeneral')
        ops.test('NormUnbalance',1e-8, 10)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1/num_incr)
        ops.analysis('Static')
        
        ops.record()
        ok = ops.analyze(num_incr)
        
        print(ops.nodeDisp(3,2), ops.getLoadFactor(1))
        # perform the analysis
        
        # data = np.zeros((num_incr+1,2))
        # for j in range(num_incr):
        #     ops.analyze(1)
        #     data[j+1,0] = ops.nodeDisp(3,2)
        #     print(ops.getLoadFactor(1))
        #     data[j+1,1] = ops.getLoadFactor(1)*P
        
        # plt.plot(data[:,0], data[:,1])
        # plt.xlabel('Horizontal Displacement')
        # plt.ylabel('Horizontal Load')
        # plt.show()
        
        ## Report analysis status
        if ok == 0: print("Analysis done.")
        else:       print("Convergence issue.")

    def plot(self):
        opsplt.plot_model('nodes','elements')
        
    
        
        
 
class Material():
    def __init__(self, f_cu=45, f_tu=80, eps_y=4e-3, eps_r=13e-3, eps_tu=7.1e-3):
        self.f_cu = f_cu # Ultimate strength
        self.f_tu = -f_tu # Tensile strength
        self.f_r = self.f_cu*0.75 # Rupture strength
        
        # self.rf_tu = self.f_tu
        # self.rf_cu = self.f_cu
        # self.rf_r = self.f_r
        
        self.E_wc = self.f_cu/eps_y # Modulus of elasticity
        self.E_wt = self.E_wc
        self.E_tc = self.f_tu/eps_tu # Modulus of elasticity
        self.E_q = (self.f_r - self.f_cu)/(eps_r-eps_y) # Softening modulus of elasticity
        
        # self.rE_wc = self.E_wc
        # self.rE_wt = self.E_wt
        # self.rE_q = self.E_q

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
    
        
m = StructuralModel(texp=0, w=100, h=200, mesh_size=[5,5])
m.define_material()
m.define_section()

nodes = [[1, 0, 0],
         [2, 0, 5000],
         [3, 5000, 5000],
         [4, 10000, 5000],
         [5, 10000, 0]]

elements = [[0, 1, 2],
            [1, 2, 3],
            [2, 3, 4],
            [3, 4, 5]]

fixities = [[1, 1, 1, 1],
            [5, 1, 1, 1]]

m.define_geometry(nodes, elements, fixities)

loads = [[3, 0, -1000, 0]]
m.define_loads(loads)

recorders = [['disp_node3.txt', 3, 2, 'disp']]
m.define_recorders(recorders)

m.material.plot(0)
m.section.plot()
m.plot()

m.analyse(100)
m = None

# # c = OpenSeesPyAnalysis(0, 100,100,1000,20,300,35,5,0.5,100,10,[1,1,1,1]) 
# # c.define_model()
  
#     ## Calculate temperature profile, define unique temperature tags and 
#     ## group section patches
#     [FiberList, MatList] = FiberSection.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, side)
    
#     ## Define material properties

          
#     # Define unique materials from temperature tags
#     for T, mattag in MatList.items():
        
#         # Update material properties given the temperature
#         [s3t,s2t,s1t,s1c,s2c,s3c], [e3t,e2t,e1t,e1c,e2c,e3c] = StressStrain(T,Ti,fcu,fcr,ftu,Ew,Eq)
        
#         # Define hysteretic material
#         ops.uniaxialMaterial('Hysteretic',mattag,s1c,e1c,s2c,e2c,s3c,e3c,s1t,e1t,s2t,e2t,s3t,e3t,1,1,1,1)
        
#     ## Define fiber section
#     sectag = 1
#     ops.section('Fiber',sectag)
    
#     for fiber in FiberList:
#         mattag, ny, nz, coord1, coord2 = fiber
#         ops.patch('rect',mattag,ny,nz,*coord1,*coord2)
    
#     ## Create nodes
#     ops.node(1,0,0)
#     ops.node(2,e,L/2)
#     ops.node(3,0.0,L)
    
#     ## Define fixities (pin-pin)
#     ops.fix(1,1,1,0)
#     ops.fix(3,1,0,0)
    
#     ## Create fiber-based elements
#     ops.geomTransf('PDelta',1)
#     numIntgrPts = 10
#     ops.beamIntegration('Lobatto',1,sectag,numIntgrPts)
    
#     ops.element('forceBeamColumn',1,1,2,1,1)
#     ops.element('forceBeamColumn',2,2,3,1,1)
    
      
#     ## Apply loads
#     ops.timeSeries('Linear', 1)
#     ops.pattern('Plain',1,1)
#     ops.load(3,0,-P,0)
    
#     # Load step increment (N)
#     dP = 1/numIncr
    
#     ## Define recorders
#     ops.recorder('Node','-file','disp_node2.out','-node',2,'-dof',1,'disp')
#     ops.recorder('Node','-file','reaction_node1.out','-node',1,'-dof',2,'reaction')
#     ops.recorder('Node','-file','reaction_node3.out','-node',3,'-dof',2,'reaction')
    
#     ## Set up analysis
#     ops.constraints('Transformation')
#     ops.numberer('RCM')
#     ops.system('BandGeneral')
#     ops.test('NormUnbalance',1e-6,100)
#     ops.algorithm('Newton')
#     ops.integrator('LoadControl',dP)
#     ops.analysis('Static')
    
#     ## Run analysis
#     ops.record()
#     ok = ops.analyze(numIncr)
    
#     ## Report analysis status
#     if ok == 0: print("Analysis done.")
#     else:       print("Convergence issue.")
            
#     ops.wipe()
    
#     ## Compile output
#     nodeDisp = abs(np.loadtxt('disp_node2.out'))
#     reaction_node1 = np.loadtxt('reaction_node1.out')
#     reaction_node3 = np.loadtxt('reaction_node3.out')
    
#     appliedLoad = reaction_node1+reaction_node3
    
#     return appliedLoad, nodeDisp


# # Function to calculate temperature dependent material properties
# #
# # ----------------------------------------------------------------------------
# # INPUT:
# #   T           = temperature of the material (Celsius)
# #   Ti          = ambient temperature (Celsius)
# #   fcu         = ultimate compressive stress (MPa)
# #   fcr         = compressive stress at rupture (MPa)
# #   ftu         = ultimate tensile stress (MPa)
# #   Ew          = modulus of elasticity (MPa)
# #   Eq          = modulus of elasticity, softening (MPa)
# # ----------------------------------------------------------------------------
# #
# # ----------------------------------------------------------------------------
# # OUTPUT:
# #   stress      = list of stress points
# #   strain      = list of strain points
# # ----------------------------------------------------------------------------
# #
