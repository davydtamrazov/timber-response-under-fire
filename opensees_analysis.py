#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openseespy.opensees as ops
import numpy as np
from timber_section import Section


class OpenSeesPyAnalysis():

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
    def __init__(self, texp, w, h, L, Ti, Tp, a, mSize, e, P, num_incr,
                           side, plot=False, disp_fiber=True):
        self.section = Section(width=w, height=h, Ti=Ti, Tp=Tp)
        
        
    def reset(self):
        ops.wipe()
    ## Initialize the 2D model with three DOFs at each node
    ops.wipe()
    ops.model('basic','-ndm',2,'-ndf',3)
    
    ## Calculate temperature profile, define unique temperature tags and 
    ## group section patches
    [FiberList, MatList] = FiberSection.GetFibers(t, b, h, Ti, Tp, a, mSize, plot, disp_fiber, side)
    
    ## Define material properties
    fcu = 45
    fcr = fcu*0.75
    ftu = -80
    Ew = fcu/4e-3
    Eq = (fcu-fcr)/(4e-3-13e-3)
          
    # Define unique materials from temperature tags
    for T, mattag in MatList.items():
        
        # Update material properties given the temperature
        [s3t,s2t,s1t,s1c,s2c,s3c], [e3t,e2t,e1t,e1c,e2c,e3c] = StressStrain(T,Ti,fcu,fcr,ftu,Ew,Eq)
        
        # Define hysteretic material
        ops.uniaxialMaterial('Hysteretic',mattag,s1c,e1c,s2c,e2c,s3c,e3c,s1t,e1t,s2t,e2t,s3t,e3t,1,1,1,1)
        
    ## Define fiber section
    sectag = 1
    ops.section('Fiber',sectag)
    
    for fiber in FiberList:
        mattag, ny, nz, coord1, coord2 = fiber
        ops.patch('rect',mattag,ny,nz,*coord1,*coord2)
    
    ## Create nodes
    ops.node(1,0,0)
    ops.node(2,e,L/2)
    ops.node(3,0.0,L)
    
    ## Define fixities (pin-pin)
    ops.fix(1,1,1,0)
    ops.fix(3,1,0,0)
    
    ## Create fiber-based elements
    ops.geomTransf('PDelta',1)
    numIntgrPts = 10
    ops.beamIntegration('Lobatto',1,sectag,numIntgrPts)
    
    ops.element('forceBeamColumn',1,1,2,1,1)
    ops.element('forceBeamColumn',2,2,3,1,1)
    
      
    ## Apply loads
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain',1,1)
    ops.load(3,0,-P,0)
    
    # Load step increment (N)
    dP = 1/numIncr
    
    ## Define recorders
    ops.recorder('Node','-file','disp_node2.out','-node',2,'-dof',1,'disp')
    ops.recorder('Node','-file','reaction_node1.out','-node',1,'-dof',2,'reaction')
    ops.recorder('Node','-file','reaction_node3.out','-node',3,'-dof',2,'reaction')
    
    ## Set up analysis
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance',1e-6,100)
    ops.algorithm('Newton')
    ops.integrator('LoadControl',dP)
    ops.analysis('Static')
    
    ## Run analysis
    ops.record()
    ok = ops.analyze(numIncr)
    
    ## Report analysis status
    if ok == 0: print("Analysis done.")
    else:       print("Convergence issue.")
            
    ops.wipe()
    
    ## Compile output
    nodeDisp = abs(np.loadtxt('disp_node2.out'))
    reaction_node1 = np.loadtxt('reaction_node1.out')
    reaction_node3 = np.loadtxt('reaction_node3.out')
    
    appliedLoad = reaction_node1+reaction_node3
    
    return appliedLoad, nodeDisp


# Function to calculate temperature dependent material properties
#
# ----------------------------------------------------------------------------
# INPUT:
#   T           = temperature of the material (Celsius)
#   Ti          = ambient temperature (Celsius)
#   fcu         = ultimate compressive stress (MPa)
#   fcr         = compressive stress at rupture (MPa)
#   ftu         = ultimate tensile stress (MPa)
#   Ew          = modulus of elasticity (MPa)
#   Eq          = modulus of elasticity, softening (MPa)
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# OUTPUT:
#   stress      = list of stress points
#   strain      = list of strain points
# ----------------------------------------------------------------------------
#

def StressStrain(T,Ti,fcu,fcr,ftu,Ew,Eq):
    
    # Strength temperature reduction factors (Eurocode 5)
    kSC = min(1,1-(T-Ti)*0.35/(100-Ti)) if T<100 else max(0.0001,0.65-(T-100)*0.65/200)
    kST = min(1,1-(T-Ti)*0.60/(100-Ti)) if T<100 else max(0.0001,0.40-(T-100)*0.40/200)
        
    # Modulus of elasticity temperature reduction factors (Eurocode 5)
    kEC = min(1,1-(T-Ti)*0.65/(100-Ti)) if T<100 else max(0.0001,0.35-(T-100)*0.35/200)
    kET = min(1,1-(T-Ti)*0.50/(100-Ti)) if T<100 else max(0.0001,0.50-(T-100)*0.50/200)
    
    # Update modulus of elasticity values
    Ew_t = Ew*kET
    Ew_c = Ew*kEC
    Eq *= kEC
    
    # Update strength values
    fcu *= kSC
    fcr *= kSC
    ftu *= kST
        
    # Stress points (t = tension, c = compression)
    stress = [0,0,ftu,fcu,fcr,0]

    # Strain points (t = tension, c = compression)
    e1c = fcu/Ew_t
    e2c = e1c + (fcr-fcu)/Eq
    e3c = e2c + 1e-8
    
    e1t = ftu/Ew_c
    e2t = e1t - 1e-8
    e3t = e2t - 1e-8
    
    strain = [e3t,e2t,e1t,e1c,e2c,e3c]
    
    return stress,strain