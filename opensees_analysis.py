#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import openseespy.opensees as ops
import openseespy.postprocessing.ops_vis as opsv
from openseespy.preprocessing import DiscretizeMember
import matplotlib.pyplot as plt
import numpy as np

from timber_section import Section
from timber_material import Material

class OpenSeesPyModel():
    '''
    Definition of the structural model in OpenSeesPy, analysed with 2nd-order 
    inelastic analysis.
    
    Attributes:
        section: Timber section object defined by width and height.
        material: Material object, defined by multiple material properties.
        nodes: List of node definitions [node_tag, coord_x, coord_y]
        elements: List of element definitions [mat_tag, node1, node2, n_disc]
        fixities: List of node fixities [node, tran_x, tran_y, rot_z]
        loads: List of point load definitions [node, Px, Py, Mz]
        
    Example usage:
        # Frame with point load at the midpoint of the beam
        nodes = [
            [1, 0, 0],
            [2, 0, 5000],
            [3, 5000, 5000],
            [4, 10000, 5000],
            [5, 10000, 0]]
        elements = [
            [0, 1, 2, 1],
            [1, 2, 3, 1],
            [2, 3, 4, 1],
            [3, 4, 5, 1]]
        fixities = [
            [1, 1, 1, 0],
            [5, 0, 1, 0]]
        loads = [
            [3, 0, -1000, 0]]
        node_recorders = [
            ['node2.txt', 2, 2, 'disp']]
        
        m = OpenSeesPyModel()
        m.define_section(200, 200, [5,5])
        m.define_fire_exposure(0, 300)
        m.section.plot(save='section.png')
        m.define_material()
        m.define_fibers()
        m.define_geometry(nodes, elements, fixities, num_integ=5)
        m.define_loads(loads)
        m.define_node_recorders(node_recorders)
        m.analyse(num_incr=1000)
        m.plot_deformed_shape(xlim=[-1000,12000], ylim=[-500,5000], scale=4, 
                              arrow_len=1, arrow_width=1.5, save='defl.png')
    '''

    def __init__(self, ndm=2, ndf=3):
        '''
        Initialise OpenSeesPy model.
        
        Args:
            ndm: An integer number of dimensions (default=2D).
            ndf: An integer number of degrees of freedom (default=3 DOF).
        '''
        
        ops.wipe()
        ops.model('basic', '-ndm', ndm, '-ndf', ndf)   
         
        self.material = None
        self.section = None
        self.nodes = None
        self.elements = None
        self.fixities = None
        self.loads = None
        
    def define_section(self, w, h, mesh_size, Ti=20, Tp=300):
        '''
        Define empty section.
        
        Args:
            w: An integer section width (z axis).
            h: An integer section height (y axis).
            mesh_size: A tuple mesh size in along z and y axes.
            Ti: An integer ambient temperature, Celsius (default=20C).
            Tp: An integer surface burning temperature, Celsius (default=300C).
        '''
        
        self.section = Section(w, h, Ti, Tp)
        self.section.discretise(mesh_size)
        ops.section('Fiber', 1)
        
    def define_fire_exposure(self, texp, T):
        '''
        Define fire exposure parameters and update the section
        
        Args:
            texp: A float time of exposure, min.
            T: An integer temperature at the surface, Celsius.
        '''
        
        self.section.set_exposure_time(texp)
        self.section.apply_temperature(T)
        
    def define_material(self, f_cu=45, f_tu=80, eps_y=4e-3, 
                        eps_r=13e-3, eps_tu=7.1e-3):
        '''
        Define a material object and assign appropriate properties to each 
        fiber in the section.
        
        Args:
            f_cu: A float ultimate compressive strength, MPa (default=45MPa).
            f_tu: A float ultimate tensile strength, MPa (default=80MPa).
            eps_y: A float yield compression strain (default=4e-3).
            eps_r: A float compressive strain at rupture (default=13e-3).
            eps_tu: A float tensile strain at rupture (default=7.1d-3).
        '''

        if self.section == None:
            raise Exception('No section is defined.')
            
        self.material = Material(f_cu, f_tu, eps_y, eps_r, eps_tu)
            
        for T, mattag in self.section.temp_dict.items():
            stress, strain = self.material.get_stress_strain(T, 
                                                             self.section.Ti)
            ops.uniaxialMaterial('ElasticMultiLinear', mattag, 0.0, 
                                 '-strain', *strain, '-stress', *stress)
            
    def define_fibers(self):
        '''
        Define fibers in the section.
        '''
        if self.material == None:
            raise Exception('No material is defined.')
        
        for p in self.section.patch_list:
            mattag, ny, nz, coord1, coord2 = p
            ops.patch('rect', mattag, ny, nz, *coord1, *coord2)
        
    def define_geometry(self, nodes, elements, fixities, num_integ=10):
        '''
        Define geometry of the structure (all dimensions are in mm).
        
        Args:
            nodes: A list of nodes in a form [node_tag, coord1, coord2].
            elements: A list of elements in a form 
                [ele_tag, node1, node2, disc].
            fixities: A list of fixities in a form [node, x, y, z].
            num_integ: Number of integration points along each element 
                (default=10)
        '''
        self.nodes = nodes
        self.elements = elements
        self.fixities = fixities
        
        if self.section == None:
            raise Exception('No section is defined.')
        
        ops.geomTransf('PDelta', 1)
        ops.beamIntegration('Lobatto', 1, 1, num_integ)
        
        for nd in self.nodes: ops.node(*nd)
            
        for el in self.elements:
            ele_tag, node1, node2, disc = el
            DiscretizeMember.DiscretizeMember(node1, node2, disc, 
                                              'forceBeamColumn', 1, 1,
                                              nodeTag=len(ops.getNodeTags())+1,
                                              eleTag=len(ops.getEleTags())+1)
        
        for fx in self.fixities: ops.fix(*fx)
        
    def define_loads(self, loads):
        '''
        Apply loads (all loads are in N)
        
        Args:
            loads: A list of point loads in a form [node, Px, Py, Pz].
        '''
        
        if self.nodes == None:
            raise Exception('No geometry is defined.')
            
        self.loads = loads
        
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        
        for ld in self.loads: ops.load(*ld)
        
    def define_node_recorders(self, node_recorders):
        '''
        Define node recorders that will track the specified results type for 
        each node.
        
        Args:
            node_recorder: A list of nodes to be tracked in the form 
                [recorder_name (str), node (int), dof(int), restype (str)].
                
                restype type can be seleceted from the following:
                    'disp' displacement,
                    'vel' velocity,
                    'accel' acceleration,
                    'incrDisp' incremental displacement,
                    'reaction' nodal reaction.
        '''
        
        for nr in node_recorders: 
            name, node, dof, restype = nr
            ops.recorder('Node','-file', name,'-closeOnWrite','-node', node, 
                         '-dof', dof, restype)
            
    def analyse(self, num_incr):
        '''
        Analyse the system.
        
        Args:
            num_incr: An integer number of load increments.
            return_node_disp: 
        '''

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
                
    def wipe_model(self):
        '''
        Clear the model after the analysis, required to output recorders.
        '''
        
        ops.wipe()

    def plot_structure(self):
        '''
        Plot nodes and elements of the model.
        '''
        
        fig, ax = plt.subplots(dpi=50)
        ax.set_axis_off()
        opsv.plot_model('nodes','elements')
        plt.show()
        
    def plot_deformed_shape(self, xlim, ylim, scale=1, arrow_len=10, 
                            arrow_width=2, save=''):
        '''
        Plot deformed shape of the model.
        
        Args:
            xlim: A list of left and right limits of the x axis.
            ylim: A list of bottom and top limits of the y axis.
            scale: A float scale of the displayed deformations (default=1).
            arrow_len: An integer length of the load arrows displayed 
                (default=10).
            arrow_width: An integer head width of the load arrows displayed 
                (default=2).
            save: A string indicating save path for the figure (default='',
                meaning that the figure will NOT be saved by default).
        '''
        
        fig, ax = plt.subplots(dpi=75)
        ax.set_axis_off()
        ax.grid(True, which='both', alpha=0.5)
        ax.axhline(y=0, color='k', lw=1)
        
        opsv.plot_defo(scale, fmt_undefo='k-', fmt_interp='k--')      
        
        ax.axis('equal')
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