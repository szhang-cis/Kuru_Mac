from __future__ import print_function
#import os, sys, gc
#from time import time
#from copy import deepcopy
#import numpy as np
#import numpy.linalg as la
from warnings import warn
#from Florence import QuadratureRule, FunctionSpace
#from Florence.Base import JacobianError, IllConditionedError
#from Florence.Utils import PWD, RSWD

#from Florence.FunctionSpace import Tri
#from Florence.FunctionSpace import Tet
#from Florence.FunctionSpace import Quad, QuadES
#from Florence.FunctionSpace import Hex, HexES

from Florence.FiniteElements.LocalAssembly.KinematicMeasures import *
#from Florence.FiniteElements.LocalAssembly._KinematicMeasures_ import _KinematicMeasures_
#from Florence import Mesh
#from Florence.MeshGeneration import vtk_writer
#from Florence.Utils import constant_camera_view

class PostProcess(object):
    """Post-process class for finite element solvers"""

    def __init__(self,ndim,nvar):

        self.domain_bases = None
        self.postdomain_bases = None
        self.boundary_bases = None
        self.ndim = ndim
        self.nvar = nvar
        self.analysis_type = None
        self.analysis_nature = None
        self.material_type = None

        self.is_scaledjacobian_computed = False
        self.is_material_anisotropic = False
        self.directions = None

        self.mesh = None
        self.sol = None
        self.recovered_fields = None

        self.formulation = None
        self.material = None
        self.fem_solver = None

        self.parallel_model = None
        self.ncpu = None

    def SetBases(self,domain=None,postdomain=None,boundary=None):
        """Sets bases for all integration points for 'domain', 'postdomain' or 'boundary'
        """

        if domain is None and postdomain is None and boundary is None:
            warn("Nothing to be set")

        self.domain_bases = domain
        self.postdomain_bases = postdomain
        self.boundary_bases = boundary

    def SetAnalysis(self,analysis_type=None, analysis_nature=None):
        self.analysis_type = analysis_type
        self.analysis_nature = analysis_nature

    def SetMesh(self,mesh):
        """Set initial (undeformed) mesh"""
        self.mesh = mesh

    def SetSolution(self,sol):
        self.sol = sol

    def SetFormulation(self,formulation):
        self.formulation = formulation

    def SetMaterial(self,material):
        self.material = material

    def SetFEMSolver(self,fem_solver):
        self.fem_solver = fem_solver

    def QuantityNamer(self, num, print_name=True):
        """Returns the quantity (for augmented solution i.e. primary and recovered variables)
            name given its number (from numbering order)

            Variables (quantity) numbering:

                quantity    mechanics 2D    mechanics 3D    electro_mechanics 2D    electro_mechanics 3D
                ----------------------------------------------------------------------------------------
                0           ux              ux              ux                      ux
                ----------------------------------------------------------------------------------------
                1           uy              uy              uy                      uy
                ----------------------------------------------------------------------------------------
                2           F_xx            uz              phi                     uz
                ----------------------------------------------------------------------------------------
                3           F_xy            F_xx            F_xx                    phi
                ----------------------------------------------------------------------------------------
                4           F_yx            F_xy            F_xy                    F_xx
                ----------------------------------------------------------------------------------------
                5           F_yy            F_xz            F_yx                    F_xy
                ----------------------------------------------------------------------------------------
                6           H_xx            F_yx            F_yy                    F_xz
                ----------------------------------------------------------------------------------------
                7           H_xy            F_yy            H_xx                    F_yx
                ----------------------------------------------------------------------------------------
                8           H_yx            F_yz            H_xy                    F_yy
                ----------------------------------------------------------------------------------------
                9           H_yy            F_zx            H_yx                    F_yz
                ----------------------------------------------------------------------------------------
                10          J               F_zy            H_yy                    F_zx
                ----------------------------------------------------------------------------------------
                11          C_xx            F_zz            J                       F_zy
                ----------------------------------------------------------------------------------------
                12          C_xy            H_xx            C_xx                    F_zz
                ----------------------------------------------------------------------------------------
                13          C_yy            H_xy            C_xy                    H_xx
                ----------------------------------------------------------------------------------------
                14          G_xx            H_xz            C_yy                    H_xy
                ----------------------------------------------------------------------------------------
                15          G_xy            H_yx            G_xx                    H_xz
                ----------------------------------------------------------------------------------------
                16          G_yy            H_yy            G_xy                    H_yx
                ----------------------------------------------------------------------------------------
                17          detC            H_yz            G_yy                    H_yy
                ----------------------------------------------------------------------------------------
                18          S_xx            H_zx            detC                    H_yz
                ----------------------------------------------------------------------------------------
                19          S_xy            H_zy            S_xx                    H_zx
                ----------------------------------------------------------------------------------------
                20          S_yy            H_zz            S_xy                    H_zy
                ----------------------------------------------------------------------------------------
                21          p_hyd           J               S_yy                    H_zz
                ----------------------------------------------------------------------------------------
                22                          C_xx            E_x                     J
                ----------------------------------------------------------------------------------------
                23                          C_xy            E_y                     C_xx
                ----------------------------------------------------------------------------------------
                24                          C_xz            D_x                     C_xy
                ----------------------------------------------------------------------------------------
                25                          C_yy            D_y                     C_xz
                ----------------------------------------------------------------------------------------
                26                          C_yz            p_hyd                   C_yy
                ----------------------------------------------------------------------------------------
                27                          C_zz                                    C_yz
                ----------------------------------------------------------------------------------------
                28                          G_xx                                    C_zz
                ----------------------------------------------------------------------------------------
                29                          G_xy                                    G_xx
                ----------------------------------------------------------------------------------------
                30                          G_xz                                    G_xy
                ----------------------------------------------------------------------------------------
                31                          G_yy                                    G_xz
                ----------------------------------------------------------------------------------------
                32                          G_yz                                    G_yx
                ----------------------------------------------------------------------------------------
                33                          G_zz                                    G_yy
                ----------------------------------------------------------------------------------------
                34                          detC                                    G_zz
                ----------------------------------------------------------------------------------------
                35                          S_xx                                    detC
                ----------------------------------------------------------------------------------------
                36                          S_xy                                    S_xx
                ----------------------------------------------------------------------------------------
                37                          S_xz                                    S_xy
                ----------------------------------------------------------------------------------------
                38                          S_yy                                    S_xz
                ----------------------------------------------------------------------------------------
                39                          S_yz                                    S_yy
                ----------------------------------------------------------------------------------------
                40                          S_zz                                    S_yz
                ----------------------------------------------------------------------------------------
                41                          p_hyd                                   S_zz
                ----------------------------------------------------------------------------------------
                42                                                                  E_x
                ----------------------------------------------------------------------------------------
                43                                                                  E_y
                ----------------------------------------------------------------------------------------
                44                                                                  E_z
                ----------------------------------------------------------------------------------------
                45                                                                  D_x
                ----------------------------------------------------------------------------------------
                46                                                                  D_y
                ----------------------------------------------------------------------------------------
                47                                                                  D_z
                ----------------------------------------------------------------------------------------
                48                                                                  p_hyd
                ----------------------------------------------------------------------------------------




            where S represents Cauchy stress tensor, E the electric field, D the electric
            displacements and p_hyd the hydrostatic pressure
        """
        namer = None
        if num > 48:
            if print_name:
                print('Quantity corresponds to ' + str(namer))
            return namer

        lines = []
        with open(__file__) as f:
            lines.append(f.readlines())
        lines = lines[0]

        line_number = len(lines)+1

        for counter, line in enumerate(lines):
            line = line.strip()
            if "quantity" in line and "mechanics" in line and "2D" in line and "3D" in line:
                line_number = counter
            if counter > line_number+1 and counter < line_number+100:
                spl = list(filter(None, line.split(" ")))
                if spl[0] == str(num):
                    if self.nvar == 2 and self.ndim==2:
                        namer = spl[-4]
                    elif self.nvar == 3 and self.ndim==2:
                        namer = spl[-2]
                    elif self.nvar == 3 and self.ndim==3:
                        namer = spl[-3]
                    elif self.nvar == 4:
                        namer = spl[-1]
                    break

        if print_name:
            print('Quantity corresponds to ' + str(namer))

        if "ux" in namer:
            namer = "u_x"
        elif "uy" in namer:
            namer = "u_y"
        elif "uz" in namer:
            namer = "u_z"
        elif "phi" in namer:
            namer = "\phi"
        return namer

    def WriteVTK(self,filename=None, quantity="all", configuration="deformed", steps=None, write_curved_mesh=True,
        interpolation_degree=10, ProjectionFlags=None, fmt="binary", equally_spaced=False, parallelise=False):
        """Writes results to a VTK file for Paraview
            quantity = "all" means write all solution fields, otherwise specific quantities
            would be written based on augmented solution numbering order
            step - [list or np.1darray of sequentially aranged steps] which time steps/increments should be written
            inputs:
                fmt:                    [str] VTK writer format either "binary" or "xml".
                                        "xml" files do not support big vtk/vtu files
                                        typically greater than 2GB whereas "binary" files can.  Also "xml" writer is
                                        in-built whereas "binary" writer depends on evtk/pyevtk module
        """
        if fmt == "xml":
            pass
        elif fmt == "binary":
            try:
                from pyevtk.hl import pointsToVTK, linesToVTK, gridToVTK, unstructuredGridToVTK
                from pyevtk.vtk import VtkVertex, VtkLine, VtkTriangle, VtkQuad, VtkTetra, VtkPyramid, VtkHexahedron
            except ImportError:
                raise ImportError("Could not import evtk. Install it using 'pip install pyevtk'")
        else:
            raise ValueError("Writer format not understood")
        formatter = fmt

        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.sol is None:
            raise ValueError("solution not set for post-processing")
        if self.formulation.fields == "electrostatics":
            configuration = "original"
            tmp = np.copy(self.sol)
            self.sol = np.zeros((self.sol.shape[0],self.formulation.ndim+1,self.sol.shape[1]))
            # self.sol[:,:self.formulation.ndim,:] = 0.
            self.sol[:,-1,:] = tmp
            quantity = self.formulation.ndim

        if isinstance(quantity,int):
            if quantity>=self.sol.shape[1]:
                self.GetAugmentedSolution(parallelise=parallelise)
                if quantity >= self.sol.shape[1]:
                    raise ValueError('Plotting quantity not understood')
            iterator = range(quantity,quantity+1)
        elif isinstance(quantity,str):
            if quantity=="all":
                self.GetAugmentedSolution(parallelise=parallelise)
                iterator = range(self.sol.shape[1])
            else:
                raise ValueError('Plotting quantity not understood')
        elif isinstance(quantity,list):
            requires_augmented_solution = False
            for i in quantity:
                if i >= self.sol.shape[1]:
                    requires_augmented_solution = True
                    break
            if requires_augmented_solution:
                self.GetAugmentedSolution(parallelise=parallelise)
            iterator = quantity
        else:
            raise ValueError('Writing quantity not understood')

        if filename is None:
            warn("file name not specified. I am going to write in the current directory")
            filename = PWD(__file__) + "/output.vtu"
        elif filename is not None:
            if isinstance(filename,str) is False:
                raise ValueError("file name should be a string")
        if ".vtu" in filename and fmt is "binary":
            filename  = filename.split('.')[0]

        C = self.mesh.InferPolynomialDegree() - 1
        if C == 0:
            write_curved_mesh = False


        # GET LINEAR MESH & SOLUTION
        lmesh, sol = self.mesh.GetLinearMesh(remap=True,solution=self.sol)


        if lmesh.element_type =='tri':
            cellflag = 5
            offset = 3
            actual_ndim = 2
        elif lmesh.element_type =='quad':
            cellflag = 9
            offset = 4
            actual_ndim = 2
        if lmesh.element_type =='tet':
            cellflag = 10
            offset = 4
            actual_ndim = 3
        elif lmesh.element_type == 'hex':
            cellflag = 12
            offset = 8
            actual_ndim = 3
        actual_ndim = lmesh.points.shape[1]

        ndim = lmesh.points.shape[1]
        if self.formulation.fields == "electrostatics":
            sol = self.sol[:lmesh.nnode,...]
            q_names = ["phi","phi","phi","phi"]
        else:
            q_names = [self.QuantityNamer(quant, print_name=False) for quant in iterator]

        LoadIncrement = self.sol.shape[2]

        increments = range(LoadIncrement)
        if steps is not None:
            increments = steps

        if write_curved_mesh is False:
            parallelise = False
        if len(increments) == 1:
            parallelise = False

        # PARALLEL MODE
        if parallelise:
            from multiprocessing import Pool, cpu_count
            from contextlib import closing
            if self.ncpu is None:
                self.ncpu = cpu_count()
            if self.parallel_model is None:
                self.parallel_model = "pool"

            increments = np.array(increments).flatten()
            partitioned_steps = np.array_split(increments,self.ncpu)

            pps = []
            for ip in range(len(partitioned_steps)):
                pp = deepcopy(self)
                # pp.sol = self.sol[:,:,partitioned_steps[ip]]
                pp.number_of_load_increments = len(partitioned_steps[ip])
                zipper_object = ParallelVTKWriterZipper(pp,
                    filename=filename, quantity=quantity, configuration=configuration,
                    steps=partitioned_steps[ip], write_curved_mesh=write_curved_mesh,
                    interpolation_degree=interpolation_degree,
                    ProjectionFlags=ProjectionFlags, fmt=fmt,
                    equally_spaced=equally_spaced)
                pps.append(zipper_object)

            # Pool
            if self.parallel_model == "pool":
                with closing(Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelWriteVTK,pps)
                    pool.terminate()

            # Thread Pool
            elif self.parallel_model == "thread_pool":
                import multiprocessing.dummy
                with closing(multiprocessing.dummy.Pool(self.ncpu)) as pool:
                    res = pool.map(ParallelWriteVTK,pps)
                    pool.terminate()

            else:
                raise ValueError("Parallel model not understood")


            # # Serial
            # for ip in range(len(partitioned_steps)):
            #     ParallelWriteVTK(pps[ip])

            return



        if write_curved_mesh:

            if lmesh.element_type =='tet':
                cellflag = 5
                tmesh = PostProcess.TessellateTets(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_on_faces=False, plot_points=True,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags,
                    EquallySpacedPoints=equally_spaced)
            elif lmesh.element_type =='hex':
                cellflag = 5
                tmesh = PostProcess.TessellateHexes(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_on_faces=False, plot_points=True,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags,
                    EquallySpacedPoints=equally_spaced)
            elif lmesh.element_type =='quad':
                cellflag = 5
                tmesh = PostProcess.TessellateQuads(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_points=True, EquallySpacedPoints=equally_spaced,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags)
            elif lmesh.element_type =='tri':
                cellflag = 5
                tmesh = PostProcess.TessellateTris(self.mesh, np.zeros_like(self.mesh.points),
                    QuantityToPlot=self.sol[:,0,0], plot_points=True, EquallySpacedPoints=equally_spaced,
                    interpolation_degree=interpolation_degree, ProjectionFlags=ProjectionFlags)
            else:
                raise ValueError('Element type not understood')

            nsize = tmesh.nsize
            if hasattr(tmesh,'nface'):
                # FOR 3D ELEMENTS E.G. TETS AND HEXES
                nface = tmesh.nface
            else:
                tmesh.smesh = self.mesh
                tmesh.faces_to_plot = tmesh.smesh.elements
                nface = tmesh.smesh.elements.shape[0]
                tmesh.smesh.GetEdges()

                connections_elements = np.arange(tmesh.x_edges.size).reshape(tmesh.x_edges.shape[1],tmesh.x_edges.shape[0])
                connections = np.zeros((tmesh.x_edges.size,2),dtype=np.int64)
                for i in range(connections_elements.shape[0]):
                    connections[i*(tmesh.x_edges.shape[0]-1):(i+1)*(tmesh.x_edges.shape[0]-1),0] = connections_elements[i,:-1]
                    connections[i*(tmesh.x_edges.shape[0]-1):(i+1)*(tmesh.x_edges.shape[0]-1),1] = connections_elements[i,1:]
                connections = connections[:(i+1)*(tmesh.x_edges.shape[0]-1),:]
                tmesh.connections = connections

            un_faces_to_plot = np.unique(tmesh.faces_to_plot)
            fail_flag = False
            try:
                ssol = self.sol[un_faces_to_plot,:,:]
            except:
                fail_flag = True

            if fail_flag is False:
                if tmesh.smesh.elements.max() > un_faces_to_plot.shape[0]:
                    ssol = self.sol
                    fail_flag = True
                    warn("Something went wrong with mesh tessellation for VTK writer. I will proceed anyway")

            if tmesh.smesh.all_edges.shape[0] > tmesh.edge_elements.shape[0]:
                tmesh.smesh.all_edges = tmesh.edge_elements
                fail_flag = True
                warn("Something went wrong with mesh tessellation for VTK writer. I will proceed anyway")


            for Increment in increments:

                extrapolated_sol = np.zeros((tmesh.points.shape[0], self.sol.shape[1]))
                for ielem in range(nface):
                    extrapolated_sol[ielem*nsize:(ielem+1)*nsize,:] = np.dot(tmesh.bases_2,
                        ssol[tmesh.smesh.elements[ielem,:],:, Increment])

                if not fail_flag:
                    svpoints = self.mesh.points[np.unique(tmesh.faces_to_plot),:] + ssol[:,:tmesh.points.shape[1],Increment]
                else:
                    svpoints = self.mesh.points + ssol[:,:tmesh.points.shape[1],Increment]

                for iedge in range(tmesh.smesh.all_edges.shape[0]):
                    ielem = tmesh.edge_elements[iedge,0]
                    edge = tmesh.smesh.elements[ielem,tmesh.reference_edges[tmesh.edge_elements[iedge,1],:]]
                    coord_edge = svpoints[edge,:]
                    if tmesh.points.shape[1] == 3:
                        tmesh.x_edges[:,iedge], tmesh.y_edges[:,iedge], tmesh.z_edges[:,iedge] = np.dot(coord_edge.T,tmesh.bases_1)
                    elif tmesh.points.shape[1] == 2:
                        tmesh.x_edges[:,iedge], tmesh.y_edges[:,iedge] = np.dot(coord_edge.T,tmesh.bases_1)

                if tmesh.points.shape[1] == 3:
                    edge_coords = np.concatenate((tmesh.x_edges.T.copy().flatten()[:,None],
                        tmesh.y_edges.T.copy().flatten()[:,None],
                        tmesh.z_edges.T.copy().flatten()[:,None]),axis=1)
                elif tmesh.points.shape[1] == 2:
                    edge_coords = np.concatenate((tmesh.x_edges.T.copy().flatten()[:,None],
                        tmesh.y_edges.T.copy().flatten()[:,None], np.zeros_like(tmesh.y_edges.T.copy().flatten()[:,None])),axis=1)
                    svpoints = np.concatenate((svpoints, np.zeros((svpoints.shape[0],1))),axis=1)

                if formatter == "xml":

                    vtk_writer.write_vtu(Verts=edge_coords,
                        Cells={3:tmesh.connections},
                        fname=filename.split('.')[0]+'_curved_lines_increment_'+str(Increment)+'.vtu')

                    vtk_writer.write_vtu(Verts=svpoints,
                        Cells={1:np.arange(svpoints.shape[0])},
                        fname=filename.split('.')[0]+'_curved_points_increment_'+str(Increment)+'.vtu')

                    for quant in iterator:
                        vtk_writer.write_vtu(Verts=tmesh.points+extrapolated_sol[:,:ndim],
                            Cells={cellflag:tmesh.elements}, pdata=extrapolated_sol[:,quant],
                            fname=filename.split('.')[0]+'_curved_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')

                elif formatter == "binary":

                    unstructuredGridToVTK(filename.split('.')[0]+'_curved_lines_increment_'+str(Increment),
                        np.ascontiguousarray(edge_coords[:,0]), np.ascontiguousarray(edge_coords[:,1]), np.ascontiguousarray(edge_coords[:,2]),
                        np.ascontiguousarray(tmesh.connections.ravel()), np.arange(0,2*tmesh.connections.shape[0],2)+2,
                        np.ones(tmesh.connections.shape[0])*3)

                    pointsToVTK(filename.split('.')[0]+'_curved_points_increment_'+str(Increment),
                        np.ascontiguousarray(svpoints[:,0]), np.ascontiguousarray(svpoints[:,1]), np.ascontiguousarray(svpoints[:,2]),
                        data=None)

                    if tmesh.points.shape[1] == 2:
                        points = np.zeros((tmesh.points.shape[0],3))
                        points[:,:2] = tmesh.points+extrapolated_sol[:,:ndim]
                    else:
                        points = tmesh.points+extrapolated_sol[:,:ndim]
                    for counter, quant in enumerate(iterator):
                        unstructuredGridToVTK(filename.split('.')[0]+'_curved_quantity_'+str(quant)+'_increment_'+str(Increment),
                            np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                            np.ascontiguousarray(tmesh.elements.ravel()), np.arange(0,3*tmesh.elements.shape[0],3)+3,
                            np.ones(tmesh.elements.shape[0])*cellflag,
                            pointData={q_names[counter]: np.ascontiguousarray(extrapolated_sol[:,quant])})

        else:

            if configuration == "original":
                for Increment in increments:
                    if formatter == "xml":
                        for quant in iterator:
                            vtk_writer.write_vtu(Verts=lmesh.points,
                                Cells={cellflag:lmesh.elements}, pdata=sol[:,quant,Increment],
                                fname=filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')
                    elif formatter == "binary":
                        # points = lmesh.points
                        if lmesh.InferSpatialDimension() == 2:
                            points = np.zeros((lmesh.points.shape[0],3))
                            points[:,:2] = lmesh.points
                        else:
                            points = lmesh.points
                        for counter, quant in enumerate(iterator):
                            unstructuredGridToVTK(filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment),
                                np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                                np.ascontiguousarray(lmesh.elements.ravel()), np.arange(0,offset*lmesh.nelem,offset)+offset,
                                np.ones(lmesh.nelem)*cellflag,
                                pointData={q_names[counter]: np.ascontiguousarray(sol[:,quant,Increment])})

            elif configuration == "deformed":
                for Increment in increments:
                    if formatter == "xml":
                        for quant in iterator:
                            vtk_writer.write_vtu(Verts=lmesh.points+sol[:,:ndim,Increment],
                                Cells={cellflag:lmesh.elements}, pdata=sol[:,quant,Increment],
                                fname=filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment)+'.vtu')
                    elif formatter == "binary":
                        if lmesh.InferSpatialDimension() == 2:
                            points = np.zeros((lmesh.points.shape[0],3))
                            points[:,:2] = lmesh.points + sol[:,:ndim,Increment]
                        else:
                            points = lmesh.points + sol[:,:ndim,Increment]

                        for counter, quant in enumerate(iterator):
                            unstructuredGridToVTK(filename.split('.')[0]+'_quantity_'+str(quant)+'_increment_'+str(Increment),
                                np.ascontiguousarray(points[:,0]), np.ascontiguousarray(points[:,1]), np.ascontiguousarray(points[:,2]),
                                np.ascontiguousarray(lmesh.elements.ravel()), np.arange(0,offset*lmesh.nelem,offset)+offset,
                                np.ones(lmesh.nelem)*cellflag,
                                pointData={q_names[counter]: np.ascontiguousarray(sol[:,quant,Increment])})

        return
