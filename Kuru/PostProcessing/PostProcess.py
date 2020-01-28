from __future__ import print_function
import os, sys, gc
#from time import time
#from copy import deepcopy
#import numpy as np
#import numpy.linalg as la
from warnings import warn
from Kuru import FunctionSpace, QuadratureRule
#from Florence.Base import JacobianError, IllConditionedError
#from Florence.Utils import PWD, RSWD

#from Florence.FunctionSpace import Tri
#from Florence.FunctionSpace import Tet
#from Florence.FunctionSpace import Quad, QuadES
#from Florence.FunctionSpace import Hex, HexES

from Kuru.FiniteElements.LocalAssembly.KinematicMeasures import *
from Kuru.FiniteElements.LocalAssembly._KinematicMeasures_ import _KinematicMeasures_
from Kuru import Mesh
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

    def StressRecovery(self, steps=None, average_derived_quantities=True):

        """
            steps:          [list,np.1darray] for which time steps/increments the data should
                            be recovered
        """

        if self.mesh is None:
            raise ValueError("Mesh not set for post-processing")
        if self.sol is None:
            raise ValueError("Solution not set for post-processing")
        if self.formulation is None:
            raise ValueError("formulation not set for post-processing")
        if self.material is None:
            raise ValueError("material not set for post-processing")
        if self.fem_solver is None:
            raise ValueError("FEM solver not set for post-processing")

        if self.sol.shape[1] > self.nvar:
            return

        det = np.linalg.det
        inv = np.linalg.inv

        mesh = self.mesh
        fem_solver = self.fem_solver
        formulation = self.formulation
        material = self.material

        # GET THE UNDERLYING LINEAR MESH
        # lmesh = mesh.GetLinearMesh()
        C = mesh.InferPolynomialDegree() - 1
        ndim = mesh.InferSpatialDimension()

        elements = mesh.elements
        points = mesh.points
        nelem = elements.shape[0]; npoint = points.shape[0]
        nodeperelem = elements.shape[1]

        # GET QUADRATURE
        norder = 2*C
        if norder == 0:
            norder=1
        # quadrature = QuadratureRule(qtype="gauss", norder=norder, mesh_type=mesh.element_type, optimal=3)
        # Domain = FunctionSpace(mesh, quadrature, p=C+1)
        Domain = FunctionSpace(mesh, p=C+1, evaluate_at_nodes=True)
        Jm = Domain.Jm
        AllGauss = Domain.AllGauss
        Bases = Domain.Bases

        # requires_geometry_update = fem_solver.requires_geometry_update
        requires_geometry_update = True # ALWAYS TRUE FOR THIS ROUTINE
        TotalDisp = self.sol[:,:]

        if hasattr(self,"number_of_load_increments"):
            LoadIncrement = self.number_of_load_increments
        else:
            LoadIncrement = fem_solver.number_of_load_increments
        increments = range(LoadIncrement)
        if steps is not None:
            LoadIncrement = len(steps)
            increments = steps

        # COMPUTE THE COMMON/NEIGHBOUR NODES ONCE
        all_nodes = np.unique(elements)
        Elss, Poss = mesh.GetNodeCommonality()[:2]

        F = np.zeros((nelem,nodeperelem,ndim,ndim))
        CauchyStressTensor = np.zeros((nelem,nodeperelem,ndim,ndim))
        # DEFINE CONSTITUENT STRESSES FOR GROWTH-REMODELING PROBLEM
        FibreStress = np.zeros((nelem,nodeperelem,5))  # 5-fibres
        Softness = np.zeros((nelem,nodeperelem,5))  # 5-fibres

        MainDict = {}
        MainDict['F'] = np.zeros((LoadIncrement,npoint,ndim,ndim))
        MainDict['CauchyStress'] = np.zeros((LoadIncrement,npoint,ndim,ndim))
        MainDict['FibreStress'] = np.zeros((LoadIncrement,npoint,5))
        MainDict['Softness'] = np.zeros((LoadIncrement,npoint,5))

        for incr, Increment in enumerate(increments):
            if TotalDisp.ndim == 3:
                Eulerx = points + TotalDisp[:,:ndim,Increment]
            else:
                Eulerx = points + TotalDisp[:,:ndim]

            # LOOP OVER ELEMENTS
            for elem in range(nelem):
                # GET THE FIELDS AT THE ELEMENT LEVEL
                LagrangeElemCoords = points[elements[elem,:],:]
                EulerELemCoords = Eulerx[elements[elem,:],:]
                # GROWTH-REMODELING VALUES FOR THIS ELEMENT
                material.MappingFieldVariables(mesh,Domain,elem)

                if material.has_low_level_dispatcher:

                    # GET LOCAL KINEMATICS
                    SpatialGradient, F[elem,:,:,:], detJ, dV = _KinematicMeasures_(Jm, AllGauss[:,0],
                            LagrangeElemCoords, EulerELemCoords, requires_geometry_update)
                    # PARAMETERS FOR INCOMPRESSIBILITY (MEAN DILATATION METHOD HU-WASHIZU)
                    if material.is_incompressible:
                        stiffness_k, pressure = _VolumetricStiffnessIntegrand_(SpatialGradient, detJ, dV, formulation.nvar)
                        material.pressure = material.kappa*pressure

                    if self.formulation.fields == "electro_mechanics":
                        # GET ELECTRIC FIELD
                        ElectricFieldx[elem,:,:] = - np.einsum('ijk,j',SpatialGradient,ElectricPotentialElem)
                        # COMPUTE WORK-CONJUGATES AND HESSIAN AT THIS GAUSS POINT
                        _D_dum ,CauchyStressTensor[elem,:,:], _ = material.KineticMeasures(F[elem,:,:,:],
                                ElectricFieldx[elem,:,:], elem=elem)
                        ElectricDisplacementx[elem,:,:] = _D_dum[:,:,0]
                    elif self.formulation.fields == "mechanics":
                        # COMPUTE WORK-CONJUGATES AND HESSIAN AT THIS GAUSS POINT
                        CauchyStressTensor[elem,:,:], _ = material.KineticMeasures(F[elem,:,:,:],elem=elem)
                        if material.has_field_variables:
                            FibreStress[elem,:,:],Softness[elem,:,:] = material.LLConstituentStress(F[elem,:,:,:],elem)

                else:
                    # GAUSS LOOP IN VECTORISED FORM
                    ParentGradientX = np.einsum('ijk,jl->kil', Jm, LagrangeElemCoords)
                    # MATERIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla_0 (N)]
                    MaterialGradient = np.einsum('ijk,kli->ijl', inv(ParentGradientX), Jm)
                    # DEFORMATION GRADIENT TENSOR [\vec{x} \otimes \nabla_0 (N)]
                    F[elem,:,:,:] = np.einsum('ij,kli->kjl', EulerELemCoords, MaterialGradient)
                    # COMPUTE REMAINING KINEMATIC MEASURES
                    StrainTensors = KinematicMeasures(F[elem,:,:,:], fem_solver.analysis_nature)

                    # GEOMETRY UPDATE IS A MUST
                    # MAPPING TENSOR [\partial\vec{X}/ \partial\vec{\varepsilon} (ndim x ndim)]
                    ParentGradientx = np.einsum('ijk,jl->kil',Jm,EulerELemCoords)
                    # SPATIAL GRADIENT TENSOR IN PHYSICAL ELEMENT [\nabla (N)]
                    SpatialGradient = np.einsum('ijk,kli->ilj',inv(ParentGradientx),Jm)
                    # COMPUTE ONCE detJ (GOOD SPEEDUP COMPARED TO COMPUTING TWICE)
                    detJ = np.einsum('i,i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)),
                        np.abs(StrainTensors['J']))

                    # COMPUTE PARAMETERS FOR MEAN DILATATION METHOD, IT NEEDS TO BE BEFORE COMPUTE HESSIAN AND STRESS
                    if material.is_incompressible:
                        dVolume = np.einsum('i,i->i',AllGauss[:,0],np.abs(det(ParentGradientX)))
                        MaterialVolume, CurrentVolume = 0.0, 0.0
                        for i in range(AllGauss.shape[0]):
                            CurrentVolume += detJ[i]
                            MaterialVolume += dVolume[i]

                        material.pressure = material.kappa*(CurrentVolume-MaterialVolume)/MaterialVolume

                    # LOOP OVER GAUSS POINTS
                    for counter in range(AllGauss.shape[0]):

                        if material.energy_type == "enthalpy":
                            # COMPUTE CAUCHY STRESS TENSOR
                            CauchyStressTensor[elem,counter,:] = material.CauchyStress(StrainTensors,elem,counter)
                            if material.has_field_variables:
                                FibreStress[elem,counter,:],Softness[elem,counter,:] = material.ConstituentStress(
                                    StrainTensors,elem,counter)

                        elif material.energy_type == "internal_energy":
                            # COMPUTE CAUCHY STRESS TENSOR
                            CauchyStressTensor[elem,counter,:] = material.CauchyStress(StrainTensors,elem,counter)
                            if material.has_field_variables:
                                FibreStress[elem,counter,:],Softness[elem,counter,:] = material.ConstituentStress(
                                    StrainTensors,elem,counter)


            if average_derived_quantities:
                for inode in all_nodes:
                    # Els, Pos = np.where(elements==inode)
                    Els, Pos = Elss[inode], Poss[inode]
                    ncommon_nodes = Els.shape[0]
                    for uelem in range(ncommon_nodes):
                        MainDict['F'][incr,inode,:,:] += F[Els[uelem],Pos[uelem],:,:]
                        MainDict['CauchyStress'][incr,inode,:,:] += CauchyStressTensor[Els[uelem],Pos[uelem],:,:]
                        MainDict['FibreStress'][incr,inode,:] += FibreStress[Els[uelem],Pos[uelem],:]
                        MainDict['Softness'][incr,inode,:] += Softness[Els[uelem],Pos[uelem],:]

                    # AVERAGE OUT
                    MainDict['F'][incr,inode,:,:] /= ncommon_nodes
                    MainDict['CauchyStress'][incr,inode,:,:] /= ncommon_nodes
                    MainDict['FibreStress'][incr,inode,:] /= ncommon_nodes
                    MainDict['Softness'][incr,inode,:] /= ncommon_nodes

            else:
                for inode in all_nodes:
                    # Els, Pos = np.where(elements==inode)
                    Els, Pos = Elss[inode], Poss[inode]
                    ncommon_nodes = Els.shape[0]
                    uelem = 0
                    MainDict['F'][incr,inode,:,:] = F[Els[uelem],Pos[uelem],:,:]
                    MainDict['CauchyStress'][incr,inode,:,:] = CauchyStressTensor[Els[uelem],Pos[uelem],:,:]
                    MainDict['FibreStress'][incr,inode,:] = FibreStress[Els[uelem],Pos[uelem],:]
                    MainDict['Softness'][incr,inode,:] = Softness[Els[uelem],Pos[uelem],:]


        self.recovered_fields = MainDict
        return

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

    @staticmethod
    def TessellateHexes(mesh, TotalDisp, QuantityToPlot=None, plot_on_faces=True,
        ProjectionFlags=None, interpolation_degree=20, EquallySpacedPoints=False,
        plot_points=False, plot_edges=True, plot_surfaces=True):

        """High order curved hexahedral surfaces mesh plots, based on high order nodal FEM.
            The equally spaced FEM points do not work as good as the Fekete points
        """



        from Kuru.QuadratureRules import GaussLobattoPointsQuad
        from Kuru.QuadratureRules.NumericIntegrator import GaussLobattoQuadrature
        from Kuru.QuadratureRules.EquallySpacedPoints import EquallySpacedPoints as EquallySpacedPointsnD
        from Kuru.MeshGeneration.NodeArrangement import NodeArrangementQuad
        from Kuru.FunctionSpace import Quad
        from Kuru.FunctionSpace.OneDimensional.Line import LagrangeGaussLobatto, Lagrange

        from copy import deepcopy
        from scipy.spatial import Delaunay

        assert mesh.element_type == "hex"

        # SINCE THIS IS A 3D PLOT
        ndim=3

        C = interpolation_degree
        p = C+1
        nsize = int((C+2)**2)
        CActual = mesh.InferPolynomialDegree() - 1
        nsize_2 = int((CActual+2)**2)

        GaussLobattoPoints = GaussLobattoPointsQuad(C)
        if EquallySpacedPoints:
            GaussLobattoPoints = EquallySpacedPointsnD(ndim,C)

        # BUILD DELAUNAY TRIANGULATION OF REFERENCE ELEMENTS
        TrianglesFunc = Delaunay(GaussLobattoPoints)
        Triangles = TrianglesFunc.simplices.copy()

        # GET EQUALLY-SPACED/GAUSS-LOBATTO POINTS FOR THE EDGES
        GaussLobattoPointsOneD = GaussLobattoQuadrature(C+2)[0].flatten()
        if EquallySpacedPoints:
            GaussLobattoPointsOneD = EquallySpacedPointsnD(ndim-1, C).flatten()

        BasesQuad = np.zeros((nsize_2,GaussLobattoPoints.shape[0]),dtype=np.float64)
        hpBases = Quad.LagrangeGaussLobatto
        for i in range(GaussLobattoPoints.shape[0]):
            BasesQuad[:,i] = hpBases(CActual,GaussLobattoPoints[i,0],GaussLobattoPoints[i,1])[:,0]

        BasesOneD = np.zeros((CActual+2,GaussLobattoPointsOneD.shape[0]),dtype=np.float64)
        for i in range(GaussLobattoPointsOneD.shape[0]):
            BasesOneD[:,i] = LagrangeGaussLobatto(CActual,GaussLobattoPointsOneD[i])[0]

        # GET ONLY THE FACES WHICH NEED TO BE PLOTTED
        if ProjectionFlags is None:
            faces_to_plot_flag = np.ones(mesh.faces.shape[0])
        else:
            faces_to_plot_flag = ProjectionFlags.flatten()

        # CHECK IF ALL FACES NEED TO BE PLOTTED OR ONLY BOUNDARY FACES
        if faces_to_plot_flag.shape[0] > mesh.faces.shape[0]:
            # ALL FACES
            corr_faces = mesh.all_faces
            # FOR MAPPING DATA E.G. SCALED JACOBIAN FROM ELEMENTS TO FACES
            face_elements = mesh.GetElementsFaceNumberingHex()

        elif faces_to_plot_flag.shape[0] == mesh.faces.shape[0]:
            # ONLY BOUNDARY FACES
            corr_faces = mesh.faces
            # FOR MAPPING DATA E.G. SCALED JACOBIAN FROM ELEMENTS TO FACES
            face_elements = mesh.GetElementsWithBoundaryFacesHex()
        else:
            # raise ValueError("I do not understand what you want to plot")
            corr_faces = mesh.all_faces
            face_elements = mesh.GetElementsFaceNumberingHex()

        faces_to_plot = corr_faces[faces_to_plot_flag.flatten()==1,:]

        if QuantityToPlot is not None and plot_on_faces:
            quantity_to_plot = QuantityToPlot[face_elements[faces_to_plot_flag.flatten()==1,0]]

        # BUILD MESH OF SURFACE
        smesh = Mesh()
        smesh.element_type = "quad"
        # smesh.elements = np.copy(corr_faces)
        smesh.elements = np.copy(faces_to_plot)
        smesh.nelem = smesh.elements.shape[0]
        smesh.points = mesh.points[np.unique(smesh.elements),:]


        # MAP TO ORIGIN
        unique_elements, inv = np.unique(smesh.elements,return_inverse=True)
        mapper = np.arange(unique_elements.shape[0])
        smesh.elements = mapper[inv].reshape(smesh.elements.shape)

        smesh.GetBoundaryEdgesQuad()
        smesh.GetEdgesQuad()
        edge_elements = smesh.GetElementsEdgeNumberingQuad()



        # GET EDGE ORDERING IN THE REFERENCE ELEMENT
        reference_edges = NodeArrangementQuad(CActual)[0]
        reference_edges = np.concatenate((reference_edges,reference_edges[:,1,None]),axis=1)
        reference_edges = np.delete(reference_edges,1,1)

        # GET EULERIAN GEOMETRY
        if TotalDisp.ndim == 3:
            vpoints = mesh.points + TotalDisp[:,:ndim,-1]
        elif TotalDisp.ndim == 2:
            vpoints = mesh.points + TotalDisp[:,:ndim]
        else:
            raise AssertionError("mesh points and displacment arrays are incompatible")

        # svpoints = vpoints[np.unique(mesh.faces),:]
        svpoints = vpoints[np.unique(faces_to_plot),:]
        del vpoints
        gc.collect()

        if plot_edges:
            # GET X, Y & Z OF CURVED EDGES
            x_edges = np.zeros((C+2,smesh.all_edges.shape[0]))
            y_edges = np.zeros((C+2,smesh.all_edges.shape[0]))
            z_edges = np.zeros((C+2,smesh.all_edges.shape[0]))

            for iedge in range(smesh.all_edges.shape[0]):
                ielem = edge_elements[iedge,0]
                edge = smesh.elements[ielem,reference_edges[edge_elements[iedge,1],:]]
                coord_edge = svpoints[edge,:]
                x_edges[:,iedge], y_edges[:,iedge], z_edges[:,iedge] = np.dot(coord_edge.T,BasesOneD)


            # PLOT CURVED EDGES
            connections_elements = np.arange(x_edges.size).reshape(x_edges.shape[1],x_edges.shape[0])
            connections = np.zeros((x_edges.size,2),dtype=np.int64)
            for i in range(connections_elements.shape[0]):
                connections[i*(x_edges.shape[0]-1):(i+1)*(x_edges.shape[0]-1),0] = connections_elements[i,:-1]
                connections[i*(x_edges.shape[0]-1):(i+1)*(x_edges.shape[0]-1),1] = connections_elements[i,1:]
            connections = connections[:(i+1)*(x_edges.shape[0]-1),:]

        # CURVED SURFACES
        if plot_surfaces:

            nface = smesh.elements.shape[0]
            nnode = nsize*nface
            nelem = Triangles.shape[0]*nface

            Xplot = np.zeros((nnode,3),dtype=np.float64)
            Tplot = np.zeros((nelem,3),dtype=np.int64)

            # FOR CURVED ELEMENTS
            for ielem in range(nface):
                Xplot[ielem*nsize:(ielem+1)*nsize,:] = np.dot(BasesQuad.T, svpoints[smesh.elements[ielem,:],:])
                Tplot[ielem*TrianglesFunc.nsimplex:(ielem+1)*TrianglesFunc.nsimplex,:] = Triangles + ielem*nsize

            if QuantityToPlot is not None:
                Uplot = np.zeros(nnode,dtype=np.float64)
                if plot_on_faces:
                    # for ielem in range(nface):
                    #     Uplot[ielem*nsize:(ielem+1)*nsize] = quantity_to_plot[ielem]
                    # NEW APPROACH FOR CELL DATA - CHECK
                    Uplot = np.zeros(nelem,dtype=np.float64)
                    for ielem in range(nface):
                        Uplot[ielem*TrianglesFunc.nsimplex:(ielem+1)*TrianglesFunc.nsimplex] = quantity_to_plot[ielem]
                else:
                    # IF QUANTITY IS DEFINED ON NODES
                    quantity = QuantityToPlot[np.unique(faces_to_plot)]
                    for ielem in range(nface):
                        Uplot[ielem*nsize:(ielem+1)*nsize] = np.dot(BasesQuad.T, quantity[smesh.elements[ielem,:]])



        # THIS IS NOT A FLORENCE MESH COMPLIANT MESH
        tmesh = Mesh()
        tmesh.element_type = "tri"
        if plot_surfaces:
            tmesh.elements = Tplot
            tmesh.points = Xplot
            if QuantityToPlot is not None:
                tmesh.quantity = Uplot
            tmesh.nelem = nelem
            tmesh.nnode = nnode
            tmesh.nface = nface
        tmesh.nsize = nsize
        tmesh.bases_1 = BasesOneD
        tmesh.bases_2 = BasesQuad.T

        tmesh.smesh = smesh
        tmesh.faces_to_plot = faces_to_plot
        tmesh.svpoints = svpoints

        if plot_edges:
            tmesh.x_edges = x_edges
            tmesh.y_edges = y_edges
            tmesh.z_edges = z_edges
            tmesh.connections = connections
            tmesh.edge_elements = edge_elements
            tmesh.reference_edges = reference_edges

        return tmesh
